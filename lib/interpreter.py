#!/usr/bin/python
# -*- coding: utf-8 -*-

from modules.pgkb_functions import *
from tqdm import tqdm
import pprint as pp
# --------------------------------------------------------------------------

class Interpreter:

	def __init__(self, p):

		self.p = p

		self.advice = {}

		self.authobj = Authenticate()

		self.p.conn.commit()


	def Genotyper(self):
		
		self.p.remakeTable("patguidelines")
		# For each drug:
			
		self.p.sql.execute("SELECT DISTINCT DrugID FROM Pairs")
		
		for (did,) in self.p.sql.fetchall():
			
			# Select guidelines involved
			
			self.p.sql.execute("SELECT DISTINCT GuID from Guidelines WHERE DrugID = ?", (did,))
			
			for (guid,) in self.p.sql.fetchall():
		
				# fetch involved genes
		
				self.p.sql.execute("SELECT DISTINCT GeneID, o.GeneSymbol FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
					
				# Fetch involved haplotypes and scores
			
				genotype = {}
		
				for (gid, genesymbol) in self.p.sql.fetchall():
			
					self.p.sql.execute("SELECT DISTINCT p.HapID, Distance1, Distance2,h.hgvs, starname FROM PatHaplotypes p JOIN Haplotypes h on h.HapID = P.HapID WHERE h.GeneID = ?", (gid,))
					
					haplotypes = self.p.sql.fetchall()
					
					self.p.sql.execute("SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
					
					options = [starname[0] for starname in self.p.sql.fetchall()]
										
					if len(haplotypes) == 0:
						
						continue
							
					genotype[genesymbol] = {"al1":

							{"starname":"A", "distance":10.0},
								"al2":
							{"starname":"A", "distance":10.0}
							}
						
					for (hapid, al1, al2, hgvs, starname) in haplotypes:

						if "[=]" in hgvs:

							ref = {"starname":starname, "al1":al1, "al2":al2}

						# FOR ALLELE IN PATIENT:
						
						for i, alleleScore in enumerate([float(al1), float(al2)]):
							
							curLoc = genotype[genesymbol]["al{}".format(i+1)]
							
							curScore = curLoc["distance"]
		
							if alleleScore < curScore and starname in options:
									
								curLoc["starname"] = starname
			
								curLoc["distance"] = alleleScore

					# -----------------------------------------------------
					
					formatted_genotype = []
					
					for i, (gene, allele) in enumerate(genotype.items()):
						
						alleles = []
		
						for subdict in allele.values():
							
							refVal = ref["al{}".format(i+1)]

							if float(subdict["distance"]) >= float(refVal):

								alleles.append(ref["starname"])

							else:	
								
								alleles.append(subdict["starname"])
						
						allele_string = "/".join(alleles)
						
						formatted_genotype.append(":".join([gene, allele_string]))
						
					string_genotype = ";".join(formatted_genotype)
														
					# Find matching advice

					uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(guid, string_genotype) 
						
					data = getJson(uri, self.authobj)
		
					if data == None:

						uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(guid, string_genotype.replace("1A", "1"))

						data = getJson(uri, self.authobj)

					sql = self.p.insertSQL("patguidelines").render(guid = guid, genotype = string_genotype, json = data)
					
					self.p.sql.executescript(sql)

				self.p.conn.commit()

				# Save to advice table 'PatGuidelines' (DrugID, GeneID, Category(Metabolizer type), Advice)


	def Annotate(self):

			# only when there is no haplotype available?
						
			self.p.remakeTable("patannotations")

			self.p.sql.execute('''
							SELECT DISTINCT       
        							a.AnID,							
        							VarName,
        							VarID,
            						MutType,
									RefPGKB,
									AltPGKB,
									RefVCF,
									AltVCF,
									PatRef,
									PatAlt,
									CallNum,
									Start
							FROM Overview
							JOIN Annotations a
							on a.VarHapID = VarID
							WHERE RefVCF = PatRef''')

							
			print "Annotating SNPs... /(* ` ^ `*/)"
			
			for (anid, rsid, varid, muttype, refP, altP, refV, altV, ref, alt, call, start) \
			in tqdm(self.p.sql.fetchall()):
				
				uri = \
				'https://api.pharmgkb.org/v1/data/clinicalAnnotation/{}?view=max' \
				.format(anid)

				data = getJson(uri, self.authobj)
				
				# --------------------------------------

				if muttype == "snp":
					
					call = call.replace("/", "")
					allele = call.replace("0", ref)
					allele = allele.replace("1", alt)
					revAllele = allele[::-1]

				elif muttype != "snp":

					allele = call.replace("0", refP.replace("-", "del"))
					allele = allele.replace("1", altP.replace("-", "del"))
					s_allele = allele.split("/")
					revAllele = s_allele[1] + "/" + s_allele[0]

				# --------------------------------------
				
				sql = self.p.insertSQL("patannotations").render(json = data, revallele = revAllele, patallele = allele)
			
				self.p.sql.executescript(sql)

			self.p.conn.commit()
			

# ---------------------------- NEXT STEP: ReportMaker --------------------------------
