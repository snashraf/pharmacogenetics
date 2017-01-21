#!/usr/bin/python
# -*- coding: utf-8 -*-

# --------------------------------------------------------------------------

class Interpreter:

	def __init__(self, p):

		self.advice = {}

		self.conn = p.conn

		self.sql = p.conn.cursor()

		self.conn.commit()


	def Genotyper(self):
		
		# For each drug:
			
		self.sql.execute("SELECT DISTINCT DrugID FROM Pairs")
		
		for (did,) in self.sql.fetchall():
			
			# Select guidelines involved
			
			self.sql.execute("SELECT DISTINCT GuID from Guidelines WHERE DrugID = ?", (did,))
			
			for (guid,) in self.sql.fetchall():
		
				# fetch involved genes
		
				self.sql.execute("SELECT DISTINCT GeneID, o.GeneSymbol FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
					
				# Fetch involved haplotypes and scores
			
				genotype = {}
		
				for (gid, genesymbol) in self.sql.fetchall():
			
					self.sql.execute("SELECT DISTINCT p.HapID, Distance1, Distance2,h.hgvs, starname FROM PatHaplotypes p JOIN Haplotypes h on h.HapID = P.HapID WHERE h.GeneID = ?", (gid,))
					
					haplotypes = self.sql.fetchall()
					
					self.sql.execute("SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
					
					options = [starname[0] for starname in self.sql.fetchall()]
										
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
					
					print genotype[genesymbol]

					formatted_genotype = []
					
					for i, (gene, allele) in enumerate(genotype.items()):
						
						alleles = []
		
						for subdict in allele.values():
							
							refVal = ref["al{}".format(i+1)]

							print float(subdict['distance']), float(refVal)

							if float(subdict["distance"]) >= float(refVal):

								alleles.append(ref["starname"].replace("1A", "1"))

							else:	
								
								alleles.append(subdict["starname"])
						
						allele_string = "/".join(alleles)
						
						formatted_genotype.append(":".join([gene, allele_string]))
						
					string_genotype = ";".join(formatted_genotype)
						
					print guid
					
					print string_genotype
						
					raw_input("Print enter to continue...")
								
					# Find matching advice
						
					# Save to advice table 'PatGuidelines' (DrugID, GeneID, Category(Metabolizer type), Advice)


	def AnnFinder(self):
		
		pass
		
		# For each drug:
			
			#  Find involved genes
			
			# FOR GENE IN LIST:
			
				# Is there an identified guideline for the patient?
			
				# IF YES:
			
					# continue with next gene
			
				# IF NO:
			
					# can i find a haplotype for this gene?
			
					# if yes:
			
						# which closest haplotype has an annotation available?
			
					# if no:
			
						# continue on
			
			# go through 'variant'-tagged annotations
			
			# Find involved annotations and match with patientvars
			
			# DO SNPS and OTHERS seperately! GA/G becomes A/- for example, we need this for the annotation
			
			# Look up annotation
			
			# Get result (jinja?) and add to PatAnnotations

	def GetSNPAnnotations(self):

			# only when there is no haplotype available?
			
			self.authobj = Authenticate()
			
			self.remakeTable("patannotations")

			self.sql.execute('''
							SELECT v.VarID, p.CallBase, a.AnID
							FROM LocPGKB l
							JOIN PatientVars p
							ON p.Start = l.Start
							JOIN Variants v on l.VarId = v.VarId
							JOIN Annotations a
							ON a.VarHapId = v.VarId
							''')
							
			print "Annotating SNPs... /(* ` ^ `*/)"
			
			for (VarID, CallBase, AnID) in tqdm(self.sql.fetchall()):
				
				uri = \
				'https://api.pharmgkb.org/v1/data/clinicalAnnotation/{}?view=max' \
				.format(AnID)

				data = getJson(uri, self.authobj)
				
				allele = CallBase.replace("/", "")
				
				sql = self.insertSQL("patannotations").render(json = data, patallele = allele)
			
				self.sql.executescript(sql)

				# FOR SOMETHING THAT CREATES OVERVIEWS (new class? GUI? webserv?)
				
				overviewQuery = '''
				select distinct d.chemname, g.genename, p.phenotype, p.patallele, loe
				from patannotations p
				join annotations a
				on p.anid = a.anid
				join pairs p
				on a.drugid = p.drugid
				join drugs d
				on p.drugid = d.drugid
				join genes g
				on g.geneid = p.geneid
				order by g.genename;
				'''

				queryHaplotypeScores = '''
				select geneid, hapid, hapname, distance1, distance2 from pathaplotypes p
				join haplotypes h
				on p.HapID = h.hapid
				join genes g
				on h.geneid = g.geneid;
				'''
			

# ---------------------------- NEXT STEP: ReportMaker --------------------------------
