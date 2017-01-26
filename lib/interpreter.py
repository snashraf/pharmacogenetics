#!/usr/bin/python
# -*- coding: utf-8 -*-

from modules.pgkb_functions import *
from tqdm import tqdm
import pprint as pp
from db import Database
import os

# --------------------------------------------------------------------------

class Interpreter(Database):

	def __init__(self, dbname):
		self.path = os.path.dirname(__file__)
		dbfolder = os.path.join(self.path, 'db')
		dbpath = os.path.join(dbfolder, '%s.db' % dbname)

		Database.__init__(self, dbpath)

		self.advice = {}

		self.authobj = Authenticate()

		self.conn.commit()


	def Haplotype(self):

		self.sql.executescript('''
						DROP TABLE IF EXISTS PatGenotypes;
						CREATE TABLE IF NOT EXISTS PatGenotypes
						(GeneID text,
						Allele1_Phyl text, Allele2_Phyl text,
						Allele1_Set text, Allele2_Set text);
		''')

		self.sql.execute('''
		SELECT DISTINCT h.GeneID, g.GeneSymbol from Haplotypes h
		JOIN Genes g on g.GeneID = h.GeneID''')

		for (gid, genesymbol) in self.sql.fetchall():

			genotype={}

			print genesymbol

			self.sql.execute("SELECT DISTINCT p.HapID, Distance1, Distance2, OldScore1, OldScore2, HapLen, h.hgvs, starname FROM PatHaplotypes p JOIN Haplotypes h on h.HapID = P.HapID WHERE h.GeneID = ?", (gid,))

			haplotypes = self.sql.fetchall()

			if len(haplotypes) == 0:

				continue

			genotype[genesymbol] = {
					"al1_PHEN":
					{"starname":"A", "distance":10.0},
					"al1_OLD":
					{"starname":"A", "score":0.0, "hapLen":0},
					"al2_PHEN":
					{"starname":"A", "distance":10.0},
					"al2_OLD":
					{"starname":"A", "score":0.0, "hapLen":0}
					}

			for (hapid, al1, al2, s1, s2, hLen, hgvs, starname) in haplotypes:

				if "[=]" in hgvs:

					ref = {"starname":starname, "al1":al1, "al2":al2}

				# FOR ALLELE IN PATIENT:

				for i, alleleDistance in enumerate([float(al1), float(al2)]):

					curLoc = genotype[genesymbol]["al{}_PHEN".format(i+1)]

					curDistance = curLoc["distance"]

					if alleleDistance < curDistance:

						curLoc["starname"] = starname

						curLoc["distance"] = alleleDistance

		# ---------------------------------- calculate old style --------------------------------------

				for i, alleleScore in enumerate([float(s1), float(s2)]):

					curLoc = genotype[genesymbol]["al{}_OLD".format(i+1)]

					curScore = curLoc["score"]

					curLen = curLoc["hapLen"]

					if (alleleScore >= curScore and hLen > curLen) and alleleScore > 1.0:

						curLoc["starname"] = starname

						curLoc["score"] = alleleScore

						curLoc["hapLen"] = hLen

			# -----------------------------------------------------

			genotypeNEW = []

			genotypeOLD = []

			for i, (gene, allele) in enumerate(genotype.items()):

				allelesNEW = []

				allelesOLD = []

				for name, subdict in allele.items():

					if "PHEN" in name:

						refVal = ref["al{}".format(i+1)]

						if float(subdict["distance"]) >= float(refVal):

							allelesNEW.append(ref["starname"])

						else:

							allelesNEW.append(subdict["starname"])

					if "OLD" in name:

						if float(subdict["score"]) > 1.0:

							allelesOLD.append(subdict["starname"])

						else:

							allelesOLD.append(ref["starname"])

				item = (gid, allelesNEW[0], allelesNEW[1], allelesOLD[0], allelesOLD[1])

				self.sql.execute("INSERT INTO PatGenotypes VALUES(?,?,?,?,?)", item)
				self.conn.commit()


	def FindGuidelines(self):

		self.remakeTable("patguidelines")
		# For each drug:

		self.sql.execute("SELECT DISTINCT GuID from Guidelines")

		for (guid,) in tqdm(self.sql.fetchall()):
			print guid
			# fetch involved genes

			self.sql.execute("SELECT DISTINCT o.GeneID, g.GeneSymbol FROM Guidelines o JOIN Genes g on g.GeneID = o.GeneID WHERE GuID = ?", (guid,))

			# Fetch involved haplotypes and scores

			genotype = {}

			for (gid, genesymbol) in self.sql.fetchall():

				self.sql.execute("SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))

				options = [starname[0] for starname in self.sql.fetchall()]

				self.sql.execute('''
				SELECT DISTINCT
				Allele1_Phyl, Allele2_Phyl,
				Allele1_Set, Allele2_Set
				FROM PatGenotypes
				WHERE GeneID = ?
				''', (gid,))

				scores = self.sql.fetchone()

				genotypeNEW = [scores[0], scores[1]]
				genotypeOLD = [scores[2], scores[3]]

# ----------------------------------------------------------------------

			string_genotype_NEW = ";".join(genotypeNEW)
			string_genotype_OLD = ";".join(genotypeOLD)

			print "NEW:", string_genotype_NEW
			print "OLD:", string_genotype_OLD

			# Find matching advice
			string_genotype = string_genotype_OLD

			uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(guid, string_genotype)

			data = getJson(uri, self.authobj)

			if data == None:

				uri = "https://api.pharmgkb.org/v1/report/guideline/{}/annotations?genotypes={}".format(guid, string_genotype.replace("1A", "1"))

				data = getJson(uri, self.authobj)

			sql = self.insertSQL("patguidelines").render(guid = guid, genotype = string_genotype, json = data)
			print sql
			self.sql.executescript(sql)

			self.conn.commit()

			# Save to advice table 'PatGuidelines' (DrugID, GeneID, Category(Metabolizer type), Advice)


	def Annotate(self):

			# only when there is no haplotype available?

			self.sql.execute('''
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
							WHERE RefVCF = PatRef
							AND a.AnID NOT IN
							(SELECT DISTINCT AnID FROM PatAnnotations)''')

			print "Annotating SNPs... /(* ` ^ `*/)"

			for (anid, rsid, varid, muttype, refP, altP, refV, altV, ref, alt, call, start) \
			in tqdm(self.sql.fetchall()):

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

				sql = self.insertSQL("patannotations").render(json = data, revallele = revAllele, patallele = allele)

				try:
					self.sql.executescript(sql)
				except:
					print sql
					continue

			self.conn.commit()


# ---------------------------- NEXT STEP: ReportMaker --------------------------------
