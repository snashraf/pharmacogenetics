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
						New1_1 text, New1_2 text, New1_3 text,
						New2_1 text, New2_2 text, New2_3 text,
						Old1_1 text, Old1_2 text, Old1_3 text,
						Old2_1 text, Old2_2 text, Old2_3 text);
		''')

		self.sql.execute('''
		SELECT DISTINCT h.GeneID, g.GeneSymbol from Haplotypes h
		JOIN Genes g on g.GeneID = h.GeneID''')

		genotype={}

		for (gid, genesymbol) in self.sql.fetchall():

			self.sql.executescript('''
			DROP VIEW IF EXISTS CurView;
			CREATE VIEW CurView AS
			SELECT DISTINCT
			h.GeneID as GeneID,
			p.HapID,
			p.Distance1 as new1, p.Distance2 as new2,
			po.score1 as old1, po.score2 as old2,
			po.MatchLen,
			h.hgvs, starname
			FROM PatHaplotypes p
			JOIN PatHaplotypes_OLD po
			on p.hapid = po.hapid
			JOIN Haplotypes h
			on h.HapID = P.HapID
			''')

			haplotypes = self.sql.fetchall()

			self.sql.execute('''
	SELECT DISTINCT starname, new1, new2, old1, old2 FROM CurView
	WHERE GeneID = ?
	and hgvs LIKE "%=%"
				''',(gid,))

			refhap = {}

			for (starname, new1, new2, old1, old2) in self.sql.fetchall():
				refhap['starname'] = starname
				refhap['new1']=float(new1)
				refhap['new2']=float(new2)
				refhap['old1']=0.0
				refhap['old2']=0.0

			allelesNEW = {1:[], 2:[]}
			allelesOLD = {1:[], 2:[]}
			allelesNEW = {1:[], 2:[]}
			allelesOLD = {1:[], 2:[]}

			# NEW STYLE

			for i in range(1,3):
				self.sql.execute('''
				SELECT DISTINCT starname, new{} FROM CurView
				WHERE GeneID = ?
				ORDER BY new{} ASC
				'''.format(i, i), (gid,))

				# shuffle for ref at the top
				for (starname, score) in self.sql.fetchall():
					if float(score) == refhap["new{}".format(i)] and starname != refhap['starname'] and refhap['starname'] not in allelesNEW[i]:
						allelesNEW[i].append(refhap['starname'])
						allelesNEW[i].append(starname)
					else:
						allelesNEW[i].append(starname)

				# OLD STYLE
				self.sql.execute('''
				SELECT DISTINCT starname, old{}, matchLen FROM CurView
				WHERE GeneID = ?
				ORDER BY old{} DESC, MatchLen DESC
				'''.format(i, i), (gid,))

				for (starname, score, mlen) in self.sql.fetchall():
					if float(score) == 0.0 and starname != refhap['starname'] and refhap['starname'] not in allelesOLD[i]:
						allelesOLD[i].append(refhap['starname'])
						allelesOLD[i].append(starname)
					else:
						allelesOLD[i].append(starname)
			try:
				item = (gid, allelesNEW[1][0],allelesNEW[1][1],allelesNEW[1][2],
				allelesNEW[2][0],allelesNEW[2][1],allelesNEW[2][2],
				allelesOLD[1][0], allelesOLD[1][1], allelesOLD[1][2],
				allelesOLD[2][0], allelesOLD[2][1], allelesOLD[2][2])
				print item
				self.sql.execute("INSERT INTO PatGenotypes VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", item)
				self.conn.commit()
			except:
				continue

	def FindGuidelines(self):

		self.remakeTable("patguidelines")
		# For each drug:

		self.sql.execute("SELECT DISTINCT GuID from Guidelines")

		for (guid,) in tqdm(self.sql.fetchall()):

			print guid
			# fetch involved genes

			self.sql.execute("SELECT DISTINCT o.GeneID, g.GeneSymbol FROM Guidelines o JOIN Genes g on g.GeneID = o.GeneID WHERE GuID = ?", (guid,))

			genes = self.sql.fetchall()

			guideGenes = {}

			if len(genes) == 0:
				continue
			# Fetch involved haplotypes and scores

			genotype = {}

			for (gid, genesymbol) in genes:

				self.sql.execute("SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))

				options = [starname[0] for starname in self.sql.fetchall()]

				self.sql.execute('''
				SELECT DISTINCT
				New1_1, New2_1,
				Old1_1, Old2_1
				FROM PatGenotypes
				WHERE GeneID = ?
				''', (genesymbol,))

				scores = self.sql.fetchone()
				print scores

				if scores == None:
					continue

				guideGenes[gid] = {}
				guideGenes[gid]['new'] =  "{}:{}/{}".format(genesymbol, scores[0], scores[1])
				guideGenes[gid]['old'] =  "{}:{}/{}".format(genesymbol, scores[2], scores[3])

		# ----------------------------------------------------------------------

				string_genotype_NEW = ";".join([guideGenes[gid]['new'] for gid in guideGenes.keys()])
				string_genotype_OLD = ";".join([guideGenes[gid]['old'] for gid in guideGenes.keys()])

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
