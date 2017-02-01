#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from jinja2 import Template, FileSystemLoader, Environment
import os
import pprint as pp
from db import Database
from patient import Patient
from modules.pgkb_functions import DeconstructGuideline

# Needed for document:

JSON_TEMPLATE ='''
# --------------------------------------------------------------------------------------

{samplename, haplotable:[ {symbol:{ new1:{1, 2, 3}, new2{1, 2, 3}, old1{1, 2, 3}, old2{1, 2, 3}, tree)}, ...],
						annotable:[ (drug, gene, guidelinelink, Ann1-2, Ann3-4), ...],
						drugs:[drugname:{drugDesc, genes:[{symbol, hapNEW, hapOLD, guideline, annotations:[:{1A,
						 	1B:{patAllele, rsid, text},
							2A, 2B, 3, 4}, ...]
						]}
							}
						} , ...]

}
'''
# ============================================================

class ReportMaker(Database):

	def __init__(self, dbname, outfile):
		print "Initiating Report..."
		self.path = os.path.dirname(__file__)
		dbfolder = os.path.join(self.path, 'db')
		dbpath = os.path.join(dbfolder, '%s.db' % dbname)
		self.outfile = outfile
		Database.__init__(self, dbpath)

	# ==========================================================

		# AnnotationView

		#self.template = self.getTemplate("latex/python_template.tex")

	def MakeJson(self):
		self.sn = raw_input("Please enter a sample name: ")

		# gather data for tables
		JSON = {"sampleName":self.sn,
					  "haplotable":[],
					  "annotable":[],
					  "pairs":[]
					  }

		self.sql.execute('''
					select g.genesymbol, g.geneid,
					new1_1, new1_2, new1_3,
					new2_1, new2_2, new2_3,
					old1_1, old1_2, old1_3,
					old2_1, old2_2, old2_3
					from
					patgenotypes p
					join genes g
					on g.GeneID = p.GeneID
								''')

		for (symbol, gid, n11, n12, n13, n21, n22, n23, o11, o12, o13, o21, o22, o23) in self.sql.fetchall():
				entry = {"symbol":symbol,
				 "new1":{"1":n11, "2":n12, "3":n13},
				 "new2":{"1":n21, "2":n22, "3":n23},
				 "old1":{"1":o11, "2":o12, "3":o13},
				 "old2":{"1":o21, "2":o22, "3":o23}}
				tree = ""
				with open(self.path + "/output/alignments/aligned/" + gid + "_aln_rooted.dnd", "rb") as f:
					self.sql.execute('SELECT DISTINCT hapid, starname FROM Haplotypes WHERE GeneID = ?', (gid,))
					hapconvert = {hapid:starname for (hapid, starname) in self.sql.fetchall()}
					content = f.readlines()
					for line in content:
						for hapid, starname in hapconvert.items():
							starname = starname.replace(",", "")
							line = line.replace(hapid, starname).replace("a1", "{}_hap1".format(self.sn)).replace("a2", "{}_hap2".format(self.sn))
						tree += line.replace(" ", "").rstrip("\n")
					entry['tree'] = tree
					JSON['haplotable'].append(entry)

				# ------------------------------

		# get gene-drug pairs

		self.sql.executescript('''DROP VIEW IF EXISTS JsonView;
								CREATE VIEW JsonView AS
								select distinct
								g.genesymbol as symbol,
								g.geneid as geneid,
								g.genename as genename,
								d.chemname as drugname,
								d.drugid as drugid,
								l.guid as guid,
								l.markdown as markdown,
								l.summary as summary,
								a.loe as loe,
								v.rsid as rsid,
								p.patallele as allele,
								p.phenotype as phenotype
								from genes g
								join variants v
								on v.geneid = g.geneid
								join annotations a
								on a.varhapid = v.varid
								join drugs d
								on a.drugid = d.drugid
								join patannotations p
								on p.anid = a.anid
								left join guidelines l
								on l.geneid = g.geneid
								and l.drugid = d.drugid
								order by d.chemname asc;
									''')

		self.sql.execute('''
		SELECT DISTINCT symbol, geneid, genename, drugname, drugid, guid from JsonView
		''')
		for (symbol, gid, genename, drugname, did, guid) in self.sql.fetchall():
			# get annotation amounts
			self.sql.execute('''
					select distinct LoE, count(*) from JsonView
					where geneid = ?
					and drugid = ?
					group by LoE;
					''', (gid, did))

			for (loe, amount) in self.sql.fetchall():
				cat1 = 0
				cat2 = 0
				if "1" in loe or "2" in loe:
					cat1 += amount
				elif "2" in loe or "3" in loe:
					cat2 += amount
				else:
					pass

			# entry for summary table 2

			annoEntry = {"drug":drugname,
								  "gene":symbol,
								  "guid":guid,
								  "annCount":
								  		{"1-2":cat1, "3-4":cat2}
								 }

			mainEntry = {}
			JSON['annotable'].append(annoEntry)

			# ------ get more data on guidelines and annotations -----
			geneEntry = {"drug":drugname,
			"genesymbol":symbol,
			"genename":genename,
			"guideline":
				{"guid":"", "summary":"", "tex":""},
			"annotations":{"1A":[], "1B":[], "2A":[], "2B":[], "3":[], "4":[]}}

			for hap in JSON['haplotable']:
				if hap["symbol"] != symbol:
					continue
				else:
					geneEntry['hapNEW'] = "{}/{}".format(hap['new1']["1"], hap['new2']["1"])
					geneEntry['hapOLD'] = "{}/{}".format(hap['old1']["1"], hap['old2']["1"])

			# ^ this goes in json[drugs][DRUGNAME][genes] through append to list ^

			# --- get the annotation data ---
			self.sql.execute('''
						select distinct
						loe, allele,
						rsid, phenotype,
						guid, markdown, summary
						from JsonView
						where geneid = ?
						and drugid = ?
						''', (gid, did))

			for (loe, allele, rsid, phenotype, guid, markdown, summary) in self.sql.fetchall():
				# get guideline info
				if guid != None:
					geneEntry['guideline']['guid']=guid
					geneEntry['guideline']['summary']=summary
					geneEntry['guideline']['tex']=DeconstructGuideline(markdown)
				# get annotation info
				annEntry = {"rsid":rsid, "patAllele":allele, "phenotype":phenotype}
				if annEntry not in geneEntry['annotations'][loe]:
					geneEntry['annotations'][loe].append(annEntry)

			JSON['pairs'].append(geneEntry)

		pp.pprint(JSON['pairs'])

		self.JSON = JSON
# -----------------------------------------------------------------------------------------

	def MakeReport(self):
		reportText = self.template.render(sampleName=self.sn, jsonlist=self.JSON)
		with open(self.outfile + self.sn + ".tex", "wb") as f:
			f.write(reportText.encode('utf-8'))

# ----------------------------------------------------------------------------
