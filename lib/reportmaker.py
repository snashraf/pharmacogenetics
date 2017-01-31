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

JSON_TEMPLATE = '''
{"drugID":DrugID, "drugName":DrugName,relatedGenes:
	[
		{"geneID":geneID, "description":desc, "geneName":geneName, "patGuideline":
												{"metaType":metaType, "phenotype":phenotype,
												"recommendation":reccomendation},
										"patAnnotations":
													[
										{"varID":varID, "patAllele":patAllele, "annotation":annotation},
										{} etc.
													]
	]
}

# --------------------------------------------------------------------------------------

{samplename, haplotable:[ {symbol, new1:{1, 2, 3}, new2{1, 2, 3}, old1{1, 2, 3}, old2{1, 2, 3}, tree}, ...],
						annotable:[ (drug, gene, guidelinelink, Ann1-2, Ann3-4), ...],
						drugs:[{drugname, genes:[{symbol, hapNEW, hapOLD, guideline, annotations:[:{1A,
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

		self.sql.executescript('''DROP VIEW IF EXISTS AnnOverview;
						CREATE VIEW AnnOverview AS
				        SELECT DISTINCT v.GeneID, g.GeneName, a.DrugID, g.GeneSymbol, v.RSID, a.LoE, p.PatAllele, p.Phenotype
		                FROM Annotations a
		                JOIN Variants v ON a.VarHapID = v.VarID
		                JOIN Genes g ON v.GeneID = g.GeneID
                        JOIN PatAnnotations p
                        ON p.AnID = a.AnID;''')

		self.template = self.getTemplate("latex/template_python.tex")

		# Hap/GuidelineView

		self.sql.executescript('''DROP VIEW IF EXISTS HapOverview;
						CREATE VIEW HapOverview AS
						SELECT DISTINCT
						GuID,
						g.GeneID,
						DrugID,
						Allele1_Phyl as new1,
						Allele2_Phyl as new2,
						Allele1_Set as old1,
						Allele2_Set as old2,
						Markdown
						FROM Guidelines g
						JOIN PatGenotypes p
						on p.GeneID = g.GeneID;''')

		# List of drugs (sorted alphabetically, maybe?)

	def MakeJson(self):

		colorchart = {
								"1A":"red",
								"1B":"orange",
								"2A":"cyan",
								"2B":"blue",
								"3":"teal",
								"4":"green"
								}

		self.sql.execute('''SELECT DISTINCT a.DrugID, d.ChemName
							FROM AnnOverview a
							JOIN Drugs d
							ON d.DrugID = a.DrugID
							ORDER BY d.ChemName ASC''')

		self.jsons = []

		for (did, name) in self.sql.fetchall():
			print did, name
			js={}
			js['drugID'] = did
			js['drugName'] = name

		# Collect involved genes (through bound variants in DrugVars)

			self.sql.execute('''
		                        SELECT DISTINCT GeneID, GeneSymbol, GeneName
		                        FROM AnnOverview
								WHERE DrugID = ?''', (did,))

		# For each gene:
			print name
			for (gid, symbol, name) in self.sql.fetchall():
				js_gene = {}
				js_gene['geneID'] = gid
				js_gene['geneName'] = symbol
				js_gene['geneDesc'] = name
				print name

				self.sql.execute('''
				        SELECT DISTINCT new1, new2, old1, old2, Markdown FROM HapOverview
						WHERE DrugID = ? AND GeneID = ?''', (did, gid))

				for (n1, n2, o1, o2, markdown) in self.sql.fetchall():
					js_guide = js_gene['patGuide'] = {}
					js_guide['hapNEW'] = "{}/{}".format(n1, n2)
					js_guide['hapOLD'] =  "{}/{}".format(o1, o2)
					tex = DeconstructGuideline(markdown)
					print tex
					js_guide['tex'] = tex

		# Print Annotations for this gene-drug combination

				self.sql.execute('''
									SELECT DISTINCT RSID, PatAllele, phenotype, LoE
									FROM AnnOverview
									WHERE DrugID = ? AND GeneID = ?
									ORDER BY LoE ASC
									''', (did, gid))

				for (rsid, allele, phenotype, loe) in self.sql.fetchall():
					item = {"varName":rsid, "patAllele":allele, "annotation":phenotype, "strength":loe, "color":colorchart[loe]}
					js_gene.setdefault("patAnnotations", []).append(item)
				js.setdefault("relatedGenes", []).append(js_gene)

			self.jsons.append(js)

	def MakeReport(self):
		sn = raw_input("Please enter a sample name: ")

		reportText = self.template.render(sampleName=sn, jsonlist=self.jsons)
		with open(self.outfile + sn + ".tex", "wb") as f:
			f.write(reportText.encode('utf-8'))

# ----------------------------------------------------------------------------
