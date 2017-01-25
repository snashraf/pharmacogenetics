#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from jinja2 import Template, FileSystemLoader, Environment
import os
import pprint as pp

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
'''
# ============================================================

class ReportMaker(object):
	def __init__(self):
		print "Initiating Report..."
		self.path = os.path.dirname(__file__)

		if len(self.path) == 0:
			self.path = os.getcwd()

		self.conn = sqlite3.connect(self.path+"/db/pgx_full.db")  # connect to db- if it doesn't exist, create it

		self.sql = self.conn.cursor()  # cursor for self.sqlite3, used to do things in database

		self.templateText = '''
		\documentclass{resume} % Use the custom resume.cls style

		\usepackage[left=0.75in,top=0.6in,right=0.75in,bottom=0.6in]{geometry} % Document margins

		\name{\\{{SampleName}}\\} % Your name
		\address{Generated on \today} % Your address

		\begin{document}

		{% for json in jsonlist %}

		\begin{rSection}{\\{{json.drugName}}\\}
		% \item DrugDesc

		{% for gene in json.relatedGenes %}
		\begin{rSubsection}{\\{{gene.geneName}}\\}{TEST}{TEST}{TEST}
		\newline\newline ------------------------------------------------------ Dosing Guideline --------------------------------------------------------
		{% if gene.patGuide %}
		\item \textbf{\\{{gene.haplotype}}\\}  \newline{{gene.patGuideline.phenotype}} \newline {{gene.patGuideline.recommendation}}
		{% endif %}
		{% if gene.patAnnotations %}
		\newline\newline ---------------------------------------------------- Clinical Annotations -----------------------------------------------------
		{% for ann in json.patAnnotations %}
		\item \textbf{\\{{ann.varName}}\\} {{ann.patAllele}} \newline {{ann.annotation}}
		{% endfor %}
		{% endif %}
		\end{rSubsection}

		{% endfor %}

		\end{rSection}

		{% endfor %}

		\end{document}
		'''
	# ==========================================================

		# AnnotationView

		self.sql.executescript('''DROP VIEW AnnOverview;
						CREATE VIEW AnnOverview AS
				        SELECT DISTINCT v.GeneID, d.DrugID, g.GeneSymbol, v.RSID, p.PatAllele, p.Phenotype
		                FROM DrugVars d
		                JOIN Variants v ON d.VarID = v.VarID
		                JOIN Genes g ON v.GeneID = g.GeneID
                        JOIN Annotations a
                    	ON a.VarHapID = d.VarID
                        JOIN PatAnnotations p
                        ON p.AnID = a.AnID;''')

		self.template = Template(self.templateText)

		# Hap/GuidelineView

		self.sql.executescript('''DROP VIEW HapOverview;
						CREATE VIEW HapOverview AS
						SELECT DISTINCT * FROM PatGuidelines p
						JOIN Guidelines g
						ON g.GuID = p.GuID;''')

		# List of drugs (sorted alphabetically, maybe?)

		self.sql.execute('''SELECT DISTINCT d.DrugID, d.ChemName
							FROM DrugVars v
							JOIN Drugs d
							ON d.DrugID = v.DrugID
							ORDER BY d.ChemName ASC''')
		jsons = []

		for (did, name) in self.sql.fetchall():
			print did, name
			js={}
			js['drugID'] = did
			js['drugName'] = name
			js['relatedGenes']=[]

		# Collect involved genes (through bound variants in DrugVars)

			self.sql.execute('''
		                        SELECT DISTINCT GeneID, GeneSymbol
		                        FROM AnnOverview
								WHERE DrugID = ?''', (did,))

		# For each gene:

			for (gid, symbol) in self.sql.fetchall():
				js_gene = {}
				js_gene['geneID'] = gid
				js_gene['geneName'] = symbol
				js_gene['patAnnotations'] = []
				js_gene['patGuide'] = {}
				js_guide = js_gene['patGuide']
				js_guide = {}

				self.sql.execute('''
				        SELECT DISTINCT  MetaCat, Strength, Term, Markdown FROM HapOverview
						WHERE DrugID = ? AND GeneID = ?''', (did, gid))

				for (metacat, strength, term, markdown) in self.sql.fetchall():
					js_guide['metaType'] = metacat
					js_guide['strength'] = strength
					if "Implication" in term:
						js_guide['implications'] = markdown
					if "Phenotype" in term:
						js_guide['phenotype'] = term
					if "Recommendation" in term:
						js_guide['recommendation'] = markdown

		# Print Annotations for this gene-drug combination

				self.sql.execute('''
									SELECT DISTINCT RSID, PatAllele, phenotype
									FROM AnnOverview
									WHERE DrugID = ? AND GeneID = ?
									''', (did, gid))

				for (rsid, allele, phenotype) in self.sql.fetchall():
					item = {"varName":rsid, "patAllele":allele, "annotation":phenotype}
					js_gene['patAnnotations'].append(item)
				js['relatedGenes'].append(js_gene)
			jsons.append(js)
		pp.pprint(jsons)
		reportText = self.template.render(sampleName="TEST", jsonlist=jsons)
		print reportText
# ----------------------------------------------------------------------------

r = ReportMaker()
