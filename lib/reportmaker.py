#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from jinja2 import Template, FileSystemLoader, Environment

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

TEMPLATE = '''
\documentclass{resume} % Use the custom resume.cls style

\usepackage[left=0.75in,top=0.6in,right=0.75in,bottom=0.6in]{geometry} % Document margins

\name{{{SampleName}}} % Your name
\address{Generated on \today} % Your address

\begin{document}

{% for json in jsonlist %}

\begin{rSection}{{{json.drugName}}}
% \item {{json.drugDesc}}

{% for gene in json.relatedGenes %}

\begin{rSubsection}{{{gene.geneName}}}{gene.haplotype}{gene.patGuideline.metaType}{}
% \item \textit{{{gene.description}}}
\newline\newline ------------------------------------------------------ Dosing Guideline --------------------------------------------------------
\item \textbf{{{gene.haplotype}}} \newline {{gene.patGuideline.recommendation}}
\newline\newline ---------------------------------------------------- Clinical Annotations -----------------------------------------------------

{% for ann in json.patAnnotations %}
\item \textbf{{{ann.varName}}} {{ann.patAllele}} \newline {{ann.annotation}}
{% endfor %}

\end{rSubsection}

{% endfor %}

\end{rSection}

{% endfor %}

\end{document}
'''
# ==========================================================

# CREATE A VIEW FOR NECESSARY DATA

# List of drugs (sorted alphabetically, maybe?)

sql.execute("SELECT DISTINCT DrugID, DrugName FROM Variants")

jsons = []

for (did, name) in self.sql.fetchall():

    js={}
    js['drugID'] = did
    js['drugName'] = name
    js['relatedGenes']=[]

# For each drug:
# --------- HEADER: Drug Name -----------------

# Collect involved genes (through bound variants in DrugVars)

    sql.execute('''
                        SELECT DISTINCT v.GeneID, g.Symbol, pg.genotype
                        FROM DrugVars d
                        JOIN Variants v ON d.VarID = v.VarID
                        JOIN Genes g ON v.GeneID = g.GeneID
                        JOIN PatGenotypes pg ON pg.GeneID = g.GeneID
                        ''')

# For each gene:

    for (gid, symbol, haplotype) in self.sql.fetchall():

        js_gene = {}
        js_gene['geneID'] = gid
        js_gene['geneName'] = symbol
        js_gene['haplotype'] = haplotype

        # Get involved Guideline if it exists

        js_guide = js_gene['patGuide']
        js_guide = {}

        sql.execute('''
                            SELECT DISTINCT GuID, name, Phenotype, Strength, Recommendation
                            FROM PatGuidelines g
                            ''')

        for (guid, name, phen, strength, rec) in self.sql.fetchall():

            js_guide['metaType'] = name
            js_guide['phenotype'] = phen
            js_guide['recommendation'] = rec
            js_guide['strength'] = strength

# Print Annotations for this gene-drug combination

        sql.execute...
        for (guid, name, phen, strength, rec) in self.sql.fetchall():

            js_guide['metaType'] = name
            js_guide['phenotype'] = phen
            js_guide['recommendation'] = rec
            js_guide['strength'] = strength

        jsons.append[js]
