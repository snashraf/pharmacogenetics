#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser

from lib.data import *
from lib.patient import *
from lib.reportmaker import *

# --------------------------------------------------------------------------

def main():

	options = [

	{"short":"n" , "long":"name", "action":"store","dest":"dbname", "default":"database",
			"help":"Specify database name. This will be used to create/access a database. Default is 'database'"},

	{"short":"b" , "long":"build", "action":"append","dest":"dbtasks", "default":[],
			"help":"Download database (specify tables, default is do whole database.\
			Options are pairs (gene-drug pairs), genes, vars (allvars/rsvars/hapvars) \
			and drugs(chemicals"},

	{"short":"a" , "long":"action", "action":"append","dest":"pattasks", "default":[],
		"help":"Which actions to perform on patient file. Current options: 'import' \
		(import important snp sites from vcf), 'haplotype' (phylogenetically haplotype patient),\
		'annotate' (annotate found changes and identified patient haplotypes)."},

	{"short":"p" , "long":"patient", "action":"store","dest":"gvcf", "default":None,
			"help":"Patient compressed vcf [g.vcf.gz] file to parse"},

	{"short":"o" , "long":"outfile", "action":"store","dest":"outfile", "default":".",
			"help":"Location to send report TEX file to."},

	{"short":"r" , "long":"reset", "action":"append","dest":"reset", "default":[],
			"help":"Database table to reset (mostly for testing purposes)"},

	{"short":"i" , "long":"interpret", "action":"append","dest":"interpret", "default":[],
			"help":'''Do something with collected data. Use 'guidelines' to find patient matching guidelines and add to DB.
			Use 'annotate' to do the same for annotations. Use 'report' after this to generate a patient report using this info.'''}]

	parser = OptionParser(usage='usage: %prog [options] filename',
						version='%prog 1.0')

	for o in options:

		parser.add_option(
			'-{}'.format(o['short']),
			'--{}'.format(o['long']),
			action=o['action'],
			dest=o['dest'],
			default=o['default'],
			help=o['help']
			)

	(options, args) = parser.parse_args()

	if len(options.dbtasks) > 0:
		CreateDB(options.dbname, options.dbtasks, options.reset)
	if options.gvcf:
		CreatePatient(options.dbname, options.gvcf, options.pattasks)
	if len(options.interpret) > 0:
		InterpretResults(options.dbname, options.interpret, options.outfile)

def CreateDB(dbname, tables, reset):

# --------------------------------------------------------------------------

	d = DataCollector(dbname)

	print reset

	if len(reset) > 0:
		for table in reset:
			d.remakeTable(reset)

# --------------------------------------------------------------------------

	options = OrderedDict ([

		("pairs", d.GetPairs),

		("drugs", d.GetDrugs),

        ("drugmatches", d.GetDrugMatches),

		("drugvars", d.GetDrugVars),

		("genes", d.GetGeneData),

		("haplotypes", d.GetHaplotypes),

		("hapvars", d.GetHapVars),

		("etcvars", d.GetNonRS),

		("indels", d.ConvertIndels),

		("annotations", d.GetAnnotations),

		("guidelines", d.GetGuidelines),

		("guideoptions", d.GetGuideOptions),

		("bed", d.BedFile)

		])

	options['vars'] = [options['drugvars'], options["hapvars"], options["etcvars"]]

	options['all'] = options.values()

# --------------------------------------------------------------------------

	# Consider putting options in DICT
	# Consider list comprehension
	# [dict[table] for table in tables]

	for table in tables:

		try:

			o = options[table]

			if not hasattr( d, "authobj"):

				d.Authenticate()

			if type(o) is not list:

				options[table]()

			elif type(o) is list:

				for item in o:

					item()

		except:

			raise

			print "Invalid option entered. \n Valid options: {}".format(", ".join(options.keys()))

	d.conn.commit()


def CreatePatient(dbname, gvcf, tables):

	p = Patient(dbname, gvcf)

	options = OrderedDict ([

	("import", p.GetPositions),

	("haplotype", p.GetHaplotypes),

	("score", p.HapScorer)

	])

	options['all'] = options.values()

# ----------------------------------------------

	for table in tables:

		try:

			o = options[table]

			if type(o) is not list:

				options[table]()

			elif type(o) is list:

				for item in o:

					item()

		except:

			raise

			print "Invalid option entered. \n Valid options: {}".format(", ".join(options.keys()))

	p.conn.commit()

def InterpretResults(dbname, opts, outfile):

	i = Interpreter(dbname)
	r = ReportMaker(dbname, outfile)

# ===========================

	options = OrderedDict ([

	("haplotype", i.Haplotype),

	("guidelines", i.FindGuidelines),

	("annotate", i.Annotate),

	("report", [r.MakeJson, r.MakeReport])

	])

	for opt in opts:

		try:

			o = options[opt]

			if type(o) is not list:

				options[opt]()

			elif type(o) is list:

				for item in o:

					item()

		except:
			raise
			print "Invalid option entered. \n Valid options: {}".format(", ".join(options.keys()))
			return


# --------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
