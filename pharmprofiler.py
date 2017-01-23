#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser

from lib.data import *
from lib.patient import *

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

	{"short":"r" , "long":"reset", "action":"store","dest":"reset", "default":None,
			"help":"Database table to reset (mostly for testing purposes)"}
]

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

		print options.reset
		CreateDB(options.dbname, options.dbtasks, options.reset)

	if options.gvcf:

		if '.gz' not in options.gvcf:

			print 'Please convert to .gz and create tabix file first.'

		else:

			CreatePatient(options.dbname, options.gvcf, options.pattasks)


def CreateDB(dbname, tables, reset):

# --------------------------------------------------------------------------

	d = DataCollector(dbname)

	print reset

	if reset is not None:

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

		("guideoptions", d.GetGuideOptions)

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

	("interpret", p.Interpret)

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

# --------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
