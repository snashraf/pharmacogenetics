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
				
	{"short":"t" , "long":"table", "action":"append","dest":"tables", "default":[],
			"help":"Download database (specify tables, default is do whole database.\
			Options are pairs (gene-drug pairs), genes, vars (allvars/rsvars/hapvars) \
			and drugs(chemicals"},
				
	{"short":"p" , "long":"patient", "action":"store","dest":"gvcf", "default":None,
			"help":"Patient compressed vcf [g.vcf.gz] file to parse"}
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

	if len(options.tables) > 0:

		CreateDB(options.dbname, options.tables)

	if options.gvcf:

		if '.gz' not in options.gvcf:

			print 'Please convert to .gz and create tabix file first.'

		else:

			CreatePatient(options.dbname, options.gvcf)


def CreateDB(dbname, tables):

# --------------------------------------------------------------------------

	d = DataCollector(dbname)

# --------------------------------------------------------------------------

	options = OrderedDict ([
		
		("pairs", d.GetPairs),
		
		("genes", d.GetGeneData),

		("haplotypes", d.GetHaplotypes),
		
		("rsvars", d.GetVarData),
		
		("hapvars", d.GetNonRS),
		
		("drugs", d.GetDrugData),
		
		("annotations", d.GetAnnotations)
		
		])

	options['vars'] = [options["rsvars"], options["hapvars"]]

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


def CreatePatient(dbname, gvcf):

	p = Patient(dbname, gvcf)

# --------------------------------------------------------------------------------

if __name__ == '__main__':

	main()