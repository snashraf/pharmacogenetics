#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser

from lib.data import DataCollector

from lib.db import Database

from lib.patient import Patient


def main():

	parser = OptionParser(usage='usage: %prog [options] filename',
						  version='%prog 1.0')

	parser.add_option(
		'-n',
		'--name',
		action='store',
		dest='dbname',
		default="database",
		help='''
		Specify database name. This will be used to create/access a database. Default is "database".
		''' ,
		)

	parser.add_option(
		'-t',
		'--maketable',
		action='append',
		dest='tables',
		default=[],
		help='''
	Download database (specify tables, default is do whole database)
	Options are pairs (gene-drug pairs), genes, vars (allvars/rsvars/hapvars) and drugs(chemicals)'
		''' ,
		)

	parser.add_option(  # optional because action defaults to "store"
		'-p',
		'--patient',
		action='store',
		dest='gvcf',
		default=None,
		help='Patient compressed vcf [g.vcf.gz] file to parse',
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

	d = DataCollector(dbname)

	if 'all' in tables:

		d.Update()

	# Consider putting options in DICT
	# Consider list comprehension
	# [dict[table] for table in tables]
	else:

		for table in tables:

			if table == 'pairs':

				d.GetPairs()

			elif table == 'genes':

				d.GetGeneData()

			elif table == 'allvars':

				d.GetVarData()

				d.GetNonRS()

			elif table == 'rsvars':

				d.GetVarData()

			elif table == 'drugs':

				d.GetDrugData()

			elif table == 'hapvars':

				d.GetNonRS()

		d.conn.commit()


def CreatePatient(dbname, gvcf):

	p = Patient(dbname, gvcf)

# --------------------------------------------------------------------------------

if __name__ == '__main__':

	main()
