#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
from variant import Variant
from gene import Gene
from drug import Drug
from db import Database
from pair import Pair
from pgkb_functions import Authenticate, hg19conv
from tqdm import tqdm
import sqlite3

# ---------------------------------------------------------------------

class DataCollector(Database):

	'''
	This class imports a design VCF and creates a database with these positions, using PharmGKB.
	'''

	def __init__(self, dbname):
		"""
		takes database object to work on!
		"""
		Database.__init__(self, dbname)

		self.authobj = Authenticate()


	def GetPairs(self):
		"""
		Create table for drug-gene pairs, to be used later to fetch drug data
		:return:
		"""

		self.remakeTable('drugpairs')

		print 'Getting gene-drug pairs... ~( * v*)/\\(^w ^)~'

		# get the uri for finding well-annotated pairs (given by pharmgkb)

		uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

			# get data and read in this json file

		data = urllib2.urlopen(uri)

		guid = 'nan'

		options = 'nan'

		for doc in tqdm(json.load(data)):

			gid = doc['gene']['id']

			symbol = doc['gene']['symbol']

			did = doc['chemical']['id']

			p = Pair(gid, symbol, did, self.authobj)

			item = (p.gid, p.did, p.guid, p.options, ",".join(p.varids))

			self.insertValues("drugpairs", item)

		self.conn.commit()


	def GetGeneData(self):
		'''
		Fetches data on given gene IDs.
		:return:
		'''

		print 'Getting gene data... (/o*)'

		self.remakeTable("genes")

		self.remakeTable("haplotypes")

		# get all unique gene ids from the variant table

		self.sql.execute('SELECT DISTINCT gid FROM drugpairs')

		genes = self.sql.fetchall()

		# go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):
			
			g = Gene(gid)

			# insert the resulting name and haplotypes into sql table

			item = (gid, g.name, g.chr, g.start,
							 g.stop)

			self.insertValues("genes", item)

			for hap in g.haplotypes:

				for (rsid, alt) in hap['rsids']:

					item = (
						hap['id'],
						gid,
						hap['starname'],
						hap['hgvs'],
						rsid,
						alt,
						)

					self.insertValues("haplotypes", item)

			self.conn.commit()


	def GetVarData(self):
		"""
		Using Variant objects, this function fetches data on variants from the PharmGKB servers,
		if not available uses the Entrez servers, possibilities are limited for now.
		:return:
		"""

		print 'Getting variant data... ~(^_^)~'

		self.remakeTable("variants")

		self.remakeTable("alias")

		# get all rsids in the design vcf

		self.sql.execute('SELECT DISTINCT a.rsid, a.gid FROM haplotypes a JOIN genes g on a.gid = g.gid where rsid LIKE "rs%" order by a.gid, a.rsid'
						 )

		# rotate through rsids and create variant objects to fetch information

		for (rsid, gid) in tqdm(self.sql.fetchall()):

			# create variant instances with the given rsid.
			
			try:

				v = Variant(rsid)

			except KeyError:

				continue

			# this results in a combination tuple of rsid and gid and aliases

			item = (
				rsid,
				v.id,
				gid,
				v.chr,
				v.begin,
				v.end,
				v.ref,
				v.alt,
				v.nbegin,
				v.nend,
				v.nref,
				v.nalt,
				v.muttype,
				)

			self.sql.execute('''INSERT INTO variants VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)'''
							 , item)

			# go through aliases, ignore duplicates and put in alias table

			for alias in v.names:

				try:

					self.insertValues("alias", (rsid, alias))
				
				except sqlite3.IntegrityError:

					# on duplicate, ignore

					continue

			self.conn.commit()


	def GetNonRS(self):

		self.sql.execute("SELECT rsid, alt, gid from haplotypes where rsid LIKE '%chr%' order by rsid"
						 )

		print 'Parsing non-rs variants...'

		for (rsid, alt, gid) in tqdm(self.sql.fetchall()):

			item = hg19conv(rsid, gid, alt)

			if item is not None:

				self.insertValues("variants", item)

		self.conn.commit()


	def GetOthers(self):
		'''
		Lalalala
		'''
		pass


	def GetDrugData(self):
		"""
		Gets info on seperate drugs, such as name and associated terms
		terms give info on what the drug does
		:return:
		"""

		print 'Getting drug data /(>_<)\\}'

		# drop and re-create table drugs

		self.remakeTable("drugs")

		# get all the important drug ids (dids) from
		# the known gene-drug connections table

		self.sql.execute('SELECT DISTINCT did FROM drugpairs')

		# fetch matching did and use to create uri for query

		for result in tqdm(self.sql.fetchall()):

			did = result[0]

			d = Drug(did)

			for item in d.terms:

				term = item['term']

				# check for duplicates with chemical name

				if d.name not in term.lower():

					item = (did, d.name, term)

					# insert into drugs table

					self.insertValues("drugs", item)

		self.conn.commit()


	def BedFile(self):

		# creates bed file for subsetting .BAM files.

		self.sql.execute('SELECT chr, start, stop, symbol FROM genes ORDER BY length(chr), chr'
						 )

		with open('PharmacogenomicGenes_PGKB.bed', 'w') as bed:

			for tup in self.sql.fetchall():

				bed.write('\t'.join(map(str, tup)) + '\n')

# -----------------------------------------------------------------