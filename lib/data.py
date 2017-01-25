#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import OrderedDict
from modules.pgkb_functions import Authenticate, hg19conv, getJson, getRef
import urllib2
import json
from tqdm import tqdm
from db import Database
import pprint as pp
import os
import sys
import csv
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
		self.path = os.path.dirname(__file__)

		if len(self.path) == 0:

			self.path = os.getcwd()

		dbfolder = os.path.join(self.path, 'db')

		self.dbpath = os.path.join(dbfolder, '%s.db' % dbname)

		Database.__init__(self, self.dbpath)


	def Authenticate(self):

		self.authobj = Authenticate()


	def GetPairs(self):
		"""
		Create table for drug-gene pairs, to be used later to fetch drug data
		:return:
		"""

		self.remakeTable('pairs')

		print 'Getting gene-drug pairs... ~( * v*)/\\(^w ^)~'

		# get the uri for finding well-annotated pairs (given by pharmgkb)

		uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

		data = getJson(uri, self.authobj)

		template = self.insertSQL("pairs")

		for response in tqdm(data):

			sql = template.render(json = response)

			self.sql.executescript(sql)

			self.conn.commit()


	def GetDrugs(self):

		self.remakeTable("drugs")

		self.drugs = []

		csv.field_size_limit(sys.maxint)

		with open(self.path+"/drugs/drugs.tsv", "rb") as f:

				tsv = csv.reader(f, delimiter='\t')

				for row in tsv:

					did = row[0]

					name = row[1]

					smiles = row[6]

					self.sql.execute("INSERT INTO DRUGS VALUES(?,?,?)",(did, name, smiles))

		self.conn.commit()


	def GetDrugMatches(self):

		print "Getting variants connected to drugs... (- w-)b"

		#self.remakeTable("drugvars")

		self.sql.execute('''
			SELECT DISTINCT DrugID FROM Drugs
			WHERE DrugID NOT IN
			(SELECT DrugID FROM	DrugVars)
		''')

		for (did,) in tqdm(self.sql.fetchall()):

			uri = 'https://api.pharmgkb.org/v1/report/connectedObjects/{}/Variant'\
				.format(did)

			data = getJson(uri, self.authobj)

			# INSERT INTO HAPVARS WITH N/A HAPLOTYPE FOR NOW

			sql = self.insertSQL("drugvars").render(json = data, did=did)

			if "INSERT" in sql:
				print sql

			self.sql.executescript(sql)

			self.conn.commit()


	def GetDrugVars(self):

		template = self.insertSQL("variants")

		self.sql.execute('''
					SELECT DISTINCT VarID from DrugVars
					WHERE VarID NOT IN (
					SELECT DISTINCT VarID from Variants
					);'''
					)

		# rotate through rsids and create variant objects to fetch information

		print "Getting more info on these variants..."

		uris = []

		for (varid,) in tqdm(self.sql.fetchall()):

			# create variant instances with the given rsid.

			uri = \
				'https://api.pharmgkb.org/v1/data/variant/{}?&view=max'.format(varid)

			uris.append(uri)

			data = getJson(uri, self.authobj)

			sql = template.render(json = data)

			self.sql.executescript(sql)

			self.conn.commit()


	def GetGeneData(self):
		'''
		Fetches data on given gene IDs.
		:return:
		'''

		print 'Getting gene data... (/o*)'

		# get all unique gene ids from the variant table

		self.sql.execute('''SELECT DISTINCT GeneID FROM Variants
		WHERE GeneID NOT IN (
			SELECT GeneID FROM Genes)''')

		for (gid,) in tqdm(self.sql.fetchall()):

			uri = 'https://api.pharmgkb.org/v1/data/gene/{}?view=max'.format(gid)

			data = urllib2.urlopen(uri)

			response = json.load(data)
			try:
				sql = self.insertSQL("genes").render(json = response)
			except:
				pp.pprint(response)
				continue

			self.sql.executescript(sql)

			self.conn.commit()


	def GetHaplotypes(self):

		self.sql.execute('''
		SELECT DISTINCT GeneID FROM Variants
		WHERE GeneID NOT IN (
		SELECT GeneID FROM Haplotypes)''')

		template = self.insertSQL("haplotypes")

		# TODO CATCH TABLE DOES NOT EXIST

		genes = self.sql.fetchall()

		# go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):

			uri = 'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max'.format(gid)

			data = getJson(uri, self.authobj)

			if not data:

				continue

			for response in data:

				sql = template.render(json = response)

				print sql

				self.sql.executescript(sql)

			self.conn.commit()


	def GetHapVars(self):
		"""
		:return:
		"""

		print 'Getting variants connected to haplotypes... ~(^_^)~'

		self.sql.execute('''
					SELECT DISTINCT VarName from HapVars
					WHERE VarName LIKE "rs%" AND VarName NOT IN
					(SELECT DISTINCT RSID FROM Variants)
					''')

		template = self.insertSQL("variants")

		# rotate through rsids and create variant objects to fetch information

		for (rsid,) in tqdm(self.sql.fetchall()):

			# create variant instances with the given rsid.

			uri = \
				'https://api.pharmgkb.org/v1/data/variant/?symbol={}&view=max'.format(rsid)

			data = getJson(uri, self.authobj)

			sql = template.render(json = data[0])

			self.sql.executescript(sql)

			self.conn.commit()


	def GetNonRS(self):

		self.sql.execute('''
				SELECT DISTINCT v.VarName, v.AltAllele, h.GeneID from HapVars v
				JOIN Haplotypes h ON v.HapID = h.HapID
				JOIN Genes g on h.GeneID = g.GeneID
				WHERE v.VarName LIKE "%chr%"
				AND v.VarName NOT LIKE "hg38"
				ORDER BY h.GeneID;
				''')

		print 'Parsing non-rs variants...'

		template = self.insertSQL("othervars")

		for (rsid, gid, alt) in tqdm(self.sql.fetchall()):

			d = hg19conv(rsid, gid, alt)

			sql = template.render(json = d)

			print sql

			self.sql.executescript(sql)

		self.conn.commit()


	def ConvertIndels(self):

		template = self.insertSQL("indels")

		self.remakeTable("indels")

		self.sql.execute('''
						SELECT DISTINCT
						l.VarID, RefGenome, Chromosome, Start, End, RefAllele, AltPGKB
						FROM LocPGKB l
						JOIN Variants v
						ON l.VarID = v.VarID
						JOIN AltAlleles a
						ON v.VarID = a.VarID
						''')

		print "Converting indels.. \(>w <)/"

		for (varid, genome, loc, start, end, ref, alt) in tqdm(self.sql.fetchall()):

			# create json for template usage

			shifted = {"varid":varid,
					"chromosome":loc,
					"genome":genome,
					"ref":ref,
					"alt":alt,
					"start":start,
					"end":end}

			# insertion or deletion?

			if alt == ref:

				continue

			if ref == "-":

				# left shift position by 1

				shifted['start'] = start - 1

				shifted['end'] = end

				# get reference nucleotide at that position

				prevbase = getRef(loc, shifted['start'], shifted['start'])

				# insertion scenario

				shifted['ref'] = prevbase

				shifted['alt'] = prevbase + alt

			elif alt == "-":

				# left shift position by 1

				shifted['start'] = start - 1

				shifted['end'] = end

				# get reference nucleotide at that position

				prevbase = getRef(loc, shifted['start'], shifted['start'])

				# deletion scenario A -

				shifted['ref'] = prevbase + ref

				salt = prevbase

				shifted['alt'] = salt

			else:

				pass

			# render sql

			sql = template.render(alt = alt, json = shifted)

			self.sql.executescript(sql)

			self.conn.commit()

# ------------------------------------------------------------------------------------------------------------

	def GetAnnotations(self):

		self.sql.execute('SELECT DISTINCT DrugID FROM Drugs \
									WHERE DrugID NOT IN (SELECT DISTINCT DrugID FROM Annotations)')

		for (DrugID,) in tqdm(self.sql.fetchall()):

			uri = \
			'https://api.pharmgkb.org/v1/data/clinicalAnnotation?relatedChemicals.accessionId={}' \
			.format(DrugID)

			data = getJson(uri, self.authobj)

			sql = self.insertSQL("annotations").render(json = data, DrugID = DrugID)

			print sql

			self.sql.executescript(sql)

			self.conn.commit()


	def GetGuidelines(self):

		template = self.insertSQL("guidelines")

		#self.remakeTable("guidelines")

		self.sql.execute('SELECT DISTINCT DrugID FROM Drugs')

		for (DrugID,) in tqdm(self.sql.fetchall()):

			uri = 'https://api.pharmgkb.org/v1/data/guideline?relatedChemicals.accessionId={}'\
					.format(DrugID,)

			data = getJson(uri, self.authobj)

			pp.pprint(data)

			if data is None:

				continue

			sql = template.render(json = data, did = DrugID)

			print sql

			self.sql.executescript(sql)

		self.conn.commit()


	def GetGuideOptions(self):

		template = self.insertSQL("guideoptions")

		self.remakeTable("guideoptions")

		self.sql.execute("SELECT DISTINCT GuID from Guidelines;")

		for (guid, ) in tqdm(self.sql.fetchall()):

			uri = "https://api.pharmgkb.org/v1/report/guideline/{}/options" \
			.format(guid)

			data = getJson(uri, self.authobj)

			if data is None:

				continue

			sql = template.render(guid = guid, json = data)

			self.sql.executescript(sql)

		self.conn.commit()


	def BedFile(self):

		# creates bed file for subsetting .BAM files.

		self.sql.execute('SELECT chromosome, start, end, genesymbol FROM genes ORDER BY chromosome ASC, start ASC'
						 )

		linesINT = []
		linesETC = []

		for (chrom, start, end, name) in self.sql.fetchall():
			chrom = chrom.lstrip("chr")
			try:
				int(chrom)
				linesINT.append('\t'.join(map(str, (chrom, start, end, name))) + '\n')
			except:
				linesETC.append('\t'.join(map(str, (chrom, start, end, name))) + '\n')

		linesINT.sort(key=lambda x: int(x.split("\t")[0]))

		linesALL = linesINT + linesETC

		with open('PharmacogenomicGenes_PGKB_full.bed', 'w') as bed:
			bed.writelines(linesALL)

# -----------------------------------------------------------------
