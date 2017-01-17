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

# ---------------------------------------------------------------------

class DataCollector(Database):

	'''
	This class imports a design VCF and creates a database with these positions, using PharmGKB.
	'''

	def __init__(self, dbname):
		"""
		takes database object to work on!
		"""
		path = os.path.dirname(__file__)

		if len(path) == 0:

			path = os.getcwd()

		dbfolder = os.path.join(path, 'db')

		dbpath = os.path.join(dbfolder, '%s.db' % dbname)

		Database.__init__(self, dbpath)


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

		for response in tqdm(data):

			sql = self.insertSQL("pairs").render(json = response)
			
			self.sql.executescript(sql)

			self.conn.commit()


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

		self.sql.execute('SELECT DISTINCT DrugID FROM Pairs')

		# fetch matching did and use to create uri for query

		for (did,) in tqdm(self.sql.fetchall()):

			uri = \
			'https://api.pharmgkb.org/v1/data/chemical/{}?view=max'.format(did)

			data = getJson(uri, self.authobj)
			
			sql = self.insertSQL("drugs").render(json = data)
			
			self.sql.executescript(sql)


		self.conn.commit()


	def GetGeneData(self):
		'''
		Fetches data on given gene IDs.
		:return:
		'''

		print 'Getting gene data... (/o*)'

		self.remakeTable("genes")

		# get all unique gene ids from the variant table

		self.sql.execute('SELECT DISTINCT GeneID FROM Pairs')
		
		# TODO CATCH TABLE DOES NOT EXIST

		genes = self.sql.fetchall()

		# go through results and creat e gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):
			
			uri = 'https://api.pharmgkb.org/v1/data/gene/{}?view=max'.format(gid)

			data = urllib2.urlopen(uri)

			response = json.load(data)

			sql = self.insertSQL("genes").render(json = response)

			self.sql.executescript(sql)

			self.conn.commit()


	def GetHaplotypes(self):

		self.remakeTable("haplotypes")

		self.sql.execute('SELECT DISTINCT GeneID FROM Genes')
		
		# TODO CATCH TABLE DOES NOT EXIST

		genes = self.sql.fetchall()

		# go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):

			uri = 'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max'.format(gid)

			data = getJson(uri, self.authobj)
			
			if not data:
				
				continue
				
			for response in data:

				sql = self.insertSQL("haplotypes").render(json = response)

				self.sql.executescript(sql)

		self.conn.commit()


	def GetVarData(self):
		"""
		Using Variant objects, this function fetches data on variants from the PharmGKB servers,
		if not available uses the Entrez servers, possibilities are limited for now.
		:return:
		"""

		print 'Getting variant data... ~(^_^)~'

		self.remakeTable("variants")

		# get all rsids in the design vcf

		self.sql.execute('''
					SELECT DISTINCT v.VarName from HapVars v
					JOIN Haplotypes h ON v.HapID = h.HapID
					JOIN Genes g on h.GeneID = g.GeneID
					WHERE v.VarName LIKE "rs%"
					ORDER BY h.GeneID;
					'''
					)

		# rotate through rsids and create variant objects to fetch information

		for (rsid,) in tqdm(self.sql.fetchall()):

			# create variant instances with the given rsid.

			uri = \
				'https://api.pharmgkb.org/v1/data/variant/?symbol={}&view=max'.format(rsid)

			data = getJson(uri, self.authobj)

			# pp.pprint(response)

			sql = self.insertSQL("variants").render(json = data[0])

			self.sql.executescript(sql)

		self.conn.commit()


	def GetNonRS(self):

		self.remakeTable("othervars")

		self.sql.execute('''
				SELECT DISTINCT v.VarName, v.AltAllele, h.GeneID from HapVars v
				JOIN Haplotypes h ON v.HapID = h.HapID
				JOIN Genes g on h.GeneID = g.GeneID
				WHERE v.VarName LIKE "%chr%"
				ORDER BY h.GeneID;
				''')

		print 'Parsing non-rs variants...'

		for (rsid, gid, alt) in tqdm(self.sql.fetchall()):

			d = hg19conv(rsid, gid, alt)

			sql = self.insertSQL("othervars").render(json = d)

			self.sql.executescript(sql)

		self.conn.commit()


	def ConvertIndels(self):

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

			sql = self.insertSQL("indels").render(alt = alt, json = shifted)

			print sql

			self.sql.executescript(sql)


# ------------------------------------------------------------------------------------------------------------


	def GetAnnotations(self):
		
		self.remakeTable("annotations")

		self.sql.execute('SELECT DISTINCT DrugID, GeneID FROM Pairs')

		for (DrugID, GeneID) in tqdm(self.sql.fetchall()):
			
			uri = \
			'https://api.pharmgkb.org/v1/report/pair/{}/{}/clinicalAnnotation?view=max' \
			.format(DrugID, GeneID)

			data = getJson(uri, self.authobj)
						
			sql = self.insertSQL("annotations").render(json = data, DrugID = DrugID)
		
			self.sql.executescript(sql)

		self.conn.commit()



	def GetGuidelines(self):
		
		self.sql.execute('SELECT DISTINCT DrugID, GeneID FROM Pairs')

		for (DrugID, GeneID) in self.sql.fetchall():
			
			sources = ['cpic', 'dpwg', 'pro']
			
			for source in sources:
				
				uri = \
					'https://api.pharmgkb.org/v1/data/guideline?source={}&relatedChemicals.accessionId={}&relatedGenes.accessionId={}&view=max' \
					.format(source, DrugID, GeneID)
	
				data = getJson(uri, authobj)
	
				pp.pprint(data)
				
				sql = self.insertSQL("guidelines").render(json = d)
	
				print sql
	
				self.sql.executescript(sql)

		self.conn.commit()


	def BedFile(self):

		# creates bed file for subsetting .BAM files.

		self.sql.execute('SELECT chr, start, stop, symbol FROM genes ORDER BY length(chr), chr'
						 )

		with open('PharmacogenomicGenes_PGKB.bed', 'w') as bed:

			for tup in self.sql.fetchall():

				bed.write('\t'.join(map(str, tup)) + '\n')

# -----------------------------------------------------------------
