#!/usr/bin/python
# -*- coding: utf-8 -*-

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

			# get data and read in this json file

		data = urllib2.urlopen(uri)

		guid = 'nan'

		options = 'nan'

		for doc in tqdm(json.load(data)):

			sql = self.insertSQL("pairs").render(json = doc)

			self.sql.executescript(sql)

		self.conn.commit()


	def GetAnnotations(self):
		pass


	def GetGuidelines(self):
		pass


	def GetGeneData(self):
		'''
		Fetches data on given gene IDs.
		:return:
		'''

		print 'Getting gene data... (/o*)'

		self.remakeTable("genes")

		# get all unique gene ids from the variant table

		self.sql.execute('SELECT DISTINCT gid FROM drugpairs')
		
		# TODO CATCH TABLE DOES NOT EXIST

		genes = self.sql.fetchall()

		# go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):

			g = Gene(gid)

			sql = self.insertSQL("genes").render(json = g.json)

			self.sql.executescript(sql)

		self.conn.commit()


	def GetHaplotypes(self):

		self.remakeTable("haplotypes")

		self.sql.execute('SELECT DISTINCT gid FROM genes')
		
		# TODO CATCH TABLE DOES NOT EXIST

		genes = self.sql.fetchall()

		# go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

		for (gid,) in tqdm(genes):

			uri = 'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId={}&view=max'.format(self.gid)

        	data = urllib2.urlopen(uri)

        	response = json.load(data)

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

		self.sql.execute("SELECT DISTINCT h.RSID FROM Haplotypes h \
						JOIN Genes g on h.gid = g.gid where rsid LIKE \
						 "rs%" order by a.gid, a.rsid'
						 )"

		# rotate through rsids and create variant objects to fetch information

		for (rsid,) in tqdm(self.sql.fetchall()):

			# create variant instances with the given rsid.

			v = Variant(rsid)

        	response = json.load(data)

			sql = self.insertSQL("variants").render(json = response)

			self.sql.executescript(sql)

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

		self.sql.execute('SELECT DISTINCT did FROM pairs')

		# fetch matching did and use to create uri for query

		for result in tqdm(self.sql.fetchall()):

			did = result[0]

			d = Drug(did)

			self.insertSQL("drugs", d.json)

			for term in d.json['terms']:

				text = item['term']
				
				item = (did, text)

				self.insertValues("drug_terms", item)

		self.conn.commit()


	def BedFile(self):

		# creates bed file for subsetting .BAM files.

		self.sql.execute('SELECT chr, start, stop, symbol FROM genes ORDER BY length(chr), chr'
						 )

		with open('PharmacogenomicGenes_PGKB.bed', 'w') as bed:

			for tup in self.sql.fetchall():

				bed.write('\t'.join(map(str, tup)) + '\n')

# -----------------------------------------------------------------
