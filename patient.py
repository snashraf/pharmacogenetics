#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import sqlite3
from collections import Counter
from operator import itemgetter

# ---------------------------------------------------------------------

class Patient:

	'''
	This class imports a design VCF and an input VCF.
	'''

	def __init__(self, f):
		"""
		f = given VCF input file
		GetDesign imports the test design and matches rs values to 
		positions
		Import imports the input vcf
		GetIDs finds changed positions and matches them to the design.
		"""
		print 'Initiating patient...'
		self.f = f
		self.reader = vcf.Reader(open(f, 'r'))
		self.genes = {}
		self.haplotypes = []
		self.starhaps = {}
		print "Loading patient data..."
		self.Load()

	def Load(self):
		"""
		Function for loading the most important startup functions
		:return:
		"""
		self.ImportData()
		self.GetIDs()
		self.conn.commit()
		#self.Hapmaker()
		#self.conn.commit()
		self.Hapmatcher()
		self.conn.commit()

	def ImportData(self):
		"""
		Import data from collected database and
		create patientvars tables.
		:return:
		"""

		print 'Trying to import database...'
		try:

			# connect to db

			self.conn = sqlite3.connect('pharmacogenetics.db')
			self.sql = self.conn.cursor()

			# drop and re-create patient tables
			# (maybe another db would be better?)

			self.sql.execute('''DROP TABLE IF EXISTS patientvars''')
			self.sql.execute('''DROP TABLE IF EXISTS patienthaps''')

			# patientvars focusus on positions and rsids,
			# patienthaps focuses on calculated haplotype storage

			self.sql.execute('''CREATE TABLE patientvars
												(chr text, 
												start int, 
												end int, 
												ref text, 
												alt text, 
												call text, 
												pid text, 
												pgt text)'''
							 )
			self.sql.execute('''CREATE TABLE patienthaps
												(gid text, 
												hapid text, 
												chr text,
												hapscore int,
												patscore int)'''
							 )

		except:

			# if error in db, reload db (testing purposes...)
			# self.d.Update()

			raise

	def GetIDs(self):
		"""
		Gets patient variables, reads into patientvars table
		:return:
		"""
		print 'Reading patient variants...'

		# create list of positions to localize in patient

		for record in self.reader:

			self.sql.execute("SELECT DISTINCT chr, start, end, ref, alt FROM variants")
			positions = self.sql.fetchall()
			for (chr, start, end, ref, alt) in positions:
				try:
					records = self.reader.fetch(str(chr), start=start-1, end=end)
				except:
					continue
				for record in records: # doctest: +SKIP
					ref = str(record.REF)
					alt = (str(record.ALT[0])).replace("<NON_REF>",".")
					for sample in record.samples:
						call = str(sample['GT'])  # 1/0 etc, phasing found here
						try:
							pid = str(sample['PID'])
							pgt = str(sample['PGT'])
						except:
							pid = "nan"
							pgt = "nan"
					item = (chr, start, end, ref, alt, call, pid, pgt)
					self.sql.execute('''INSERT INTO patientvars VALUES(?,?,?,?,?,?,?,?)'''
							, item)

	def Hapmaker(self):
		'''
		IN PROGRESS:
		Create HGVS haplotypes from rsids. Should implement phasing...
		:return:
		'''

		print "Starting hapmaker,..."

		# get gids that are present in patient by joining tables

		self.sql.execute("""SELECT DISTINCT v.gid from variants v
							JOIN patientvars p ON v.start = p.start""")

		# make an easy usuable list for later on, remove tuples

		gids = [tup[0] for tup in self.sql.fetchall()]

		# loop through gids

		for gid in gids:

			# get rsid, alias for that rsid, and bases at this position in chr1 and 2

			self.sql.execute('''SELECT v.rsid, a.alias, p.ref, p.alt, p.call
								FROM patientvars p
								JOIN variants v ON p.start = v.start
								JOIN alias a ON a.rsid = v.rsid
								WHERE v.gid = ?;''',
							 (gid, ))

			results = self.sql.fetchall()

			# create dictionary for storing genomic positions, amb is for unphased variants

			pos = {'amb': [], 'chr1': [], 'chr2': []}

			# find the version of the alias notation that is the highest by using counter and max

			counts = Counter(elem[1].split('.')[1] for elem in results
							 if 'NC' in elem[1] and 'g.' in elem[1])

			chosen_ver = str(max(counts.keys()))

			# loop through results

			for (rsid, alias, ref, alt, num) in results:
				# select on "NC" and genomic position

				if 'NC' in alias and 'g.' in alias:

					if '>' in alias and chosen_ver in alias:

						seq_id = alias.split('.')[0]

						gen_pos = alias.split('.')[2]

						if '|' in num:

							if '|1' in num:

								pos['chr2'].append(gen_pos)

							if '1|' in num:

								pos['chr1'].append(gen_pos)

						elif '/' in num:

							pos['amb'].append(gen_pos)

			# check for unchanged positions - make these "=" as hgvs reccs

			for (k, v) in pos.items():

				if len(v) == 0:

					pos[k] = ['=']

				else:

					# sort based on genomic position (split the genomic position on base change, leaving the number)

					v.sort(key=lambda x: x.split('ACTG')[0])

			# check if the haplotype notation should be the phased version or not

			if '/' in num:

				haplotype = '%s.%s.%s' % (seq_id, chosen_ver,
						';'.join(pos['amb']))

			elif '|' in num:

				haplotype = '%s%s.[%s];[%s]' % (seq_id, chosen_ver,
						';'.join(pos['chr1']), ';'.join(pos['chr2']))

			item = (gid, haplotype, 'N/A', 'N/A')

			self.sql.execute('''INSERT INTO patienthaps VALUES(?,?,?,?)'''
							 , item)

	def Hapmatcher(self):
		
		"""
		Matches rsids per gene to known haplotypes for that gene.
		Results are stored in a list...
		:return:
		"""

		self.sql.execute("DROP TABLE IF EXISTS patienthaps")

		self.sql.execute("CREATE TABLE IF NOT EXISTS patienthaps(gid text, hapid text, score1 int, score2 int)")

		self.sql.execute("SELECT DISTINCT gid from genes")

		gids = [tup[0] for tup in self.sql.fetchall()]

		# get list of all gids

		for gid in gids:

			self.sql.execute("SELECT rsid, alt from alleles where hgvs like '%=%' and rsid like '%rs%' and gid=?", (gid,))

			reference = { rsid : alt for (rsid, alt) in self.sql.fetchall() }
			
			# get list of all hapids for this gene
			
			self.sql.execute("SELECT DISTINCT hapid, starname, hgvs from alleles where gid=?", (gid,))
			
			hapids = self.sql.fetchall()
			
			for (hapid, starhap, hgvs) in hapids:
					
				# get haplotype alleles and create complete dictionary
				
				self.sql.execute("SELECT rsid, alt from alleles where hapid=? and rsid like '%rs%'", (hapid,))
				
				haprsids = { rsid : alt for (rsid, alt) in self.sql.fetchall()}

				haprsids = dict(reference, **haprsids)
				
				# get patient rsids
				
				self.sql.execute("select distinct v.rsid, p.alt from variants v join patientvars p on p.start=v.start join alleles a on v.gid=a.gid where p.call = '1/1' and a.hapid = ?", (hapid,))
				
				patrsids_base = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
				
				patrsids_al1 = dict(reference, **patrsids_base)
				
				modified = {k : (patrsids_al1[k], haprsids[k]) for k in patrsids_al1.keys() if patrsids_al1[k] != haprsids[k]}

				score1 = len(modified)

				# --------------------------------------------------------------------------------------------------
				
				self.sql.execute("select distinct v.rsid, p.alt from variants v join patientvars p on p.start=v.start join alleles a on v.gid=a.gid where p.call = '0/1' and a.hapid = ?", (hapid,))

				patrsids_add = { rsid : alt for (rsid, alt) in self.sql.fetchall()}

				patrsids_al2 = dict(patrsids_base, **patrsids_add)
				
				modified = {k : (patrsids_al2[k], haprsids[k]) for k in patrsids_al2.keys() if patrsids_al2[k] != haprsids[k]}
									
				score2 = len(modified)

				if len(patrsids_base) == 0 and len(patrsids_add) == 0:
								
					continue
				
				else:
					item = (gid, hapid, score1, score2)
								
					self.sql.execute("INSERT INTO patienthaps VALUES(?,?,?,?)", item)

		# -------------------------------------------------------------------------------------------------------------------------------------------------------

		self.conn.commit()
		
		
pat = Patient('data/test.g.vcf.gz')

