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
		self.HapScorer()

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

		# get gids from variants that are present in patient, fetch from design using pos

		self.sql.execute("""SELECT DISTINCT v.gid from variants v
								JOIN patientvars p ON p.start = v.start""")

		gids = [tup[0] for tup in self.sql.fetchall()]

		#print gids
		# go through gids

		for gid in gids:

			# check for searchability in pharmgkb, if not skip for now (should put an else later)
			# get rsids and bases @ location for this gene that are present in patient,
			# save this result in patientrsids, used for comparison with haplotype later on

			self.sql.execute('''SELECT v.rsid, p.call, p.ref, 
									   p.alt, d.starhaps
								FROM patientvars p
								JOIN variants v ON p.start = v.start
								JOIN drugpairs d on v.gid = d.gid
								WHERE v.gid = ? AND v.ref = p.ref 
								AND p.alt LIKE v.alt;''',
							 (gid, ))

			results = self.sql.fetchall()

			# storage for rsids per chromosome

			patient_rsids = {'chr1': [], 'chr2': [], "amb":[]}

			# sort rsids in right chromosome list

			for (rsid, call, ref, alt, starhaps) in results:

				# sort rsids based on phased-ness (phasing assumed here, as used in pharmgkb star typing)
				# include alt base in comparison, will lower results but is more accurate perhaps

				if "1/" in call:

					patient_rsids['chr1'].append(rsid)

				if "/1" in call:

					patient_rsids['chr2'].append(rsid)

				else:

					pass

			# for given gid, find the known haplotypes (pharmgkb id and "star name" nomenclature)

			self.sql.execute('SELECT DISTINCT hapid, hgvs, starname FROM alleles WHERE gid = ?'
							 , (gid, ))

			haps = self.sql.fetchall()

			# loop through these haplotypes, they will be compared to patient

			chr_matches = {}

			for i, (hapid, hgvs, starname) in enumerate(haps):

				# check if it is not reference allele (should add a check for "normal activity" names)

				if "=" in hgvs:

					self.sql.execute('SELECT rsid, alt FROM alleles WHERE hapid = ?'
							, (hapid, ))

					refrsids = self.sql.fetchall()

					refrsids = [tup[0] for tup in refrsids if "rs" in tup[0]]

				elif "[0]" in hgvs:

					print "deletion haplotype, skipping..."

					continue

				else:
					# get rsids for haplotype

					self.sql.execute('SELECT rsid, alt FROM alleles WHERE hapid = ?'
							, (hapid, ))

					haprsids = self.sql.fetchall()

					haprsids = [tup[0] for tup in haprsids if "rs" in tup[0]]

					if len(haprsids) == 0:

						continue

					# compare patient with rsidlist and calculate match score per chromosome

					for chr, rsids in patient_rsids.items():

						if len(rsids) == 0:

							continue

						comparison = set(rsids) & set(haprsids)

						# calculate match score

						match_hap = float(len(comparison)) \
							/ len(haprsids)
							
						match_pat = float(len(comparison))/len(rsids)
							
						item = (gid, hapid, chr, match_hap, match_pat)

						self.sql.execute('''INSERT INTO patienthaps VALUES(?,?,?,?,?)'''
							 , item)
						# insert THIS INTO PATIENTHAPS!!!
						
						#use this to fetch, do another join on patienthaps?
						#select distinct a.hapid, a.starname from alleles a join drugpairs d on d.gid = a.gid where d.starhaps like ("%" || a.starname || "%") and d.starhaps like "%*%";
					   

	def HapScorer(self):
		
		chrs = ["chr1", "chr2"]
		
		self.sql.execute("select distinct gid, symbol from genes")
		
		results = self.sql.fetchall()

		for (gid, symbol) in results:

			for chr in chrs:

				self.sql.execute("select distinct p.hapid, a.starname, p.hapscore, p.patscore from patienthaps p join alleles a \
				on a.gid=p.gid where p.chr=? and p.gid=? and p.hapscore != 0 and p.patscore != 0 order by hapscore, patscore limit 5", (chr,gid))
				
				results = self.sql.fetchall()

				if len(results) == 0:
				
					continue
				
				print symbol
				
				print chr
				
				print "---------------------"
				
				for (hapid, starname, hapscore, patscore) in results:
					
					lastval = "None"
					
					self.sql.execute("select distinct starhaps, guid from drugpairs where gid = ?", (gid,))
					
					vals = self.sql.fetchone()
					
					starhaps = vals[0]
					
					guid = vals[1]
					
					match = "| Partial |"
					
					if hapscore == 1:
					
						match = "| Full |"				
					
					if starname in starhaps and "nan" not in guid:
					
						lastval = guid 
					
					print hapid, "/", starname, "|", hapscore, patscore, "|", lastval, match
					
				
				print "---------------------"
		
						

# ----------------------------------------------------------------------------------------------------------

pat = Patient('data/test.g.vcf.gz')
pat.HapScorer()
