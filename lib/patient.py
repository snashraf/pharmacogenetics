#!/usr/bin/python
# -*- coding: utf-8 -*-

from data import *
import os
import vcf
from tqdm import tqdm
from modules.pgkb_functions import *
import subprocess as s
import pprint as pp
from Bio import Phylo

# ---------------------------------------------------------------------

class Patient(Database):

	'''
	This class imports a design VCF and an input VCF.
	'''

	def __init__(self, dbname, f):

		"""
		f = given VCF input file
		GetDesign imports the test design and matches rs values to
		positions
		Import imports the input vcf
		GetIDs finds changed positions and matches them to the design.
		"""

		print 'Initiating patient...'

		self.path = os.path.dirname(__file__)

		dbfolder = os.path.join(self.path, 'db')
		
		dbpath = os.path.join(dbfolder, '%s.db' % dbname)
		
		Database.__init__(self, dbpath)

# -------------------------------------------------------------
	
		self.f = f

		self.reader = vcf.Reader(open(f, 'r'))

		self.Load()


	def Load(self):

		"""
		Function for loading the most important startup functions
		:return:
		"""

		#self.GetSNPs()
	
		self.conn.commit()
	
		#self.AnnotateSNPs()
	
		self.Hapmatcher()
	
	   #self.conn.commit()

# --------------------------------- SNP parsing ----------------------------------------

	def GetSNPs(self):
	
		"""
		Gets patient variables, reads into patientvars table
		:return:
		"""
	
		self.remakeTable("patientvars")
	
	
			# create list of positions to localize in patient
	
		self.sql.execute('''
						SELECT DISTINCT l.Chromosome, l.Start, l.End, l.RefAllele, v.MutType
						FROM LocPGKB l
						JOIN Variants v
						ON v.VarID = l.VarID
						WHERE v.MutType = "snp"
						''')
	
		positions = self.sql.fetchall()
	
		print "Fetching patient variants... o(* ^ *o)"
		
		for (loc, start, end, ref, muttype) in tqdm(positions):
								
			records = self.reader.fetch(str(loc.lstrip("chr")), start-1, end=end)
			
			# TODO PLEASE DO NOT DO THIS

			for record in records: # doctest: +SKIP

				sql = self.insertSQL("patientvars").render(record = record)
			
				self.sql.executescript(sql)
		

	def AnnotateSNPs(self):

		self.authobj = Authenticate()
		
		self.remakeTable("patannotations")

		self.sql.execute('''
						SELECT v.VarID, p.CallBase, a.AnID
						FROM LocPGKB l
						JOIN PatientVars p
						ON p.Start = l.Start
						JOIN Variants v on l.VarId = v.VarId
						JOIN Annotations a
						ON a.VarHapId = v.VarId
						''')
						
		print "Annotating SNPs... /(* ` ^ `*/)"
		
		for (VarID, CallBase, AnID) in tqdm(self.sql.fetchall()):
			
			uri = \
			'https://api.pharmgkb.org/v1/data/clinicalAnnotation/{}?view=max' \
			.format(AnID)

			data = getJson(uri, self.authobj)
			
			allele = CallBase.replace("/", "")
			
			sql = self.insertSQL("patannotations").render(json = data, patallele = allele)
		
			self.sql.executescript(sql)

			# FOR SOMETHING THAT CREATES OVERVIEWS (new class? GUI? webserv?)
			
			overviewQuery = '''
			select distinct d.chemname, g.genename, p.phenotype, p.patallele, loe
			from patannotations p
			join annotations a
			on p.anid = a.anid
			join pairs p
			on a.drugid = p.drugid
			join drugs d
			on p.drugid = d.drugid
			join genes g
			on g.geneid = p.geneid
			order by g.genename;
			'''

# ---------------------------- Indel Parsing -------------------------------------



# ---------------------------- Haplotype parsing -------------------------------


	def Hapmatcher(self):

			"""
			Matches rsids per gene to known haplotypes for that gene.
			Results are stored in a list...
			:return:
			"""
	
			self.remakeTable('pathaplotypes')
	
			self.sql.execute("SELECT DISTINCT GeneID from Genes")
	
			gids = [tup[0] for tup in self.sql.fetchall()]

			# create view with everything
			self.sql.executescript('''
						DROP VIEW Overview;

						CREATE VIEW Overview AS
						    SELECT hap.geneid AS geneid,
						           gen.genesymbol AS genesymbol,
						           h.hapid AS hapid,
						           varname AS varname,
						           h.altallele AS hapallele,
						           pat.start AS start,
						           pat.refallele AS pat_ref,
						           pat.altallele AS pat_alt,
						           pat.callNum AS call,
						           hap.hapname AS hapname,
						           hap.starname AS starname,
						           hap.hgvs AS hgvs
						      FROM hapvars h
						           JOIN
						           variants v ON v.rsid = h.VarName
						           JOIN
						           locvcf loc ON loc.varid = v.varid
						           JOIN
						           patientvars pat ON pat.start = loc.start
						           JOIN
						           haplotypes hap ON hap.hapid = h.HapID
						           JOIN
						           genes gen ON gen.geneid = hap.geneid;
								''')
	
			# get list of all gids
	
			print "Haplotyping patient... (\\' n')\\*scribble*"

			for gid in tqdm(gids):
	
				sequences = []

				self.sql.execute('''
								select varname, hapname, hapid, hapallele
								from overview where geneid = ?
								and hgvs like "%=%"
								order by start asc
								''', (gid,))
	
				refrsids = self.sql.fetchall()

				if len(refrsids) == 0:
					continue

				refid = refrsids[0][2]
				
				rsidorder = []
				reference = {}

				for (rsid, hap_name, hapid, hap_ref) in refrsids:

					rsidorder.append(rsid)
					
					reference[rsid] = hap_ref

				refseq = seqMaker(rsidorder, reference, reference)
	
				sequences.append(">{}\n{}\n".format(refid, refseq))
	
				# get list of all hapids for this gene
	
				self.sql.execute('''
								SELECT DISTINCT HapID, hapname, HGVS
								from Overview 
								where GeneID = ? 
								and hgvs not like "%=%"
								''', 
								(gid,))

				selection = self.sql.fetchall()
	
				for (hapid, starhap, hgvs) in selection:
	
					# get haplotype alleles and create complete dictionary
	
					self.sql.execute('''SELECT VarName, hapallele 
										from Overview 
										where HapID = ?  
										''', (hapid,))
	
					haprsids = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
	
					if len(haprsids) == 0:
	
						continue
	
					else:
	
						haprsids = dict(reference, **haprsids)
	
						hapseq = seqMaker(rsidorder, reference, haprsids)
	
						modified = {k : (reference[k], haprsids[k]) for k in reference.keys() if reference[k] != haprsids[k]}
	
						score_ref = len(modified)
	
						if score_ref == 0:
	
							continue
	
						else:
	
							# get patient rsids

							self.sql.execute('''SELECT VarName, pat_alt 
												from Overview 
												where HapID = ?
												and call = "1/1"  
											''', (hapid,))
	
							patrsids_base = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
	
							patrsids_al1 = dict(reference, **patrsids_base)
	
							patseq1 = seqMaker(rsidorder, reference, patrsids_al1)
	
							modified = {k : (patrsids_al1[k], haprsids[k]) for k in patrsids_al1.keys() if patrsids_al1[k] != haprsids[k]}
	
							score_al1 = len(modified)
	
							# --------------------------------------------------------------------------------------------------
							
							self.sql.execute('''SELECT VarName, pat_alt 
												from Overview 
												where HapID = ?
												and call = "0/1"  
											''', (hapid,))
	
							patrsids_add = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
	
							patrsids_al2 = dict(patrsids_base, **patrsids_add)
	
							patseq2 = seqMaker(rsidorder, reference, patrsids_al2)
	
							modified = {k : (patrsids_al2[k], haprsids[k]) for k in patrsids_al2.keys() if patrsids_al2[k] != haprsids[k]}
	
							score_al2 = len(modified)
	
							if len(patrsids_base) == 0 and len(patrsids_add) == 0:
	
								continue
	
							else:
	
								hapline = ">%s\n %s\n" %(hapid, hapseq)
	
								if hapline not in sequences and hapseq != "":
	
									sequences.append(hapline)
	
								else:
	
									continue

				sequences.append(">Patient_allele1\n" + patseq1 + "\n")
	
				sequences.append(">Patient_allele2\n" + patseq2 + "\n")
				
				output = "/output/alignments/"
	
				fn = self.path + output + gid + "_aln.fasta"
	
				with open(fn, "w") as f:
	
					f.writelines(sequences)
		
				try:
	
					self.HapScorer(fn, "phylo", refid)
	
				except:
	
					raise


	def HapScorer(self, fn, mode, refid):

		if mode == "phylo":

			# phylogenetic tree

			of = fn.replace("alignments/", "alignments/aligned/")

			tn = of.strip(".fasta")+"_tree.dnd"

			with open(fn, "rb") as infile, open(of, "wb") as outfile:

				s.check_call("{}/plugins/clustalo -i {} -o {} --auto --force --guidetree-out={}"
					.format(self.path, fn, of, tn),  shell=True)

			tree = Phylo.read(tn, "newick")

			names = []

			# modified

			tree.root_with_outgroup({refid})
			
			for clade in tree.find_clades():

				if clade.name and clade.name not in names:

						names.append(clade.name)

			for hap in names:

				distances = {"hapid":hap}

				if "Patient" in hap:

					continue

				else:

					for i in range(1,3):

						dist = tree.distance("Patient_allele%i" %i, hap)

						distances["al%i" %i] = dist

					sql = self.insertSQL("pathaplotypes").render(json = distances)
				
					self.sql.executescript(sql)

					# commit to db

			self.conn.commit()

		else:
			# standard mode
			pass

		self.conn.commit()
