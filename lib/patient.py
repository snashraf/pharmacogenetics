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
from collections import OrderedDict
from interpreter import Interpreter

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

# --------------------------------- SNP parsing --------------------------

    def GetPositions(self):
        """
        Gets patient variables, reads into patientvars table
        :return:
        """

        self.remakeTable("patvars")

        # create list of positions to localize in patient

        self.sql.execute('''
			SELECT DISTINCT l.Chromosome, l.Start, l.End, l.RefAllele, v.MutType
			FROM LocVCF l
			JOIN Variants v
			ON v.VarID = l.VarID
			''')

        positions = self.sql.fetchall()

        print "Fetching patient variants... o(* ^ *o)"

        for (loc, start, end, ref, muttype) in tqdm(positions):

            records = self.reader.fetch(
                str(loc.lstrip("chr")), start - 1, end=end)

            # TODO PLEASE DO NOT DO THIS

            for record in records:  # doctest: +SKIP

                sql = self.insertSQL("patvars").render(record=record)

                self.sql.executescript(sql)


# ---------------------------- Haplotype parsing -------------------------

    def GetHaplotypes(self):
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
					SELECT  DISTINCT
					
					hap.GeneID as GeneID,
					hap.HapID as HapID,
					hap.HGVS as HGVS,
					hap.starname as starname,
					VarName, MutType,
					h.AltAllele as HapAllele,
					locb.RefAllele as RefPGKB,
					a.AltPGKB as AltPGKB,
					locv.RefAllele as RefVCF,
					a.AltVCF as AltVCF,
					pat.RefAllele as PatRef,
					pat.AltAllele as PatAlt,
					CallNum,
					locv.start as start
					
					FROM hapvars h
					   JOIN
					   variants v ON v.rsid = h.VarName
					   JOIN
					   locvcf locv ON locv.varid = v.varid
					   JOIN
					   patientvars pat ON pat.start = locv.start
					   JOIN
					   haplotypes hap ON hap.hapid = h.HapID
					   JOIN
					   genes gen ON gen.geneid = hap.geneid
					   JOIN
					   locpgkb locb ON locb.varid = locv.varid
					   JOIN
					   altalleles a on a.VarID = locv.VarID
					   ORDER BY locv.start asc;
							''')
<<<<<<< HEAD

        # get list of all gids

        print "Haplotyping patient... (\\' n')\\*scribble*"

        for gid in tqdm(gids):

            print "--------------%s--------------" % gid
            sequences = []

            # First: collect SNPs

            self.sql.execute('''
					SELECT DISTINCT VarName
=======
	
			# get list of all gids
	
			print "Haplotyping patient... (\\' n')\\*scribble*"

			for gid in tqdm(gids):
	
					
				varValues = OrderedDict([])

				# First: collect SNPs and patient positions at those SNPs
				
				self.sql.execute('''
					SELECT DISTINCT hapid, VarName, PatAlt, CallNum
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
					FROM Overview
					WHERE GeneID = ?
					AND HGVS LIKE "%[=]%"
					ORDER BY Start ASC
					''', (gid,))
<<<<<<< HEAD

            rsidorder = [rsid for rsid in self.sql.fetchall()]

            print rsidorder

            # Fetch SNPs and add to dictionary

            self.sql.execute('''
=======
					
				refRsids = self.sql.fetchall()
	
				rsidorder = [rsid for (hapid, rsid, patAlt, CallNum) in refRsids]
				
				if len(rsidorder) == 0:
					
					continue
					
				#print rsidorder
				
				patrsids_hom = { rsid : patAlt for hapid, rsid, patAlt, CallNum in refRsids if CallNum == "1/1"}
				
				patrsids_het = { rsid : patAlt for hapid, rsid, patAlt, CallNum in refRsids if CallNum == "0/1"}
							
				refid = refRsids[0][0]
					
				# Fetch SNPs and add to dictionary
				
				self.sql.execute('''
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
					SELECT DISTINCT VarName, HapAllele
					from overview
					where geneid = ?
					and hgvs like "%[=]%"
					and muttype = "snp"
					order by start asc
					''', (gid,))

<<<<<<< HEAD
            snps = {rsid: allele for (rsid, allele) in self.sql.fetchall()}

            print snps

            # Fetch indels that match RefPGKB and bind them to RefVCF

            self.sql.execute('''
=======
				snps = {rsid:allele for (rsid, allele) in self.sql.fetchall()}
								
				# Fetch indels that match RefPGKB and bind them to RefVCF
				
				self.sql.execute('''
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
					SELECT DISTINCT VarName, RefVCF
					FROM overview where geneid = ?
					and hgvs like "%[=]%"
					and muttype != "snp"
					and hapallele = refpgkb
					order by start asc
					''', (gid,))
<<<<<<< HEAD

            indels_ref = {rsid: allele for (
                rsid, allele) in self.sql.fetchall()}

            print indels_ref

            # Fetch indels that match AltPGKB and bind them to ALtVCF

            self.sql.execute('''
=======
				
				indels_ref  = {rsid:allele for (rsid, allele) in self.sql.fetchall()}
								
				# Fetch indels that match AltPGKB and bind them to ALtVCF
				
				self.sql.execute('''
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
					SELECT DISTINCT VarName, AltVCF
					FROM overview where geneid = ?
					and hgvs like "%[=]%"
					and muttype != "snp"
					and hapallele = altpgkb
					order by start asc
					''', (gid,))
<<<<<<< HEAD

            indels_alt = {rsid: allele for (
                rsid, allele) in self.sql.fetchall()}

            print indels_alt

            # Join these three dictionaries and feed to seqMaker

            reference = merge_dicts(snps, indels_ref, indels_alt)

            print reference

            raw_input("Press enter to continue...")

            continue

            # SeqMaker makes a fictive piece of DNA consisting of the haplotype
            # defining RSIDS

            refseq = seqMaker(rsidorder, reference, reference)

            sequences.append(">{}\n{}\n".format(refid, refseq))

            # get list of all hapids for this gene

            self.sql.execute('''
					SELECT DISTINCT HapID, hapname, HGVS
=======
				
				indels_alt  = {rsid:allele for (rsid, allele) in self.sql.fetchall()}
								
				# Join these three dictionaries and feed to seqMaker
				
				varValues[refid] = merge_dicts(snps, indels_ref, indels_alt)
												
				varValues['a1'] = dict(varValues[refid], **patrsids_hom)

				varValues['a2'] = dict(varValues[refid], **patrsids_het)

				# SeqMaker makes a fictive piece of DNA consisting of the haplotype defining RSIDS
				
				
				# ------------------------- create patient sequences -----------------------------
	
				# get list of all hapids for this gene
	
				self.sql.execute('''
					SELECT DISTINCT HapID, starname
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
					from Overview
					where GeneID = ?
					and hgvs not like "%=%"
					''',
<<<<<<< HEAD
                             (gid,))

            selection = self.sql.fetchall()

            for (hapid, starhap, hgvs) in selection:

                # get haplotype alleles and create complete dictionary

                self.sql.execute('''SELECT VarName, AltVCF
=======
					(gid,))

				selection = self.sql.fetchall()
	
				for (hapid, starname) in selection:
	
					# get haplotype alleles and create complete dictionary
	
					self.sql.execute('''SELECT VarName, AltVCF
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
							from Overview
							where HapID = ?
							and HapAllele = AltPGKB
							''', (hapid,))
<<<<<<< HEAD

                haprsids = {rsid: alt for (rsid, alt) in self.sql.fetchall()}

                if len(haprsids) == 0:

                    continue

                else:

                    haprsids = dict(reference, **haprsids)

                    hapseq = seqMaker(rsidorder, reference, haprsids)

                    modified = {k: (reference[k], haprsids[k]) for k in reference.keys(
                    ) if reference[k] != haprsids[k]}

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

                        patrsids_base = {rsid: alt for (
                            rsid, alt) in self.sql.fetchall()}

                        patrsids_al1 = dict(reference, **patrsids_base)

                        patseq1 = seqMaker(rsidorder, reference, patrsids_al1)

                        modified = {k: (patrsids_al1[k], haprsids[k]) for k in patrsids_al1.keys(
                        ) if patrsids_al1[k] != haprsids[k]}

                        score_al1 = len(modified)

                        # --------------------------------------------------------------------------------------------------

                        self.sql.execute('''SELECT VarName, pat_alt
									from Overview
									where HapID = ?
									and call = "0/1"
									''', (hapid,))

                        patrsids_add = {rsid: alt for (
                            rsid, alt) in self.sql.fetchall()}

                        patrsids_al2 = dict(patrsids_base, **patrsids_add)

                        patseq2 = seqMaker(rsidorder, reference, patrsids_al2)

                        modified = {k: (patrsids_al2[k], haprsids[k]) for k in patrsids_al2.keys(
                        ) if patrsids_al2[k] != haprsids[k]}

                        score_al2 = len(modified)

                        if len(patrsids_base) == 0 and len(patrsids_add) == 0:

                            continue
=======
	
					haprsids = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
					
					if len(haprsids) == 0:
	
						continue
	
					else:
	
						varValues[hapid] = dict(varValues[refid], **haprsids)
			
				sequences = []
		
				#print sequences
								
				output = "/output/alignments/"
	
				fn = self.path + output + gid + "_aln.fasta"
				
				prev_seqs = []
				
				self.sql.executescript('''
					DROP VIEW HapView;
					
					CREATE VIEW HapView AS
					
					SELECT DISTINCT * FROM GuideOptions o
					JOIN Genes g
					on g.GeneSymbol = o.GeneSymbol
					JOIN Haplotypes h
					On G.GeneID = h.GeneID
					AND o.Starname = h.Starname
					''')
				
				self.sql.execute("SELECT DISTINCT HapID FROM HapView WHERE GeneID = ?", (gid,))
				
				options = [hapid for hapid in self.sql.fetchall()]
				
				with open(fn, "w") as f:
	
					refseq = seqMaker(rsidorder, varValues[refid], varValues[refid])

					prev_seqs.append(refseq)
					
					f.write(">{}\n{}\n".format(refid, refseq))

					for var, values in varValues.items():
						
						seq = seqMaker(rsidorder, varValues[refid], values)

						if (seq not in prev_seqs and var not in options) or (var == 'a1' or var == 'a2'):
								
							f.write(">{}\n{}\n".format(var, seq))

						prev_seqs.append(seq)
						
				self.HapScorer(fn, "phylo", refid)
	
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

                        else:

                            hapline = ">%s\n %s\n" % (hapid, hapseq)

                            if hapline not in sequences and hapseq != "":

                                sequences.append(hapline)

                            else:

                                continue

            sequences.append(">Patient_allele1\n" + patseq1 + "\n")

            sequences.append(">Patient_allele2\n" + patseq2 + "\n")

<<<<<<< HEAD
            output = "/output/alignments/"

            fn = self.path + output + gid + "_aln.fasta"
=======
			tree = Phylo.read(tn, "newick")
			
			names = []
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

            with open(fn, "w") as f:

<<<<<<< HEAD
                f.writelines(sequences)
=======
			tree.root_with_outgroup({refid})
			
			Phylo.draw_ascii(tree)

			for clade in tree.find_clades():
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

            try:

<<<<<<< HEAD
                self.HapScorer(fn, "phylo", refid)

            except:
=======
						names.append(clade.name)
			for hap in names:
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

                raise

    def HapScorer(self, fn, mode, refid):

        if mode == "phylo":

            # phylogenetic tree

            of = fn.replace("alignments/", "alignments/aligned/")

<<<<<<< HEAD
            tn = of.strip(".fasta") + "_tree.dnd"
=======
						dist = tree.distance(target1="a%i" %i, target2=hap)
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

            with open(fn, "rb") as infile, open(of, "wb") as outfile:

<<<<<<< HEAD
                s.check_call("{}/plugins/clustalo -i {} -o {} --auto --force --guidetree-out={}"
                             .format(self.path, fn, of, tn),  shell=True)

            tree = Phylo.read(tn, "newick")
=======
					sql = self.insertSQL("pathaplotypes").render(json = distances)
				
					self.sql.executescript(sql)
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

            names = []

            # modified

            tree.root_with_outgroup({refid})

<<<<<<< HEAD
            for clade in tree.find_clades():

                if clade.name and clade.name not in names:

                    names.append(clade.name)

            for hap in names:

                distances = {"hapid": hap}

                if "Patient" in hap:

                    continue

                else:

                    for i in range(1, 3):

                        dist = tree.distance("Patient_allele%i" % i, hap)

                        distances["al%i" % i] = dist

                    sql = self.insertSQL(
                        "pathaplotypes").render(json=distances)

                    print sql

                    self.sql.executescript(sql)

                    # commit to db

            self.conn.commit()

        else:
            # standard mode
            pass

        self.conn.commit()
=======
		self.conn.commit()
		
	def Interpret(self):
		
		# Check for reference alleles
		
		i = Interpreter(self)
		
		i.Genotyper()
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

# ------------------------ annotations ------------------------------

    def GetSNPAnnotations(self):

        self.authobj = Authenticate()

<<<<<<< HEAD
        self.remakeTable("patannotations")
=======
			# only when there is no haplotype available?
			
			self.authobj = Authenticate()
			
			self.remakeTable("patannotations")
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac

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

            sql = self.insertSQL("patannotations").render(
                json=data, patallele=allele)

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

            queryHaplotypeScores = '''
				select geneid, hapid, hapname, distance1, distance2 from pathaplotypes p
				join haplotypes h
				on p.HapID = h.hapid
				join genes g
				on h.geneid = g.geneid;
				'''

<<<<<<< HEAD
    def GetHapAnnotations(self):
        self.sql.executescript('''
=======

	def GetHapAnnotations(self):
		self.sql.executescript('''
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
			CREATE VIEW HapAnnotations AS
			    SELECT DISTINCT *
			      FROM pathaplotypes h
			           JOIN
			           annotations a ON a.varhapid = h.hapid
			           JOIN
			           haplotypes t ON t.hapid = h.HapID
			           JOIN
			           genes g ON g.geneid = t.geneid
			     WHERE hgvs NOT LIKE "%[0]%"
			     ORDER BY geneid;
						''')
