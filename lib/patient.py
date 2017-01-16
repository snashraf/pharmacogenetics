#!/usr/bin/python
# -*- coding: utf-8 -*-

from data import *
import os
import vcf
from tqdm import tqdm
from modules.pgkb_functions import *

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

        path = os.path.dirname(__file__)

        dbfolder = os.path.join(path, 'db')
		
	dbpath = os.path.join(dbfolder, '%s.db' % dbname)
	
	print dbpath

        Database.__init__(self, dbpath)

# -------------------------------------------------------------
	
        self.f = f

        self.reader = vcf.Reader(open(f, 'r'))

        print "Loading patient data..."

        self.Load()


    def Load(self):

        """
        Function for loading the most important startup functions
        :return:
        """

	       # self.GetSNPs()
	
	     #   self.conn.commit()
	
        self.AnnotateSNPs()
	
	 #       self.Hapmatcher()
	
	   #     self.conn.commit()


    def GetSNPs(self):
	
	        """
	        Gets patient variables, reads into patientvars table
	        :return:
	        """
	
	        self.remakeTable("patientvars")
	
	
		print 'Reading patient variants...'
	
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
						WHERE p.AltAllele NOT LIKE "<NON_REF>"
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


    def Hapmatcher(self):

	        """
	        Matches rsids per gene to known haplotypes for that gene.
	        Results are stored in a list...
	        :return:
	        """
	
	        self.remakeTable('patienthaps')
	
	        self.sql.execute("SELECT DISTINCT GeneID from Genes")
	
	        gids = [tup[0] for tup in self.sql.fetchall()]
	
	        # get list of all gids
	
	        for gid in gids:
	
	            sequences = []
	
	            self.sql.execute('''SELECT h.VarName, h.AltAllele, h.HapID
	                    from HapVars h
	                    join Variants v
	                    on h.VarName = v.RSID
	                    join Haplotypes t
	                    ON t.HapID = h.HapID
	                    where t.HGVS like '%=%'
	                    and h.GID=?
	                    order by v.start''', (gid,))
	
	            refrsids = self.sql.fetchall()
	
	            rsidorder = [rsid for (rsid, alt, hapid) in refrsids]
	
	            reference = { rsid : alt for (rsid, alt, hapid) in refrsids}
	
	            refseq = seqMaker(rsidorder, reference, reference)
	
	            newline =  ">{}\n{}\n".format(refid, refseq)
	
	            sequences.append(newline)
	
	            # get list of all hapids for this gene
	
	            self.sql.execute("SELECT DISTINCT hapid, starname, hgvs from alleles where gid=? and hgvs not like '%=%'", (gid,))
	
	            for (hapid, starhap, hgvs) in self.sql.fetchall():
	
	                # get haplotype alleles and create complete dictionary
	
	                self.sql.execute("SELECT rsid, alt from alleles where hapid=? and rsid like '%rs%'", (hapid,))
	
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
	
	                        self.sql.execute("select distinct v.rsid, p.alt from variants v \
	                                        join patientvars p on p.start = v.start \
	                                        join alleles a on v.gid = a.gid \
	                                        where p.call = '1/1' \
	                                        and a.hapid = ? \
	                                        and v.rsid like '%rs%'", (hapid,))
	
	                        patrsids_base = { rsid : alt for (rsid, alt) in self.sql.fetchall()}
	
	                        patrsids_al1 = dict(reference, **patrsids_base)
	
	                        patseq1 = seqMaker(rsidorder, reference, patrsids_al1)
	
	                        modified = {k : (patrsids_al1[k], haprsids[k]) for k in patrsids_al1.keys() if patrsids_al1[k] != haprsids[k]}
	
	                        score_al1 = len(modified)
	
	                        # --------------------------------------------------------------------------------------------------
	
	                        self.sql.execute("select distinct v.rsid, p.alt from variants v \
	                                            join patientvars p on p.start=v.start join alleles a on v.gid=a.gid \
	                                            where p.call = '0/1' and a.hapid = ? and v.rsid like '%rs%'", (hapid,))
	
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
	
	            path = "data/alignments/"
	
	            fn = path + gid + "_aln.fasta"
	
	            with open(fn, "w") as f:
	
	                f.writelines(sequences)
	
	            print fn
	
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

                s.check_call("plugins/clustalo -i %s -o %s --auto --force --guidetree-out=%s" %(fn, of, tn),  shell=True)

            tree = Phylo.read(tn, "newick")

            names = []

            for clade in tree.find_clades():

                if clade.name and clade.name not in names:

                        names.append(clade.name)

            tree.root_with_outgroup({refid})

            Phylo.draw_ascii(tree)

            distances = {}

            for hap in names:

                if "Patient" in hap:

                    continue

                else:

                    for i in range(1,3):

                        matches = []

                        dist = tree.distance("Patient_allele%i" %i, hap)

                        distances["al%i" %i] = dist

                    item = (hap, distances["al1"], distances["al2"])

                    self.insertValues("patienthaps", item)

            self.conn.commit()


        else:
            # standard mode
            pass

        self.conn.commit()
