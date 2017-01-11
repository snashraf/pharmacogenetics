#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import sqlite3
from collections import Counter, OrderedDict
from modules.pgkb_functions import seqMaker
from db import Database
import subprocess as s
from Bio import Phylo
import re

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

        Database.__init__(self, dbname)

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

        self.Hapmatcher()

        self.conn.commit()


    def ImportData(self):
        
        """
        Import data from collected database and
        create patientvars tables.
        :return:
        """

        self.remakeTable("patientvars")


    def GetIDs(self):
        
        """
        Gets patient variables, reads into patientvars table
        :return:
        """
        
        print 'Reading patient variants...'

        # create list of positions to localize in patient

        for record in self.reader:

            self.sql.execute("SELECT DISTINCT loc, start, end, ref, alt, muttype FROM variants")

            positions = self.sql.fetchall()

            for (loc, start, end, ref, alt, muttype) in positions:

                records = self.reader.fetch(str(loc), start-1, end=end)

                # TODO PLEASE DO NOT DO THIS

                for record in records: # doctest: +SKIP

                    ref = str(record.REF)

                    start = record.POS

                    end = start + 1

                    alt = (str(record.ALT[0])).replace("<NON_REF>",".")

                    for sample in record.samples:

                        call = str(sample['GT'])  # 1/0 etc, phasing found here

                        try:

                            pid = str(sample['PID'])

                            pgt = str(sample['PGT'])

                        except KeyError:

                            pid = "nyan"

                            pgt = "nyan"

                    item = (loc, start, end, ref, alt, call, pid, pgt)

                    self.insertValues("patientvars", item)


    def Hapmatcher(self):

        """
        Matches rsids per gene to known haplotypes for that gene.
        Results are stored in a list...
        :return:
        """

        self.remakeTable('patienthaps')

        self.sql.execute("SELECT DISTINCT gid from genes")

        gids = [tup[0] for tup in self.sql.fetchall()]

        # get list of all gids

        for gid in gids:

            sequences = []

            self.sql.execute("SELECT a.rsid, a.alt, a.hapid from alleles a join variants v on a.rsid=v.rsid where a.hgvs like '%=%' and a.gid=? order by v.start", (gid,))

            refrsids = self.sql.fetchall()

            if not refrsids[0][2]:

                continue

            else:

                refid = refrsids[0][2]

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
