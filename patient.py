#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import sqlite3
import json
from collections import Counter
from requests_oauthlib import OAuth2Session
import urllib
import urllib2
import ast
from operator import itemgetter
from data import DataCollector


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

class Patient:

    '''
    This class imports a design VCF and an input VCF.
    '''

    def __init__(self, f):
        """
        f = given VCF input file
        GetDesign imports the test design and matches rs values to positions
        Import imports the input vcf
        GetIDs finds changed positions and matches them to the design.
        """

        self.f = f
        self.reader = vcf.Reader(open(f, 'r'))
        self.genes = {}
        self.haplotypes = []
        self.starhaps = {}
        self.Load()
        print 'Initiating patient...'

    def Load(self):
        """
        Function for loading the most important startup functions
        :return:
        """

        self.d = DataCollector('config/corrected_design.vcf')
        self.ImportData()
        self.GetIDs()  # commit db changes
        self.conn.commit()
        self.Hapmaker()
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
                                                (pos int PRIMARY KEY, num text, chr1 text, chr2 text)'''
                             )
            self.sql.execute('''CREATE TABLE patienthaps
                                                (gid text PRIMARY KEY, hgvs text, hapid text, starhap text)'''
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

        # create list for storing changed positions / rs#s

        for record in self.reader:

            # check for non-None positions, these have changed

            if '[None]' not in str(record.ALT):
                ref = record.REF
                alt = str(record.ALT[0])

                # go through record and fetch relevant data

                for sample in record.samples:
                    call = str(sample['GT'])  # 1/0 etc, phasing found here
                try:

                    # filter on call - implement phasing here

                    if call == '0/1' or call == '0|1':
                        chr1 = ref
                        chr2 = alt
                    elif call == '1/0' or call == '1|0':
                        chr1 = alt
                        chr2 = ref
                    elif call == '1/1' or call == '1|1':
                        chr1 = alt
                        chr2 = alt

                    # create sql entry from these variables

                    item = (record.POS, call, chr1, chr2)
                    self.sql.execute('''INSERT INTO patientvars VALUES(?,?,?,?)'''
                            , item)
                except KeyError:

                    # if no match is found, let user know
                    # please create design to aviod this

                    print "couldn't match", record.POS
                    pass

    def Hapmaker(self):
        '''
        IN PROGRESS:
        Create HGVS haplotypes from rsids. Should implement phasing...
        :return:
        '''

        # get gids that are present in patient by joining tables

        self.sql.execute("""SELECT DISTINCT v.gid from variants v
                            JOIN design d ON v.rsid = d.rsid
                            JOIN patientvars p ON d.pos = p.pos""")

        # make an easy usuable list for later on, remove tuples

        gids = [tup[0] for tup in self.sql.fetchall()]

        # loop through gids

        for gid in gids:

            # get rsid, alias for that rsid, and bases at this position in chr1 and 2

            self.sql.execute('''SELECT d.rsid, a.alias, p.chr2, p.num
                                FROM patientvars p
                                JOIN design d ON d.pos = p.pos
                                JOIN variants v ON d.rsid = v.rsid
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

            for (rsid, alias, chr2, num) in results:

                # select on "NC" and genomic position

                if 'NC' in alias and 'g.' in alias:
                    if '>' + chr2 in alias and chosen_ver in alias:
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
                                JOIN design d ON v.rsid = d.rsid
                                JOIN patientvars p ON d.pos = p.pos""")

        gids = [tup[0] for tup in self.sql.fetchall()]

        # go through gids

        for gid in gids:

            # check for searchability in pharmgkb, if not skip for now (should put an else later)

            if 'PA' not in gid:
                continue

            # get rsids and bases @ location for this gene that are present in patient,
            # save this result in patientrsids, used for comparison with haplotype later on

            self.sql.execute('''SELECT d.rsid, p.num, p.chr1, p.chr2
                                FROM patientvars p
                                JOIN design d ON d.pos = p.pos
                                JOIN variants v ON d.rsid = v.rsid
                                WHERE v.gid = ?;''',
                             (gid, ))

            results = self.sql.fetchall()

            # storage for rsids per chromosome

            patient_rsids = {'chr1': [], 'chr2': []}

            # sort rsids in right chromosome list

            for result in results:

                rsid = result[0]
                num = result[1]
                chr1 = result[2]
                chr2 = result[3]

                # sort rsids based on phased-ness (phasing assumed here, as used in pharmgkb star typing)
                # include alt base in comparison, will lower results but is more accurate perhaps

                if '1/' in num or '1|' in num:
                    patient_rsids['chr1'].append(rsid)

                if '/1' in num or '|1' in num:
                    patient_rsids['chr2'].append(rsid)

            # for given gid, find the known haplotypes (pharmgkb id and "star name" nomenclature)

            self.sql.execute('SELECT DISTINCT hapid, starname FROM alleles WHERE gid = ?'
                             , (gid, ))
            haps = self.sql.fetchall()

            # loop through these haplotypes, they will be compared to patient

            chr_matches = ([], [])

            for (i, hap) in enumerate(haps):

                hapid = hap[0]  # haplotype pharmgkb id
                starname = hap[1]  # star or other name, general name

                # check if it is not reference allele (should add a check for "normal activity" names)

                if i == 0:
                    self.sql.execute('SELECT rsid, alt FROM alleles WHERE hapid = ?'
                            , (hapid, ))
                    refrsids = self.sql.fetchall()
                    refrsids = [tup[0] for tup in refrsids]
                elif i >= 1:

                    # get rsids for haplotype

                    self.sql.execute('SELECT rsid, alt FROM alleles WHERE hapid = ?'
                            , (hapid, ))
                    haprsids = self.sql.fetchall()
                    haprsids = [tup[0] for tup in haprsids]

                    for (k, v) in patient_rsids.items():
                        preselect = list(set(v) & set(refrsids))
                        v = preselect

                    # compare patient with rsidlist and calculate match score per chromosome

                    comparison_chr1 = set(patient_rsids['chr1']) \
                        & set(haprsids)
                    comparison_chr2 = set(patient_rsids['chr2']) \
                        & set(haprsids)

                    # store these comparisons for use in a loop later on

                    comparisons = (comparison_chr1, comparison_chr2)

                    # go through the two comparisons, calculate match score

                    for (i, comparison) in enumerate(comparisons):

                        # calculate match score

                        match_score = float(len(comparison)) \
                            / len(haprsids)

                        # with a full match, save to corresponding list in dictionary

                        if match_score == 1.0:

                            # fetch useful data on resulting haplotypes

                            result = (hapid, starname, len(haprsids))

                            # check if results are already present, if not add

                            if result not in chr_matches[i]:
                                chr_matches[i].append(result)

                                # to do here... ref to hapmaker perhaps? can use two modes, one for starname and one for new HGVS

            for (i, result) in enumerate(chr_matches):

                # get non-empty results

                if len(result) == 0:
                    continue

                # find longest matching haplotype

                highest_hit = max(result, key=itemgetter(2))
                hapid = highest_hit[0]
                star_alt = highest_hit[1]

                if 'rs' not in star_alt:
                    starhap = num.replace('1', star_alt).replace('0',
                            '*1')
                elif 'rs' in star_alt:
                    numlist = num.split('/|')
                    starhap = num.replace('1', star_alt)

                # get haplotype pharmgkb id for this hit

                hapid = highest_hit[0]
                self.sql.execute('''SELECT DISTINCT gid FROM alleles WHERE hapid = ?'''
                                 , (hapid, ))
                gid = self.sql.fetchone()[0]
                self.sql.execute('UPDATE patienthaps SET hapid = ?, starhap = ? WHERE gid = ?'
                                 , (hapid, starhap, gid))