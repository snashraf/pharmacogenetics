#!/usr/bin/python
# -*- coding: utf-8 -*-

import urllib2
import json
from variant import Variant
from gene import Gene
from drug import Drug
import sqlite3
from pgkb_functions import PGKB_connect, getRef
from tqdm import tqdm

# ---------------------------------------------------------------------

class DataCollector(object):

    '''
....This class imports a design VCF and creates a database with these positions, using PharmGKB.
....'''

    def __init__(self):
        """
........f = given VCF input file
........GetDesign imports the test design and matches rs values to positions
........Import imports the input vcf
........GetIDs finds changed positions and matches them to the design.
........All of this is exported
........"""

        self.conn = sqlite3.connect('pharmacogenetics.db')  # connect to db- if it doesn't exist, create it
        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

    def Update(self):
        """
........This function rebuilds the database, use for updating the db periodically?
........:return:
........"""

        self.GetPairs()
        self.conn.commit()
        self.GetGeneData()
        self.conn.commit()
        self.GetVarData()
        self.conn.commit()
        self.GetNonRS()
        self.conn.commit()
        self.GetChemData()
        self.conn.commit()

    def GetPairs(self):
        """
........Create table for drug-gene pairs, to be used later to fetch drug data
........:return:
........"""

        authobj = Authenticate()

        self.sql.execute('''DROP TABLE IF EXISTS drugpairs''')
        self.sql.execute('''CREATE TABLE drugpairs(did text, gid text, guid text, starhaps text, rsids text)'''
                         )

        print 'Getting gene-drug pairs... ~( * v*)/\\(^w ^)~'

        # get the uri for finding well-annotated pairs (given by pharmgkb)

        uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

        # get data and read in this json file

        data = urllib2.urlopen(uri)

        guid = 'nan'

        options = 'nan'

        for doc in tqdm(json.load(data)):

            gid = doc['gene']['id']

            symbol = doc['gene']['symbol']

            did = doc['chemical']['id']

            # get annotated variants first

            results = PGKB_connect(authobj, 'clinAnno', did, gid)

            varids = 'nan'

            if results is not None:

                varids = []

                for doc in results:

                    rsid = doc['location']['displayName']

                    varids.append(rsid)

            if type(varids) == list:

                varids = ','.join(list(set(varids)))

            results = PGKB_connect(authobj, 'clinGuide', did, gid)

            if results is not None:

                guid = results['guid']

                optionlist = results['options']

                if optionlist == None:

                    options = 'nan'

                else:

                    for gene in optionlist['data']:

                        if symbol in gene['symbol']:

                            options = ','.join(gene['options'])

            item = (did, gid, str(guid), options, varids)

            # insert results in table drugpairs

            self.sql.execute('''INSERT INTO drugpairs VALUES(?,?,?,?,?)'''
                             , item)

        self.conn.commit()

    def GetGeneData(self):
        '''
........Fetches data on given gene IDs.
........:return:
........'''

        print 'Getting gene data... (/o*)'

        # drop already existing tables genes and alleles

        self.sql.execute('''DROP TABLE IF EXISTS genes''')

        self.sql.execute('''DROP TABLE IF EXISTS alleles''')

        # (re)create tables in database

        self.sql.execute('''CREATE TABLE genes
                                            (gid text UNIQUE, symbol text, chr text, start text, stop text)'''
                         )

        self.sql.execute('''CREATE TABLE alleles
                                            (hapid text, gid text, starname text,
                                            hgvs text, rsid text, alt text,
                                            UNIQUE(hapid, rsid, starname)
                                            ON CONFLICT REPLACE)''')

        # get all unique gene ids from the variant table

        self.sql.execute('SELECT DISTINCT gid FROM drugpairs')

        genes = [tup[0] for tup in self.sql.fetchall()]

        genes = list(set(genes))

        # go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

        for gid in tqdm(genes):
            
            g = Gene(gid)

            # insert the resulting name and alleles into sql table

            self.sql.execute('''INSERT INTO genes VALUES(?,?,?,?,?)'''
                             , (gid, g.name, g.chr, g.start,
                             g.stop))

            for allele in g.alleles:

                for (rsid, alt) in allele['rsids']:

                    self.sql.execute('''INSERT INTO alleles VALUES(?,?,?,?,?,?)'''
                            , (
                        allele['id'],
                        gid,
                        allele['starname'],
                        allele['hgvs'],
                        rsid,
                        alt,
                        ))


    def DbPrep(self):

        # drop tables if they exist already to reset them

        self.sql.execute('''DROP TABLE IF EXISTS variants''')

        self.sql.execute('''DROP TABLE IF EXISTS alias''')

        self.sql.execute('''DROP TABLE IF EXISTS transtable''')

        # create variant and alias table. Variant should have an unique combo of rsid and gid,
        # and alias should be unique in the alias table.

        self.sql.execute('''CREATE TABLE variants
                                    (rsid text, varid text, gid text, chr text, 
                    pbegin int, pend int, pref text, palt text,
                    vbegin int, vend int, vref text, valt text,
                    type text,
                                    UNIQUE(rsid,gid)
                                    ON CONFLICT REPLACE)''')

        self.sql.execute('''CREATE TABLE alias
                                            (rsid text, alias varchar(255) PRIMARY KEY)'''
                         )

        self.conn.commit()


    def GetVarData(self):
        """
........Using Variant objects, this function fetches data on variants from the PharmGKB servers,
........if not available uses the Entrez servers, possibilities are limited for now.
........:return:
........"""

        print 'Getting variant data... ~(^_^)~'

        # get all rsids in the design vcf

        self.sql.execute('SELECT DISTINCT a.rsid, a.gid FROM alleles a JOIN genes g on a.gid = g.gid where rsid LIKE "rs%" order by a.gid, a.rsid'
                         )

        # rotate through rsids and create variant objects to fetch information

        for (rsid, gid) in tqdm(self.sql.fetchall()):

            # create variant instances with the given rsid.

            try:

                v = Variant(rsid, 'pharmgkb')
            
            except KeyError:

                continue

            # this results in a combination tuple of rsid and gid and aliases

            item = (
                rsid,
                v.id,
                gid,
                v.chr,
                v.begin,
                v.end,
                v.ref,
                v.alt,
                v.begin,
                v.end,
                v.ref,
                v.alt,
                v.type,
                )

            self.sql.execute('''INSERT INTO variants VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)'''
                             , item)

            # go through aliases, ignore duplicates and put in alias table

            for alias in v.names:

                try:
                    self.sql.execute('''INSERT INTO alias VALUES(?,?)'''
                            , (rsid, alias))
                except sqlite3.IntegrityError:

                    # on duplicate, ignore

                    continue

        # create entry in transtable

        if v.type != 'snp':

            # left shift position by 1

            self.TranslateIndel(v)


    def TranslateIndel(self, v):

        v.nbegin = v.begin - 1

        v.nend = v.end

        # get reference nucleotide at that position

        prevbase = getRef(v.chr, v.nbegin, v.nbegin)

        # set defaults

        alts = []

        if v.ref == '-':

        # scenario 1: insertion (REF - ALT A)

            v.nref = prevbase

            for alt in v.alt.split(','):

                alt = prevbase + alt

                alts.append(alt)

            v.nalt = ', '.join(alts)

        elif '-' in v.alt:

        # scenario 2: deletion ( REF A, ALT -, A)

            v.nref = prevbase + v.ref

            for alt in v.alt.split(','):

                if alt == '-':

                    alt = prevbase

                if alt == v.ref:

                    alt = v.nref
                else:

                    alts.append(alt)

        elif '(' in v.ref:

            # TA repeats

            # manual for now

            v.nref = prevbase + 'TA'

            # subtract ref TAs

            alts.append(prevbase)

            alts.append(prevbase + 'TATA')

            alts.append(prevbase + 'TATATA')

        # create and insert table item (5 columns)........

        if len(alts) > 0:

            for alt in alts:

                self.sql.execute('UPDATE variants SET vbegin = ?, vend = ?, vref = ?, valt = ? WHERE rsid = ?'
                        , (v.nbegin, v.nend, v.nref, alt, rsid))
        else:

            return


    def GetNonRS(self):

        self.sql.execute("SELECT rsid, alt, gid from alleles where rsid LIKE '%chr%' order by rsid"
                         )

        print 'Parsing non-rs variants...'

        for (rsid, alt, gid) in tqdm(self.sql.fetchall()):

                    # find genome version (hg19)

            ver = rsid[rsid.find('(') + 1:rsid.find(')')]

                         # find chromosome number

            loc = rsid.split(':')[0]

                        # find location

            begin = int(rsid[rsid.find(':') + 1:rsid.find('(')])

            end = begin

                        # get reference position

            ref = getRef(chr, begin, begin)

            if ref == alt or 'delGENE' in alt:

                continue

            if len(ref) == len(alt):

                muttype = 'snp'

                nbegin = begin

                nend = nbegin

                nref = ref

                nalt = alt

            if 'ins' in alt:

                nalt = alt.replace('ins', ref)

                muttype = 'in-del'

                nend = nbegin + len(alt)

            elif 'del' in alt:

                nbegin = nbegin - 1

                prev = getRef(chr, nbegin, nbegin)

                if alt == 'del':

                    nref = prev + nref
                
                else:

                    nref = prev + alt.lstrip('del')

                nalt = prev

                muttype = 'in-del'

                nend = nbegin + len(ref)
            
            else:

                type = 'unknown'

            item = (
                rsid,
                rsid,
                gid,
                loc,
                begin,
                end,
                ref,
                alt,
                nbegin,
                nend,
                nref,
                nalt,
                muttype,
                )

            self.sql.execute('''INSERT INTO variants VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)'''
                             , item)

    def GetChemData(self):
        """
........Gets info on seperate drugs, such as name and associated terms
........terms give info on what the drug does
........:return:
........"""

        print 'Getting drug data /(>_<)\\}'

        # drop and re-create table drugs

        self.sql.execute('''DROP TABLE IF EXISTS drugs''')

        self.sql.execute('''CREATE TABLE drugs
                                (did text, name text, terms text)'''
                         )

        # get all the important drug ids (dids) from
        # the known gene-drug connections table

        self.sql.execute('SELECT DISTINCT did FROM drugpairs')

        # fetch matching did and use to create uri for query

        for result in tqdm(self.sql.fetchall()):

            did = result[0]

            d = Drug(did)

            for item in d.terms:

                term = item['term']

                # check for duplicates with chemical name

                if name not in term.lower():

                    item = (did, str(name), term)

                    # insert into drugs table

                    self.sql.execute('''INSERT INTO drugs VALUES(?,?,?)'''
                            , item)

    def BedFile(self):

        # creates bed file for subsetting .BAM files.

        self.sql.execute('SELECT chr, start, stop, symbol FROM genes ORDER BY length(chr), chr'
                         )

        with open('PharmacogenomicGenes_PGKB.bed', 'w') as bed:

            for tup in self.sql.fetchall():

                bed.write('\t'.join(map(str, tup)) + '\n')