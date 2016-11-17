#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import urllib2
from variant import Variant
from gene import Gene
import sqlite3
import json


class DataCollector:

    '''
    This class imports a design VCF and creates a database with these positions, using PharmGKB.
    '''

    def __init__(self, f):
        """
        f = given VCF input file
        GetDesign imports the test design and matches rs values to positions
        Import imports the input vcf
        GetIDs finds changed positions and matches them to the design.
        All of this is exported
        """

        self.f = f  # filename of design
        self.conn = sqlite3.connect('pharmacogenetics.db')  # connect to db- if it doesn't exist, create it
        self.sql = self.conn.cursor()  # cursor for sqlite3, used to do things in database

    def Update(self):
        """
        This function rebuilds the database, use for updating the db periodically?
        :return:
        """

        self.GetVCFData()
        self.GetVarData()
        self.GetGeneData()
        self.GetChemData()
        self.conn.commit()

    def GetVCFData(self):
        """
        This function imports the design VCF
        :return:
        """

        print 'Importing VCF...'

        # create sql table if it doesn't exist yet, otherwise recreate it

        self.sql.execute('''DROP TABLE IF EXISTS design''')

        # this can be used to pair positions to rsids - essential for patient analysis later on
        # the combination of position and rsid should be unique.

        self.sql.execute('''CREATE TABLE design
                            (pos int, rsid int, num text, ref text, alt text,UNIQUE(pos, rsid) ON CONFLICT REPLACE)'''
                         )

        # f is the filename of the design, given when creating the datacollector object

        filename = self.f

        # create a VCF reader object, using the pyvcf module. Open file for reading.

        self.reader = vcf.Reader(open(filename, 'r'))

        # fetch data from the columns - if you want to add more data do that here.

        for record in self.reader:
            for sample in record.samples:
                call = str(sample['GT'])  # call is the 1/1 or 1|1 value showing zygosity
                break
            for alt in record.ALT:
                item = (record.POS, record.ID, call, record.REF,
                        str(alt))  # create an item to insert into sql table with relevant information
                self.sql.execute('''INSERT INTO design VALUES(?,?,?,?,?)'''
                                 , item)

    def GetVarData(self):
        """
        Using Variant objects, this function fetches data on variants from the PharmGKB servers,
        if not available uses the Entrez servers, possibilities are limited for now.
        :return:
        """

        print '--- Getting variant data ---'

        # drop tables if they exist already to reset them

        self.sql.execute('''DROP TABLE IF EXISTS variants''')
        self.sql.execute('''DROP TABLE IF EXISTS alias''')

        # create variant and alias table. Variant should have an unique combo of rsid and gid,
        # and alias should be unique in the alias table.

        self.sql.execute('''CREATE TABLE variants
                                    (rsid text, gid text,
                                    UNIQUE(rsid,gid)
                                    ON CONFLICT REPLACE)'''
                         )
        self.sql.execute('''CREATE TABLE alias
                                            (rsid text, alias varchar(255) PRIMARY KEY)'''
                         )

        # get all rsids in the design vcf

        self.sql.execute('SELECT DISTINCT rsid FROM design')

        # rotate through rsids and create variant objects to fetch information

        for result in self.sql.fetchall():
            rsid = result[0]
            print rsid

            # create variant instances with the given rsid.
            # if there is no pgkb id, try entrez instead.

            try:
                v = Variant(rsid, 'pharmgkb')
            except urllib2.HTTPError:
                v = Variant(rsid, 'entrez')

            # this results in a combination tuple of rsid and gid and aliases

            for tup in v.nameid:
                gid = tup[1]

                # use this tuple to fill sql insertion

                self.sql.execute('''INSERT INTO variants VALUES(?,?)'''
                                 , (rsid, gid))

                # go through aliases, ignore duplicates and put in alias table

                for alias in v.names:
                    try:
                        self.sql.execute('''INSERT INTO alias VALUES(?,?)'''
                                , (rsid, alias))
                    except sqlite3.IntegrityError:

                    # on duplicate, ignore

                        continue

    def GetGeneData(self):
        '''
        Fetches data on given gene IDs.
        :return:
        '''

        print '--- Getting gene data ---'

        # drop already existing tables genes and alleles

        self.sql.execute('''DROP TABLE IF EXISTS genes''')
        self.sql.execute('''DROP TABLE IF EXISTS alleles''')

        # (re)create tables in database

        self.sql.execute('''CREATE TABLE genes
                                            (gid text UNIQUE, symbol text)'''
                         )
        self.sql.execute('''CREATE TABLE alleles
                                            (hapid text, gid text, starname text,
                                            hgvs text, rsid text, alt text,
                                            UNIQUE(hapid, rsid)
                                            ON CONFLICT REPLACE)'''
                         )

        # get all unique gene ids from the variant table

        self.sql.execute('SELECT DISTINCT gid FROM variants')
        genes = [tup[0] for tup in self.sql.fetchall()]
        genes = list(set(genes))

        # go through results and create gene objects for each GID with PA (so it can be found on pharmgkb)

        for gid in genes:
            print gid
            if 'PA' in gid:

                g = Gene(gid, 'pharmgkb')

                # insert the resulting name and alleles into sql table

                self.sql.execute('''INSERT INTO genes VALUES(?,?)''',
                                 (gid, g.name))

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
        else:
            pass

    def GetPairs(self):
        """
        Create table for drug-gene pairs, to be used later to fetch drug data
        :return:
        """

        self.sql.execute('''DROP TABLE IF EXISTS drugpairs''')
        self.sql.execute('''CREATE TABLE drugpairs
                        (did text, gid text, UNIQUE(did, gid)
                        ON CONFLICT REPLACE)''' )
        print 'Getting gene-drug pairs...'

        # get the uri for finding well-annotated pairs (given by pharmgkb)

        uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

        # get data and read in this json file

        data = urllib2.urlopen(uri)
        self.json = json.load(data)
        for doc in self.json:
            gid = doc['gene']['id']
            did = doc['chemical']['id']

            # insert results in table drugpairs

            self.sql.execute('''INSERT INTO drugpairs VALUES(?,?)''',
                             (did, gid))

    def GetChemData(self):
        """
        Gets info on seperate drugs, such as name and associated terms
        terms give info on what the drug does
        :return:
        """

        print '--- Getting drug data ---'

        # drop and re-create table chemicals

        self.sql.execute('''DROP TABLE IF EXISTS chemicals''')
        self.sql.execute('''CREATE TABLE chemicals
                                (did text, name text, terms text)'''
                         )

        # get all the important drug ids (dids) from
        # the known gene-drug connections table

        self.sql.execute('SELECT DISTINCT did FROM drugpairs')

        # fetch matching did and use to create uri for query

        for result in self.sql.fetchall():
            did = result[0]
            print did
            uri = \
                'https://api.pharmgkb.org/v1/data/chemical/%s?view=max' \
                % did

            # get data and read in this json file

            data = urllib2.urlopen(uri)
            self.json = json.load(data)
            name = self.json['name']
            terms = self.json['terms']
            for item in terms:
                term = item['term']

                # check for duplicates with chemical name

                if name not in term.lower():

                    item = (did, str(name), term)

                    # insert into chemicals table

                    self.sql.execute('''INSERT INTO chemicals VALUES(?,?,?)'''
                            , item)


if __name__ == '__main__':
    data = DataCollector('config/corrected_design.vcf')
    data.Update()
else:
    pass