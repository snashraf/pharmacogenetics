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
        self.variants = []
        self.genes = []
        self.chemicals = []
        self.geneids = []
        self.genenames = []
        self.nohap = []
        self.data = {}
        self.changes = {}
        self.conn = sqlite3.connect('pharmacogenetics.db')
        self.sql = self.conn.cursor()

    def Update(self):
        self.GetVCFData()
        self.GetVarData()
        self.GetGeneData()
        self.GetChemData()
        self.conn.commit()

    def GetVCFData(self):
        print 'Importing VCF...'
        self.sql.execute('''DROP TABLE IF EXISTS design''')
        self.sql.execute('''CREATE TABLE design
                            (pos int, rsid int, num text, ref text, alt text,UNIQUE(pos, rsid) ON CONFLICT REPLACE)''')
        # change design source file here

        filename = self.f

        # dictionary of pos-to-rs values

        self.contab = {}
        self.reader = vcf.Reader(open(filename, 'r'))

        # create a dictionary for this

        for record in self.reader:
            self.contab[record.POS] = record.ID
            for sample in record.samples:
                call = str(sample['GT'])
                break
            for alt in record.ALT:
                item = (record.POS, record.ID, call, record.REF, str(alt))
                self.sql.execute('''INSERT INTO design VALUES(?,?,?,?,?)''', item)

    def GetVarData(self):
        print 'Getting variant data...'
        self.sql.execute('''DROP TABLE IF EXISTS variants''')
        self.sql.execute('''CREATE TABLE variants
                                    (rsid text, gid text, UNIQUE(rsid,gid) ON CONFLICT REPLACE)''')
        self.sql.execute('''DROP TABLE IF EXISTS alias''')
        self.sql.execute('''CREATE TABLE alias
                                            (rsid text, alias varchar(255) PRIMARY KEY)''')
        self.sql.execute("SELECT DISTINCT rsid FROM design")
        for result in self.sql.fetchall():
            rsid = result[0]
            print rsid
            try:
                if 'rs' in rsid:
                    v = Variant(rsid, 'pharmgkb')
                elif 'cv' in rsid:
                    v = Variant(rsid, 'clinvar')
            except urllib2.HTTPError:
                v = Variant(rsid, 'entrez')
            for tup in v.nameid:
                symbol = tup[0]
                gid = tup[1]
                self.sql.execute('''INSERT INTO variants VALUES(?,?)''', (rsid, gid))
                for alias in v.names:
                    try:
                        self.sql.execute('''INSERT INTO alias VALUES(?,?)''', (rsid, alias))
                    except sqlite3.IntegrityError:
                        pass

    def GetGeneData(self):
        print 'Getting gene data...'
        self.sql.execute('''DROP TABLE IF EXISTS genes''')
        self.sql.execute('''CREATE TABLE genes
                                            (gid text UNIQUE, symbol text)''')
        self.sql.execute('''DROP TABLE IF EXISTS alleles''')
        self.sql.execute('''CREATE TABLE alleles
                                            (hapid text, gid text, starname text, hgvs text, rsid text, UNIQUE(hapid, rsid) ON CONFLICT REPLACE)''')
        self.GetPairs()
        genes = {'pgkb': [], 'other': []}
        self.sql.execute("SELECT DISTINCT gid FROM variants")
        for result in self.sql.fetchall():
            gid = result[0]
            print gid
            if "PA" in gid:
                g = Gene(gid, 'pharmgkb')
                self.sql.execute('''INSERT INTO genes VALUES(?,?)''', (gid, g.name))
                for allele in g.alleles:
                    for rsid in allele['rsids']:
                        self.sql.execute('''INSERT INTO alleles VALUES(?,?,?,?,?)''', (allele['id'], gid, allele['starname'], allele['hgvs'], rsid))
        else:
            pass

    def GetPairs(self):
        self.sql.execute('''DROP TABLE IF EXISTS drugpairs''')
        self.sql.execute('''CREATE TABLE drugpairs
                        (did text, gid text, UNIQUE(did, gid) ON CONFLICT REPLACE)''')
        print 'Getting gene-drug pairs...'
        self.giddid = {}
        self.druglist = []
        uri = 'https://api.pharmgkb.org/v1/report/selectPairs'

            # get data and read in this json file

        data = urllib2.urlopen(uri)
        self.json = json.load(data)
        for doc in self.json:
            gid = doc['gene']['id']
            did = doc['chemical']['id']
            self.sql.execute('''INSERT INTO drugpairs VALUES(?,?)''', (did, gid))
        pairs = len(self.json)
        print '%i annotated pairs found!' % pairs

    def GetChemData(self):
        print 'Getting drug data...'
        self.sql.execute('''DROP TABLE IF EXISTS chemicals''')
        self.sql.execute('''CREATE TABLE chemicals
                                (did text, name text, terms text)''')
        self.sql.execute("SELECT DISTINCT did FROM drugpairs")
        for result in self.sql.fetchall():
            did = result[0]
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
                if name not in term.lower():
                    prev_term = term
                    item = (did, str(name), term)
                    print item
                    self.sql.execute('''INSERT INTO chemicals VALUES(?,?,?)''', item)


if __name__ == '__main__':
    data = DataCollector('config/corrected_design.vcf')
    data.Update()
else:
    pass