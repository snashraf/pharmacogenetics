#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import sqlite3
import json
import requests
from requests_oauthlib import OAuth2Session
import urllib
import urllib2
import ast
from data import DataCollector


# -------------------------------------------------------------------------

def Authenticate():
    print 'Authenticating...'
    req = ''
    with open('config/auth.txt', 'r') as f:
        req = json.load(f)
    url = 'https://api.pharmgkb.org/v1/auth/oauthSignIn'
    data = urllib.urlencode(req)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    str_response = response.read()
    token = ast.literal_eval(str_response)
    client = OAuth2Session(token=token)
    return client

def getJson(uri, client):
    r = client.get(uri)
    data = r.json()
    if type(data) == dict:
        return None
    elif type(data) == list:
        return data

def PGKB_connect(authobj, mode, a, did):
    uri = \
        'https://api.pharmgkb.org/v1/report/pair/%s/%s/clinicalAnnotation' \
        % (a, did)
    result = getJson(uri, authobj)
    results = []
    if result is not None:
        for doc in result:
            for phen in doc['allelePhenotypes']:
                if mode == "rsid":
                    results.append(phen)
                elif mode == "haplotype":
                    pass
    else:
        return
    return results

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
        self.vars = []
        self.genes = {}
        self.giddrugs = {}
        self.alleles = {}
        self.adviceperdrug = []
        self.Load()
        print 'Initiating patient...'

    def Load(self):
        self.d = DataCollector('config/corrected_design.vcf')
        self.ImportData()
        self.GetIDs()
        # self.Hapmaker()

    def ImportData(self):
        print 'Trying to import database...'
        try:
            self.conn = sqlite3.connect('pharmacogenetics.db')
            self.sql = self.conn.cursor()
            self.sql.execute('''DROP TABLE IF EXISTS patientvars''')
            self.sql.execute('''CREATE TABLE patientvars
                                                (pos int PRIMARY KEY, num text, ref text, alt text, allele text)''')
        except:
            raise
            #self.d.Update()

    def GetIDs(self):
        print 'Reading patient variants...'

        # create list for storing changed positions / rs#s

        for record in self.reader:

            # check for non-None positions, these have changed

            if '[None]' not in str(record.ALT):
                ref = record.REF
                alt = record.ALT
                for sample in record.samples:
                    call = str(sample['GT'])
                try:
                    # match position to rs # and add rs to storage
                    if '/' in call:
                        if call == '0/1' or call == '1/0':
                            allele = ref + str(alt[0])
                        elif call == '1/1':
                            allele = str(alt[0]) * 2
                    elif '|' in call:
                        pass
                    item = (record.POS, call, record.REF, str(alt[0]), allele)
                    self.sql.execute('''INSERT INTO patientvars VALUES(?,?,?,?,?)''', item)

                except KeyError:

                    # if no match is found, let user know

                    print "couldn't match", record.POS
                    pass

    def Hapmaker(self):
        self.sql.execute("SELECT pos, rsid FROM variants AS V JOIN patientvars as P ON V.rsid=P.rsid")
        for result in self.sql.fetchall():
            print result
            pos = result[0]
            self.sql.execute("SELECT rsid FROM design WHERE pos = ?", (pos,))
            rsid = var['rsid']
            alt = var['call']['alt']
            self.sql.execute("SELECT alias FROM alias WHERE rsid = ?", (rsid,))
            for result in self.sql.fetchall():
                alias = result[0]
                if 'NC' in alias and "g." in alias:
                    for nucl in alt:
                        if '>' + str(nucl) in alias:
                            pass
        for (gene, rsids) in self.genes.items():
            print 'Haplotype for %s:' % gene
            for rsid in rsids:
                pass

    def Hapmatcher(self):
        genes = []
        self.sql.execute("SELECT pos, allele FROM patientvars")
        self.sql.execute("SELECT rsid, hapid, allele FROM alleles AS A JOIN patientvars AS P ON V.pos=P.pos")
        for result in self.sql.fetchall():
            print result
            rsid = result[0]
            self.sql.execute("SELECT DISTINCT gid FROM variants WHERE rsid = ?", (rsid,))
            for result in self.sql.fetchall():
                gid = result[0]
                genes.append(gid)
        for gid in genes:
            self.sql.execute("SELECT DISTINCT hapid FROM alleles WHERE gid = ?", (gid,))
            haplotypes = self.sql.fetchall()
            if len(haplotypes) > 0:
                for hap in haplotypes:
                    if "*1" not in hap[0]:
                        print hap


    def DrugAdvice(self, mode):
        authobj = Authenticate()
        if mode == "rsid":
            self.sql.execute("SELECT rsid, allele FROM patientvars")
            for result in self.sql.fetchall():
                rsid = result[0]
                allele = result[1]
                self.sql.execute("SELECT DISTINCT gid FROM variants WHERE rsid = ?", (rsid,))
                for result in self.sql.fetchall():
                    gid = result[0]
                    self.sql.execute("SELECT DISTINCT did FROM drugpairs WHERE gid = ?", (gid, ))
                    for result in self.sql.fetchall():
                        did = result[0]
                        results = PGKB_connect(authobj, mode, rsid, did)
                        if results is not None:
                            for phen in results:
                                if allele in phen['allele']:
                                    entry = {
                                        'did': did,
                                        'rsid': rsid,
                                        'phenotype': phen['phenotype'],
                                    }
                                    self.adviceperdrug.append(entry)
        elif mode == "haplotype":
            self.Hapmatcher()
            pass
    def PerDrug(self):
        pass

tom = Patient('data/hpc/test.vcf')
#tom.DrugAdvice("haplotype")
#print tom.adviceperdrug

