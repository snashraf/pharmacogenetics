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
from operator import itemgetter
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
        self.conn.commit()
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
        self.sql.execute("SELECT pos, rsid from variants AS V JOIN patientvars as P ON V.rsid=P.rsid")
        for result in self.sql.fetchall():
            pos = result[0]
            self.sql.execute("SELECT rsid FROM design WHERE pos = ?", (pos,))
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
        self.sql.execute("SELECT DISTINCT gid FROM variants, patientvars WHERE pos=pos")
        genes = self.sql.fetchall()
        patienthaps = {}
        self.nohaps = []
        for result in genes:
            gid = result[0]
            if "PA" not in gid:
                self.nohaps.append(gid)
                continue
            self.sql.execute("SELECT symbol FROM genes WHERE gid = ?", (gid,))
            patienthaps[gid] = []
            self.sql.execute("SELECT DISTINCT rsid,allele FROM patientvars,variants WHERE pos=pos AND gid = ?", (gid,))
            patientrsids = self.sql.fetchall()
            self.sql.execute("SELECT DISTINCT hapid, starname FROM alleles WHERE gid = ?", (gid,))
            haps = self.sql.fetchall()
            for hap in haps:
                hapid = hap[0]
                starname = hap[1]
                if "*1" not in starname:
                    # get rsids for haplotype
                    self.sql.execute("SELECT DISTINCT rsid FROM alleles WHERE hapid = ?", (hapid,))
                    # USE HAPLOTYPE FUNCTION
                    haprsids = self.sql.fetchall()
                    # compare patient with rsidlist
                    comparison = set(patientrsids) & set(haprsids)
                    match_score = float(len(comparison))/len(haprsids)
                    if match_score == 1.0:
                        result = (hapid, starname, len(haprsids))
                        if result not in patienthaps[gid]:
                            patienthaps[gid] . append(result)
                elif "*1" in starname:
                    reference = hapid
                else:
                    self.nohaps.append(gid)
        self.haplotypes = []

        for gid, result in patienthaps.items():
            if len(result) > 0:
                highest_hit = max(result,key=itemgetter(2))
                hapid = highest_hit[0]
                self.haplotypes.append(hapid)
        self.nohaps = list(set(self.nohaps))

    def Profiler(self, mode):
        self.authobj = Authenticate()
        self.Hapmatcher()
        output = []
        if mode == "rsid":
            self.sql.execute("SELECT rsid, allele FROM patientvars")
            output = self.sql.fetchall()
        elif mode == "haplotype":
            output += self.haplotypes
            for gid in self.nohaps:
                self.sql.execute("SELECT DISTINCT rsid FROM variants WHERE gid = ?", (gid,))
                for result in self.sql.fetchall():
                    rsid = result[0]
                    self.sql.execute("SELECT pos FROM design WHERE rsid = ?", (rsid,))
                    pos = self.sql.fetchone()[0]
                    self.sql.execute("SELECT allele FROM patientvars WHERE pos = ?", (pos,))
                    allele = self.sql.fetchall()
                    print allele
                    output += rsid
        print mode, output
        #self.Drugadvice(output)

    def DrugAdvice(self, input):
            for var in input:
                if type(var) == tuple:
                    rsid = result[0]
                    allele = result[1]
                if mode == "haplotype":
                    self.sql.execute("SELECT DISTINCT gid FROM alleles WHERE hapid = ?", (var,))
                elif mode == "rsid":
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

    def PerDrug(self):
        pass

tom = Patient('data/hpc/test.vcf')
tom.Profiler("haplotype")

#print tom.adviceperdrug

