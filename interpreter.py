#!/usr/bin/python
# -*- coding: utf-8 -*-
from patient import Patient
from data import DataCollector
import json
import sqlite3
from pgkb_functions import *

class Interpreter:


    def __init__(self, patientObj):

        self.p = patientObj

        self.advice = {}

        self.conn = sqlite3.connect('pharmacogenetics.db')

        self.sql = self.conn.cursor()

        #self.DrugAnnotations()

        self.conn.commit()

    def DrugAnnotations(self):

        self.sql.execute("DROP TABLE IF EXISTS patientadvice")

        self.sql.execute('''CREATE TABLE patientadvice
                                            (did text, gid text, varid text, allele text, advice text, dose text)''')
        authobj = Authenticate()

        self.hapdrug = {}

        self.sql.execute('''SELECT v.rsid, p.ref, p.alt, p.call
                            FROM patientvars p
                           JOIN variants v ON p.start = v.start''')
        var_alleles={}
        
        for (rsid, ref, alt, call) in self.sql.fetchall():
        
            if call == "0/0":
        
                allele = ref+ref
        
            elif call == "0/1":
        
                allele = ref+alt
        
            elif call == "1/1":
        
                allele = alt+alt

            var_alleles[rsid]=allele

        for var in var_alleles.keys():
        
            # check for haplotype or rs nomenclature

            self.sql.execute("SELECT DISTINCT gid FROM variants WHERE rsid = ?", (var,))

            gids = [tup[0] for tup in self.sql.fetchall()]

            # find gene ids associated with this variant or haplotype
        
            # print gids

            for gid in gids:

                # for each gene, find associated drugs

                self.sql.execute("SELECT DISTINCT did FROM drugpairs WHERE gid = ?", (gid,))

                dids = [tup[0] for tup in self.sql.fetchall()]

                for did in dids:

                    allele = var_alleles[var]
              
                    self.sql.execute("SELECT varid from variants where rsid=?", (var,))
              
                    varid = self.sql.fetchone()[0]

                    results = PGKB_connect(authobj, "clinAnno", varid, did)

                    if results is not None:

                        print "-----------", allele, "-----------"

                        for phen in results:

                            #print allele

                            if allele in phen['allele']:

                                advice = phen['phenotype']

                                print advice

                                item = (did, gid, var, allele, advice, "N/A") # tuple with 5 items

                                self.sql.execute('''INSERT INTO patientadvice VALUES(?,?,?,?,?,?)''', item)

    def DoseGuideCheck(self):
        authobj = Authenticate()
        doseInfo = None
        doseGuide = None
        # input: haplotypes and drug matches
        self.sql.execute("DROP TABLE IF EXISTS guidelines")
        self.sql.execute("CREATE TABLE IF NOT EXISTS druginfo(hapid text, did text, source text, guide text, specific text)")
        # first do variant annotations
        self.sql.execute("select distinct v.varid, d.did, p.ref, p.alt, p.call, p.start from variants v join patientvars p on p.start = v.start join drugpairs d on d.gid = v.gid where d.rsids like ('%' || v.rsid || '%') group by d.did")
        for (varid, did, ref, alt, call, start) in self.sql.fetchall():
            if call == "0/0":
                allele = ref+ref
            elif call == "0/1":
                allele = ref+alt
            elif call == "1/1":
                allele = alt+alt
            print allele
            results = PGKB_connect(authobj, "clinAnno",varid, did)
            if results is not None:
                for phen in results:
                    for al in phen['allelePhenotypes']:
                        if allele in al['allele']:
                            advice = al['phenotype']
                            print advice
                        else:
                            print "no match", al['allele']
        """
        self.sql.execute("SELECT DISTINCT gid, hapid FROM alleles")
        results = self.sql.fetchall()
        sources = ["cpic", "dpwg", "pro"]
        """
        self.conn.commit()

pat = Patient('data/test.g.vcf.gz')

print "patient created"

app = Interpreter(pat)
app.DoseGuideCheck()
