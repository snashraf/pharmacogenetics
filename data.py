#!/usr/bin/python
# -*- coding: utf-8 -*-
import vcf
import urllib2
import pprint
from variant import Variant
from gene import Gene
import numpy as np
import pickle
import pymongo
import json

# -------------------------------------------------------------------------

def Convert(dict):
    n = {}
    for k, v in dict.items():
        if isinstance(k, unicode):
            for i in ['utf-8', 'iso-8859-1']:
                try:
                    k = k.encode(i)
                except (UnicodeEncodeError, UnicodeDecodeError):
                    continue
        if isinstance(v, np.int64):
            print("k is %s , v is %s" % (k, v))
            v = int(v)
            print("V is %s" % v)
        if isinstance(v, unicode):
            for i in ['utf-8', 'iso-8859-1']:
                try:
                    v = v.encode(i)
                except (UnicodeEncodeError, UnicodeDecodeError):
                    continue
        n[k] = v
        print n
    return n

def Pickle( dict, name):
    print "Exporting pickle file..."
    with open(name+".pickle", "wb") as file:
        pickle.dump(dict, file)

def toMongo(collection, dict):
    print "Exporting to MongoDB..."
    for item in dict:
        try:
            collection.insert_one(item)
        except pymongo.errors.InvalidDocument:
            n = Convert(item)
            collection.insert(n)

# --------------------------------------------------------------------------------------

class DataCollector:
    '''
    This class imports a design VCF and an input VCF.
    '''
    def __init__(self,f):
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
        self.genenames=[]
        self.nohap = []
        self.data = {}
        self.changes = {}

    def Update(self):
        self.GetVCFData()
        self.GetVarData()
        self.GetGeneData()
        self.GetChemData()

    def GetVCFData(self):
        print "Importing VCF..."
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
            self.variants.append({"hg19pos":record.POS, "rsid":record.ID, "call":{"num": call, "ref":record.REF, "alt":record.ALT}})

    def GetVarData(self):
        print "Getting variant data..."
        for doc in self.variants:
            print doc["rsid"]
            try:
                if "rs" in doc["rsid"]:
                    v = Variant(doc["rsid"], 'pharmgkb')
                elif "cv" in doc["rsid"]:
                    v = Variant(doc["rsid"], "clinvar")
            except urllib2.HTTPError:
                v = Variant(doc["rsid"], 'entrez')
            doc["genes"] = []
            for tup in v.nameid:
                pair={"symbol":tup[0], "id":tup[1]}
                doc["genes"].append(pair)
            doc["alias"] = {"genomic":[],"coding":[],"protein":[],"nucleotide":[], "rsid":[]}
            for name in v.names:
                if ":g." in name:
                    doc["alias"]["genomic"].append(name)
                elif ":c." in name:
                    doc["alias"]["coding"].append(name)
                elif ":p." in name:
                    doc["alias"]["protein"].append(name)
                elif ":n." in name:
                    doc["alias"]["nucleotide"].append(name)
                elif "rs" in name:
                    doc["alias"]["rsid"].append(name)
                else:
                    pass
        Pickle(self.variants, "variants")

    def GetGeneData(self):
        self.GetPairs()
        print "Getting gene data..."
        genes = {"pgkb":[], "other":[]}
        gid_sym = {}
        for doc in self.variants:
            genelist = doc["genes"]
            for subdoc in genelist:
                sym = subdoc["symbol"]
                gid = subdoc["id"]
                gid_sym[gid] = sym
            if 'PA' in gid and gid not in genes["pgkb"]:
                genes["pgkb"].append(gid)
            elif 'PA' not in gid and gid not in genes["other"]:
                genes["other"].append(gid)
        for gid in genes["pgkb"]:
            print gid
            g = Gene(gid, 'pharmgkb')
            try:
                drugs = self.giddid[gid]
            except:
                drugs = []
            self.genes.append({"symbol":gid_sym[gid],"id":gid, "alleles":g.alleles, "drugs":drugs})
        Pickle(self.genes, "genes")

    def GetPairs(self):
        print "Getting gene-drug pairs..."
        self.giddid = {}
        self.druglist = []
        uri = \
                'https://api.pharmgkb.org/v1/report/selectPairs'
            # get data and read in this json file
        data = urllib2.urlopen(uri)
        self.json = json.load(data)
        for doc in self.json:
            gid = doc["gene"]["id"]
            did = doc["chemical"]["id"]
            self.druglist.append(did)
            if gid in self.giddid.keys():
                self.giddid[gid].append(did)
            else:
                self.giddid[gid] = [did]
        pairs = len(self.json)
        print "%i annotated pairs found!" %pairs

    def GetChemData(self):
        print "Getting drug data..."
        for did in self.druglist:
            print did
            uri = \
                'https://api.pharmgkb.org/v1/data/chemical/%s?view=max' %did
            # get data and read in this json file
            data = urllib2.urlopen(uri)
            self.json = json.load(data)
            name = self.json["name"]
            terms = self.json["terms"]
            self.chemicals.append({"id":did, "name":name, "terms":terms})
        Pickle(self.chemicals, "chemicals")

if __name__ == "__main__":
    data = DataCollector('config/corrected_design.vcf')
    data.Update()
else:
    pass