#!/usr/bin/python
# -*- coding: utf-8 -*-
import vcf
import pickle
from data import DataCollector

#-------------------------------------------------------------------------

def Find(dict, k, v):
    match = (item for item in dict if item[k] == v).next()
    return match

#-------------------------------------------------------------------------


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
        self.patientvars = []
        self.patientgenes = {}
        self.reader = vcf.Reader(open(f, 'r'))
        self.data = DataCollector('config/design.vcf')
        self.data.GetVCFData()
        self.contab = self.data.contab
        self.giddrugs = {}
        self.rsdrugs= []
        self.GetIDs()
        self.ImportData()
        self.GetGenes()
        self.GetChems()

    def ImportData(self):
        try:
            with open("variants.pickle", "rb") as f:
                self.variants = pickle.load(f)
            with open("genes.pickle", "rb") as f:
                self.genes = pickle.load(f)
            with open("chemicals.pickle", "rb") as f:
                self.chemicals = pickle.load(f)
        except:
            "Redownloading data..."
            print "Collecting data on variants..."
            data.GetVarData()
            print "Collecting data on gene haplotypes..."
            data.GetGeneData()
            print "Getting data on connected chemicals..."
            data.GetChemData()
            self.ImportData()

    def GetIDs(self):
        # create list for storing changed positions / rs#s
        for record in self.reader:
            # check for non-None positions, these have changed
            if "[None]" not in str(record.ALT):
                ref = record.REF
                alt = record.ALT
                for sample in record.samples:
                    call = str(sample['GT'])
                try:
                    # match position to rs # and add rs to storage
                    rs = self.contab[record.POS]
                    self.patientvars.append({"pos":record.POS,"rsid":rs,"call":{"num":call,"ref":ref,"alt":alt}})
                except KeyError:
                    # if no match is found, let user know
                    print "couldn't match", record.POS
                    pass

    def GetGenes(self):
        for var in self.patientvars:
            rsid = var['rsid']
            match = Find(self.variants, "rsid", rsid)
            gene = match["gene"]["id"]
            if gene in self.patientgenes.keys():
                self.patientgenes[gene].append(rsid)
            elif gene not in self.patientgenes.keys():
                self.patientgenes[gene] = [rsid]

    def GetChems(self):
        for gid, values in self.patientgenes.items():
            match = ""
            try:
                match = Find(self.genes, "id", gid)
                match = match['drugs']
            except:
                print "Gene %s not found!" %gid
            if len(match) > 0:
                self.giddrugs[gid] = match
                self.rsdrugs.append((self.patientgenes[gid], match))