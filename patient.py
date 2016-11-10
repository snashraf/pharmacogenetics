#!/usr/bin/python
# -*- coding: utf-8 -*-

import vcf
import pickle
from data import DataCollector


# -------------------------------------------------------------------------

def Find(lst, k, v):
    matches = []
    for dct in lst:
        matches += [item for item in lst if item[k] == v]
    if len(matches) > 0:
        return matches
    else:
        return None


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
        self.patientvars = []
        self.patientgenes = {}
        self.giddrugs = {}
        self.rsdrugs = []
        self.alleles = {}
        print 'Initiating patient...'

    def Load(self):
        self.d = DataCollector('config/corrected_design.vcf')
        self.d.GetVCFData()
        self.GetIDs()
        self.ImportData()
        self.GetGenes()
        self.GetChems()
        self.GetAlleles()

        # self.Hapmaker()

    def ImportData(self):
        print 'Trying to import database...'
        try:
            with open('variants.pickle', 'rb') as f:
                self.variants = pickle.load(f)
            with open('genes.pickle', 'rb') as f:
                self.genes = pickle.load(f)
            with open('chemicals.pickle', 'rb') as f:
                self.chemicals = pickle.load(f)
        except:
            self.Update()

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

                    rs = self.d.contab[record.POS]
                    print rs
                    self.patientvars.append({'pos': record.POS,
                            'rsid': rs, 'call': {'num': call,
                            'ref': ref, 'alt': alt}})
                except KeyError:

                    # if no match is found, let user know

                    print "couldn't match", record.POS
                    pass

    def GetGenes(self):
        print 'Finding associated genes...'
        for var in self.patientvars:
            rsid = var['rsid']
            matches = Find(self.variants, 'rsid', rsid)
            for match in matches:
                for gene in match['genes']:
                    gid = gene['id']
                    if gid in self.patientgenes.keys():
                        if rsid not in self.patientgenes[gid]:
                            self.patientgenes[gid].append(rsid)
                        else:
                            pass
                    elif gid not in self.patientgenes.keys():
                        self.patientgenes[gid] = [rsid]

    def GetChems(self):
        print 'Finding matching drugs...'
        for (gid, values) in self.patientgenes.items():
            match = Find(self.genes, 'id', gid)
            if match is not None:
                match = match[0]['drugs']
            else:
                continue
            if len(match) > 0:
                self.giddrugs[gid] = match
                self.rsdrugs.append((self.patientgenes[gid], match))

    def GetAlleles(self):
        print 'Determining patient alleles...'
        for var in self.patientvars:
            allele = ''
            rsid = var['rsid']
            call = var['call']['num']
            ref = var['call']['ref']
            alt = var['call']['alt']
            if '/' in call:
                if call == '0/1' or call == '1/0':
                    allele = ref + str(alt[0])
                elif call == '1/1':
                    allele = str(alt[0]) * 2
            elif '|' in call:
                pass
            self.alleles[rsid] = allele

    def Hapmaker(self):
        hapdictionary = {}
        for var in self.patientvars:
            filt_names = []
            rsid = var['rsid']
            alt = var['call']['alt']
            matches = Find(self.variants, 'rsid', rsid)
            names = matches[0]['alias']
            gen_names = names['genomic']
            for name in gen_names:
                if 'NC' in name:
                    for nucl in alt:
                        if '>' + str(nucl) in name:
                            filt_names.append(name)
            hapdictionary[rsid] = filt_names
        for (gene, rsids) in self.patientgenes.items():
            print 'Haplotype for %s:' % gene
            for rsid in rsids:
                print rsid
                print hapdictionary[rsid]

    def Hapmatcher(self):
        pass


tom = Patient('data/hpc/test.vcf')
tom.Load()
tom.Hapmaker()