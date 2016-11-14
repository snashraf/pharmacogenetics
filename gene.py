#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
from lxml import etree
import urllib2
from collections import Counter
import re


# -------------------------------------------------------------------------

class Gene:

    """
This class gets information on a given gene based on gene ID. Can use either PharmGKB or Entrez.
    """

    def __init__(self, geneid, mode):

        # initialize with geneID

        self.gid = geneid
        self.mode = mode
        self.alleles = []
        self.genes = []
        self.Load()  # loads data on gene

        # self.GetDrugs()  # loads info on associated drugs
        # self.GetDesc()  # gets summary of gene

        self.GetHaps()

    def Load(self):
        if self.mode == 'pharmgkb':
            uri = 'https://api.pharmgkb.org/v1/data/gene/%s?view=max' \
                % self.gid
            data = urllib2.urlopen(uri)
            self.json = json.load(data)
            self.name = self.json["symbol"]
        elif self.mode == 'entrez':
            uri = \
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%i&retmode=XML' \
                % int(self.gid)
            f = urllib2.urlopen(uri)
            data = f.read()
            f.close()
            self.tree = etree.XML(data)
            for elem in self.tree[0][3].iter():
                if 'ref_locus' in elem.tag:
                    self.name = elem.text

    def GetDesc(self):
        if self.mode == 'pharmgkb':
            print self.json
        elif self.mode == 'entrez':
            for elem in self.tree[0].iter():
                if 'summary' in elem.tag:
                    self.desc = elem.text

    def GetHaps(self):
        if self.mode == 'entrez':
            return
        else:

        # get a list of known haplotype IDs from gene ID

            uri = \
                'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId=%s&view=max' \
                % self.gid
            try:
                data = urllib2.urlopen(uri)
                self.json = json.load(data)
                self.haps = []
                for doc in self.json:
                    starname = doc['name']
                    hgvs = doc['hgvs']
                    alleles = doc['alleles']
                    hapid = doc['id']
                    copynum = doc['copyNumber']
                    rsids = []
                    for allele in alleles:
                        try:
                            change = allele['allele']
                            rsid = allele['location']['displayName']
                            rsids.append(rsid)
                        except:
                            continue
                    all = {
                        'starname': starname,
                        'hgvs': hgvs,
                        'id': hapid,
                        'copynum': copynum,
                        'rsids': rsids,
                        }
                    self.alleles.append(all)
            except:
                return