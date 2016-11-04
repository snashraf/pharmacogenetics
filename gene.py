#!/usr/bin/python
# -*- coding: utf-8 -*-
import json
from lxml import etree
import urllib2
from collections import Counter
import re

class Gene:
    """
This class gets information on a given gene based on gene ID. Can use either PharmGKB or Entrez.
    """
    def __init__(self, geneid):
        # initialize with geneID
        self.geneid = geneid
        self.Load()  # loads data on gene
        self.GetDrugs()  # loads info on associated drugs
        self.GetDesc()  # gets summary of gene

    def Load(self):
        if 'PA' in self.geneid:
            uri = 'https://api.pharmgkb.org/v1/data/gene/%s?view=max' \
                % self.geneid
            data = urllib2.urlopen(uri)
            self.json = json.load(data)[0]
        elif 'PA' not in self.geneid and 'placeholder' \
            not in self.geneid:
            uri = \
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%i&retmode=XML' \
                % int(self.geneid)
            f = urllib2.urlopen(uri)
            data = f.read()
            f.close()
            self.tree = etree.XML(data)
            for elem in self.tree[0][3].iter():
                if 'ref_locus' in elem.tag:
                    self.name = elem.text

    def GetDesc(self):
        for elem in self.tree[0].iter():
            if 'summary' in elem.tag:
                self.desc = elem.text

    def GetDrugs(self):
        self.drugs = []
        for elem in self.tree[0].iter():
            if 'Gene-commentary_text' in elem.tag:
                c = Counter(re.findall(r"\w+", elem.text))
                if 'interacts with' in elem.text and c[self.name] == 1:
                    cut = elem.text.split(' ')
                    drug = cut[len(cut) - 1].strip('.')
                    self.drugs.append(drug)
        print self.drugs

    def GetHaps(self):
        # get a list of known haplotype IDs from gene ID
        uri = \
            'https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId=%s&view=max' \
            % self.geneid
        self.haps = []

