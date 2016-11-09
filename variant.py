#!/usr/bin/python
# -*- coding: utf-8 -*-
import json
from lxml import etree
import urllib2
from wikitools import wiki, page

# -------------------------------------------------------------------------

class Variant:
    """
    This class stores information on variants. It can take info from pharmkgb,
     entrez or snpedia (for additional information))
    """
    def __init__(self, rs, mode):
        # initialize
        self.mode = mode  # mode of variant search
        self.rs = rs  # rs number of variant
        self.Load()  # load data (either json, wikitools orxml)
        self.GetGene()  # get gene matching the rs#
        self.GetAlias()  # get HGVS alias'

    def Load(self):
        # check for which mode is being used and adjust query accordingly
        if self.mode == 'pharmgkb':
            uri = \
                'https://api.pharmgkb.org/v1/data/variant/?symbol=%s&view=max' \
                % self.rs
            # get data and read in this json file
            data = urllib2.urlopen(uri)
            self.json = json.load(data)[0]
        elif self.mode == 'entrez':
            uri = \
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=%i&report=XML' \
                % int(self.rs.lstrip('rs'))
            f = urllib2.urlopen(uri)
            data = f.read()
            f.close()
            # create xml tree from downloaded file
            self.tree = etree.XML(data)

    def GetGene(self):
        # again check for which mode is used'
        if self.mode == 'pharmgkb':
            # get gene name and PHARMGKB gene Id
            self.genename = self.json['relatedGenes'][0]['symbol']
            self.geneid = self.json['relatedGenes'][0]['id']
        elif self.mode == 'entrez':
            # get gene name and ENTREZ gene id (only used if pharmgkb fails)
            for elem in self.tree.iter('*'):
                if 'symbol' in elem.attrib:
                    self.genename = elem.attrib['symbol']
                    self.geneid = elem.attrib['geneId']

    def GetAlias(self):
        # get HGVS alias' for variant, depending on mode again (use later)
        if self.mode == 'pharmgkb':
            self.names = self.json['altNames']['synonym']
        elif self.mode == 'entrez':
            self.names = []
            for elem in self.tree.iter():
                if 'hgvs' in elem.tag:
                    self.names.append(unicode(elem.text))