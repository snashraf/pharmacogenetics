#!/usr/bin/python
# -*- coding: utf-8 -*-

import json
from lxml import etree
import urllib2


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
            return
        elif self.mode == 'entrez':
            uri = \
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=%i&report=XML' \
                % int(self.rs.lstrip('rs'))
            f = urllib2.urlopen(uri)
            data = f.read()
            f.close()

            # create xml tree from downloaded file

            self.tree = etree.XML(data)
            return
        elif self.mode == 'clinvar':
            cvid = self.rs.split('cv')[1]
            uri = \
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=%s&retmode=json' \
                % cvid
            data = urllib2.urlopen(uri)
            self.json = json.load(data)
            return

    def GetGene(self):

        # again check for which mode is used'

        if self.mode == 'pharmgkb':

            # get gene name and PHARMGKB gene Id

            self.nameid = [(doc['symbol'], doc['id']) for doc in
                           self.json['relatedGenes']]
        elif self.mode == 'entrez':

            # get gene name and ENTREZ gene id (only used if pharmgkb fails)

            for elem in self.tree.iter('*'):
                if 'symbol' in elem.attrib:
                    name = elem.attrib['symbol']
                    gid = elem.attrib['geneId']
                    self.nameid = [(name, gid)]
        elif self.mode == 'clinvar':

            self.nameid = [('BCHE', 'PA25294')]

    def GetLocation(self):
        if self.mode == 'pharmgkb':
            try:
                if 'GRCh37' in self.json['location']['name']:
                    self.chrom = self.json['location']['name'].split(']'
                            )[1].split(':')[0].strip('chr')
                    self.begin = self.json['location']['begin']
                    self.end = self.json['location']['end']
                    self.ref = self.json['location']['reference']
                    self.alt = ','.join(self.json['location']['variants'
                            ])
            except:
                print 'cant find location for', self.rs
        elif self.mode == 'entrez':

            return

    def GetAlias(self):

        # get HGVS alias' for variant, depending on mode again (use later)

        if self.mode == 'pharmgkb':
            try:
                self.names = self.json['altNames']['synonym']
            except:
                self.names = []
                for doc in self.json['alternateLocations']:
                    if 'RefSeq DNA' in doc['sequence']['resource']:
                        xref = doc['sequence']['xrefId']
                        pos = doc['begin']
                        ref = doc['reference']
                        alt = doc['variants'][0]
                        name = xref + ':g.' + str(pos) + ref + '>' + alt
                        self.names.append(name)
        elif self.mode == 'entrez':
            self.names = []
            for elem in self.tree.iter():
                if 'hgvs' in elem.tag:
                    self.names.append(unicode(elem.text))
            return
        elif self.mode == 'clinvar':
            self.names = []
            return