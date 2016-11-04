#!/usr/bin/python
# -*- coding: utf-8 -*-
import urllib2
from data import Data
from variant import Variant
from gene import Gene
import re

class Haplotyper(Data):
    def __init__(self, vcffile):
        print "Initiating..."
        Data.__init__(self, vcffile)
        print "Loading variants..."
        self.GetVarData()
        #print "Loading genes..."
        #self.GetGeneData()

    def GetVarData(self):
        self.geneids = []
        self.genenames=[]
        self.data = {}
        self.hgvs = {}
        self.changes = {}
        for rs in self.rss:
            self.changes[rs[0]] = rs[1]
            try:
                v = Variant(rs[0], 'pharmgkb')
            except urllib2.HTTPError:
                    v = Variant(rs[0], 'entrez')
            self.hgvs[rs[0]] = v.names
            if v.genename in self.data.keys():
                self.data[v.genename].append(rs[0])
            elif v.genename not in self.data.keys():
                self.data[v.genename] = [rs[0]]
            if v.geneid not in self.geneids:
                self.geneids.append(v.geneid)
                self.genenames.append(v.genename)

    def HapCollector(self):
        self.hgvshaps = {}
        for name in self.genenames:
            hap = []
            for rs in self.data[name]:
                filtNames = []
                alt = self.changes[rs]
                baseChange = ">%s" %alt[0]
                for alias in self.hgvs[rs]:
                    if "NC" in alias and baseChange in alias:
                        filtNames.append(alias)
                hap.append(filtNames)
            self.hgvshaps[name] = hap

    def HapMaker(self, mode):
        self.haplotypes = []
        if mode == 'hgvs':
            self.HapCollector()
            for gene,snps in self.hgvshaps.items():
                vers = []
                posMap = {}
                positions = []
                for snp in snps:
                    for name in snp:
                        splitName = name.split(".")
                        ver = splitName[1].strip(":g")
                        vers.append(ver)
                    for ver in vers:
                        result = vers.count(ver)
                        if result == len(snps):
                            finalVer = ver
                for snp in snps:
                    for name in snp:
                        splitName = name.split(".")
                        preamble = ".".join(splitName[0:2])
                        afteramble =splitName[2]
                        splitafter = re.split("([ACTG]+)", afteramble)
                        pos = int(splitafter[0])
                        val = "".join(splitafter[1::])
                        posMap[pos] = val
                        ver = splitName[1].strip(":g")
                        if ver == finalVer:
                            positions.append(pos)
                positions.sort()
                haplotype = preamble + "[" + ";".join([str(pos)+posMap[pos] for pos in positions]) + "]"
                self.haplotypes.append((gene, haplotype))
        elif mode == 'rs':
            for gene, snps in self.data.items():
                print gene, snps

    def GetKnownHaps(self):
        pass

    def StarHap(self):
        pass
