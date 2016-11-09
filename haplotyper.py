#!/usr/bin/python
# -*- coding: utf-8 -*-
import urllib2
from data import Data
from variant import Variant
from gene import Gene
import re
import json

class Haplotyper:
    def __init__(self, vcffile):
        print "Initiating..."
        self.d = Data(vcffile)
        print "Loading variants..."
        self.GetVarData()
        #print "Loading genes..."
        #self.GetGeneData()



    def HapMaker(self, mode):
        self.haplotypes = []
        if mode == 'hgvs':
            self.NameCollector()
            for gene, rss in self.rstoname.items():
                for item in rss:
                    rs = item[0]
                    names = item[1]
                    print rs, names
                    vers = []
                    posMap = {}
                    positions = []
                    for name in names:
                        split_name = name.split(".")
                        ver = split_name[1].strip(":g")
                        vers.append(ver)
                for ver in vers:
                    result = vers.count(ver)
                    print result
                    print len(names)
                    if result == len(names):
                        finalVer = ver
                    else:
                        print "not all have same ver available"
                for rs in rss:
                    for name in names:
                        split_name = name.split(".")
                        preamble = ".".join(split_name[0:2])
                        afteramble =split_name[2]
                        splitafter = re.split("([ACTG]+)", afteramble)
                        pos = int(splitafter[0])
                        val = "".join(splitafter[1::])
                        posMap[pos] = val
                        ver = split_name[1].strip(":g")
                        if ver == finalVer:
                            positions.append(pos)
                positions.sort()
                haplotype = preamble + "[" + ";".join([str(pos)+posMap[pos] for pos in positions]) + "]"
                self.haplotypes.append((gene, haplotype))
        elif mode == 'rs':
            self.GetKnownHaps()
            for gene, snps in self.data.items():
                print gene, snps

    def GetKnownHaps(self):
        pharmgenes = []
        for _ in self.geneids:
             if "PA" in _:
                 pharmgenes.append(_)
        for gene in pharmgenes:
            uri = "https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId=%s&view=max" %gene
            data = urllib2.urlopen(uri)
            storage = json.load(data)[0]
            for hap in storage:
                print hap

    def StarHap(self):
        pass
