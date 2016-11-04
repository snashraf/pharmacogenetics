#!/usr/bin/python
# -*- coding: utf-8 -*-
import vcf


class Data:
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
        self.GetDesign()
        self.Import()
        self.GetIDs()

    def Import(self):
        # import file with the pyvcf module"
        self.reader = vcf.Reader(open(self.f, 'r'))

    def GetDesign(self):
        # change design source file here
        filename = "config/design.vcf"
        # dictionary of pos-to-rs values
        self.contab = {}
        self.reader = vcf.Reader(open(filename, 'r'))
        # create a dictionary for this
        for record in self.reader:
            self.contab[record.POS] = record.ID

    def GetIDs(self):
        # create list for storing changed positions / rs#s
        self.rss = []
        for record in self.reader:
            # check for non-None positions, these have changed
            if "[None]" not in str(record.ALT):
                try:
                    # match position to rs # and add rs to storage
                    rs = self.contab[record.POS]
                    self.rss.append((rs, record.ALT))
                except KeyError:
                    # if no match is found, let user know
                    print "couldn't match", record.POS
                    pass