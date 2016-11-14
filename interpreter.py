#!/usr/bin/python
# -*- coding: utf-8 -*-
from patient import Patient, Find
from data import DataCollector
import json
import requests
from requests_oauthlib import OAuth2Session
import urllib
import urllib2
import ast


# -------------------------------------------------------------------------




# -------------------------------------------------------------------------

class Interpreter:

    def __init__(self, f):
        self.p = Patient(f)
        print 'Loading patient:', f
        self.adviceperdrug = []
        self.p.Load()
        self.DrugAdvice('rsid')



tom = Interpreter('data/hpc/test.vcf')
print tom.adviceperdrug