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

def Authenticate():
    print 'Authenticating...'
    req = ''
    with open('config/auth.txt', 'r') as f:
        req = json.load(f)
    url = 'https://api.pharmgkb.org/v1/auth/oauthSignIn'
    data = urllib.urlencode(req)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    str_response = response.read()
    token = ast.literal_eval(str_response)
    client = OAuth2Session(token=token)
    return client


def getJson(uri, client):
    r = client.get(uri)
    data = r.json()
    if type(data) == dict:
        return None
    elif type(data) == list:
        return data


# -------------------------------------------------------------------------

class Interpreter:

    def __init__(self, f):
        self.p = Patient(f)
        print 'Loading patient:', f
        self.adviceperdrug = []
        self.p.Load()
        self.DrugAdvice('rsid')

    def DrugAdvice(self, mode):
        c = Authenticate()
        if mode == 'rsid':
            for (k, v) in self.p.rsdrugs:
                for rsid in k:
                    for did in v:
                        uri = \
    'https://api.pharmgkb.org/v1/report/pair/%s/%s/clinicalAnnotation' \
                            % (rsid, did)
                        result = getJson(uri, c)
                        if result is not None:
                            for doc in result:
                                for phen in doc['allelePhenotypes']:
                                    if self.p.alleles[rsid] in phen['allele']:
                                        matches = \
    Find(self.p.chemicals, 'id', did)
                                        name = matches[0]['name']
                                        entry = {
                                            'did': did,
                                            'drugname': name,
                                            'rsid': rsid,
                                            'phenotype': phen['phenotype'],
                                            }
                                        self.adviceperdrug.append(entry)
                                        break
                        else:
                            continue
        elif mode == 'haplotype':
            pass

    def PerDrug(self):
        pass


tom = Interpreter('data/hpc/test.vcf')
print tom.adviceperdrug