from patient import Patient
import json
import requests
from requests_oauthlib import OAuth2Session
import urllib
import urllib2
import ast

# -------------------------------------------------------------------------

def Authenticate():
    req = ""
    with open("config/auth.txt", "r") as f:
        req = json.load(f)
    url = "https://api.pharmgkb.org/v1/auth/oauthSignIn"
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

class Interpreter(Patient):

    def __init__(self, f):
        Patient.__init__(self, f)
        self.alleles = {}

    def GetAlleles(self):
        for var in self.patientvars:
            allele = ""
            rsid = var["rsid"]
            call = var["call"]["num"]
            ref = var["call"]["ref"]
            alt =  var["call"]["alt"]
            if call == "0/1" or call == "1/0":
                allele = ref + str(alt[0])
            elif call == "1/1":
                allele = str(alt[0]) * 2
            self.alleles[rsid] = allele

    def GetHaplotypes(self):
        pass

    def DrugAdvice(self, mode):
        c = Authenticate()
        if mode == "rsid":
            for k, v in self.rsdrugs:
                for rsid in k:
                    for did in v:
                        uri = "https://api.pharmgkb.org/v1/report/pair/%s/%s/clinicalAnnotation"  %(rsid, did)
                        result = getJson(uri, c)
                        if result is not None:
                            for doc in result:
                                for phen in doc['allelePhenotypes']:
                                    if self.alleles[rsid] in phen['allele']:
                                        print "---", rsid, self.alleles[rsid], "---"
                                        print phen['phenotype']
                                        break
                        else:
                            continue
        elif mode == "haplotype":
            pass

tom = Interpreter('data/hpc/test.vcf')
tom.GetAlleles()
tom.DrugAdvice("rsid")
