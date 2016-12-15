#!/usr/bin/python
# -*- coding: utf-8 -*-
from patient import Patient
from data import DataCollector
import json
import sqlite3
from pgkb_functions import *

class Interpreter:


	def __init__(self, patientObj):

		self.p = patientObj

		self.advice = {}

		self.conn = sqlite3.connect('pharmacogenetics.db')

		self.sql = self.conn.cursor()

		self.conn.commit()

	def HapScorer(self):
			
			it = ["p.al1", "p.al2"]
			
			self.sql.execute("select distinct gid, symbol from genes")
			
			results = self.sql.fetchall()

			for (gid, symbol) in results:

				print symbol
						
				print "---------------------"
					
				for i in it:

					print i
					
					self.sql.execute("select distinct p.hapid, a.starname, p.ref, %s from patienthaps p join alleles a \
					on a.hapid=p.hapid where p.gid=? order by %s, p.ref asc limit 10" %(i,i), (gid,))
					
					results = self.sql.fetchall()
										
					if len(results) == 0:
						
						continue
					
					else:
						
						for (hapid, starname, refscore, alscore) in results:
																					
							lastval = "None"
							
							self.sql.execute("select distinct starhaps, guid from drugpairs where gid = ?", (gid,))
							
							vals = self.sql.fetchone()
							
							starhaps = vals[0]
							
							guid = vals[1]
							
							match = "| Partial |"
							
							if alscore == 0:
							
								match = "| Full |"				
							
							if starname in starhaps and "nan" not in guid:
							
								lastval = guid 
							
							print i
							
							print hapid, "/", starname, "|", refscore, "|", alscore, "|", lastval, match
						
					
					print "---------------------"


	def DoseGuideCheck(self):
		
	
		authobj = Authenticate()
	
		doseInfo = None
	
		doseGuide = None
	
		self.sql.execute("drop table if exists annotations")
		
		self.sql.execute("create table if not exists annotations(varid text, did text, allele text, advice text, loe text)")

		self.sql.execute("select distinct v.varid, d.did, p.ref, p.alt, p.call, p.start from variants v join patientvars p on p.start = v.start join drugpairs d on d.gid = v.gid where d.rsids like ('%' || v.rsid || '%') group by d.did")
		
		for (varid, did, ref, alt, call, start) in self.sql.fetchall():
	
			if call == "0/0":
	
				allele = ref+ref
	
			elif call == "0/1":
	
				allele = ref+alt
	
			elif call == "1/1":
	
				allele = alt+alt
	
			print allele
	
			results = PGKB_connect(authobj, "clinAnno",varid, did) 
	
			if results is not None:
	
				for phen in results:
										
					for al in phen['allelePhenotypes']:
	
						if allele in al['allele']:
	
							advice = al['phenotype']
							
							loe = phen['levelOfEvidence']['term']
							
							item = (varid, did, allele, advice, loe)
							
							print "Level of evidence:", loe
							print advice
							
							self.sql.execute("insert into annotations values(?,?,?,?,?)", item)
								
						else:
	
							pass

		self.conn.commit()

pat = Patient('data/test.g.vcf.gz')

print "patient created"

app = Interpreter(pat)
app.HapScorer()
