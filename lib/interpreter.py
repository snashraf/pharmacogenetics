<<<<<<< HEAD
#!/usr/bin/python
# -*- coding: utf-8 -*-

# --------------------------------------------------------------------------


class Interpreter:

    def __init__(self, patientObj):

        self.p = patientObj

        self.advice = {}

        self.conn = sqlite3.connect('pharmacogenetics.db')

        self.sql = self.conn.cursor()

        self.conn.commit()

    def HapScorer(self):

        it = ['p.al1', 'p.al2']

        self.sql.execute('select distinct gid, symbol from genes')

        results = self.sql.fetchall()

        for (gid, symbol) in results:

            print symbol

            print '---------------------'

            for i in it:

                print i

                # get exact matches:

                self.sql.execute("select distinct a.starname, p.hapid, %s from patienthaps p \
											join alleles a on a.hapid=p.hapid \
											where gid=? order by %s asc"
                                 % (i, i), (gid, ))

                results = self.sql.fetchall()

                if len(results) == 0:

                    continue
                else:

                    for (starname, hapid, alscore) in results:

                        lastval = []

                        self.sql.execute(
                            'select distinct starhaps, guid from drugpairs where gid = ?', (gid, ))

                        vals = self.sql.fetchall()

                        for (starhaps, guids) in vals:

                            if starname in starhaps:

                                lastval.append(guids)

                        print hapid, '/', starname, '|', alscore, '|', \
                            list(set(lastval))

                print '---------------------'

    def DoseGuideCheck(self):

        authobj = Authenticate()

        doseInfo = None

        doseGuide = None

        self.sql.execute('drop table if exists annotations')

        self.sql.execute('create table if not exists annotations(varid text, did text, allele text, advice text, loe text)'
                         )

        self.sql.execute("select distinct v.varid, d.did, p.ref, p.alt, p.call, p.start from variants v join patientvars p on p.start = v.start join drugpairs d on d.gid = v.gid where d.rsids like ('%' || v.rsid || '%') group by d.did"
                         )

        for (
            varid,
            did,
            ref,
            alt,
            call,
            start,
        ) in self.sql.fetchall():

            if call == '0/0':

                allele = ref + ref
            elif call == '0/1':

                allele = ref + alt
            elif call == '1/1':

                allele = alt + alt

            print allele

            results = PGKB_connect(authobj, 'clinAnno', varid, did)

            if results is not None:

                for phen in results:

                    for al in phen['allelePhenotypes']:

                        if allele in al['allele']:

                            advice = al['phenotype']

                            loe = phen['levelOfEvidence']['term']

                            item = (varid, did, allele, advice, loe)

                            print 'Level of evidence:', loe
                            print advice

                            self.sql.execute(
                                'insert into annotations values(?,?,?,?,?)', item)
                        else:

                            pass

        self.conn.commit()


pat = Patient('data/test.g.vcf.gz')

print 'patient created'

app = Interpreter(pat)
app.HapScorer()
=======
#!/usr/bin/python
# -*- coding: utf-8 -*-

# --------------------------------------------------------------------------

class Interpreter:

    def __init__(self, p):

        self.advice = {}

        self.conn = p.conn

        self.sql = p.conn.cursor()

        self.conn.commit()


    def Genotyper(self):
    	
    	# For each drug:
    		
    	self.sql.execute("SELECT DISTINCT DrugID FROM Pairs")
    	
    	for (did,) in self.sql.fetchall():
    		
    		# Select guidelines involved
    		
    		self.sql.execute("SELECT DISTINCT GuID from Guidelines WHERE DrugID = ?", (did,))
    		
    		for (guid,) in self.sql.fetchall():
    	
    			# fetch involved genes
    	
			self.sql.execute("SELECT DISTINCT GeneID, o.GeneSymbol FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
		    		
		    	# Fetch involved haplotypes and scores
			
			genotype = {}
		
			for (gid, genesymbol) in self.sql.fetchall():
		
				self.sql.execute("SELECT DISTINCT p.HapID, Distance1, Distance2,h.hgvs, starname FROM PatHaplotypes p JOIN Haplotypes h on h.HapID = P.HapID WHERE h.GeneID = ?", (gid,))
				
				haplotypes = self.sql.fetchall()
				
				self.sql.execute("SELECT DISTINCT Starname FROM GuideOptions o JOIN Genes g on g.GeneSymbol = o.GeneSymbol WHERE GuID = ?", (guid,))
				
				options = [starname[0] for starname in self.sql.fetchall()]
				
				print options
				
				if len(haplotypes) == 0:
					
					continue
						
				genotype[genesymbol] = {"al1":
						{"starname":"A", "distance":10.0},
						    "al2":
						{"starname":"A", "distance":10.0}}
					
				for (hapid, al1, al2, hgvs, starname) in haplotypes:
							    			
		    		# FOR ALLELE IN PATIENT:
		    		
		    			for i, alleleScore in enumerate([float(al1), float(al2)]):
		    				
		    				curLoc = genotype[genesymbol]["al{}".format(i+1)]
		    				
		    				curScore = curLoc["distance"]
		
						if curScore > alleleScore and starname in options:
								
							curLoc["starname"] = starname
		
							curLoc["distance"] = alleleScore
			    		
		    		formatted_genotype = []
		    		
		    		for gene, allele in genotype.items():
		    			
		    			alleles = []
		
		    			for subdict in allele.values():
		    				
		    				alleles.append(subdict["starname"])
		    			
		    			allele_string = "/".join(alleles)
		    			
		    			formatted_genotype.append(":".join([gene, allele_string]))
		    			
		    		string_genotype = ";".join(formatted_genotype)
		    		
		    		print guid
		    		
		    		print string_genotype
		    		
		    		raw_input("Print enter to continue...")
		    				
		    		# Find matching advice
		    		
		    		# Save to advice table 'PatGuidelines' (DrugID, GeneID, Category(Metabolizer type), Advice)


    def AnnFinder(self):
    	
    	pass
    	
    	# For each drug:
    		
    		#  Find involved genes
    		
    		# FOR GENE IN LIST:
    		
    			# Is there an identified guideline for the patient?
    		
    			# IF YES:
    		
    				# continue with next gene
    		
    			# IF NO:
    		
    				# can i find a haplotype for this gene?
    		
    				# if yes:
    		
    					# which closest haplotype has an annotation available?
    		
    				# if no:
    		
    					# continue on
    		
    		# go through 'variant'-tagged annotations
    		
    		# Find involved annotations and match with patientvars
    		
    		# DO SNPS and OTHERS seperately! GA/G becomes A/- for example, we need this for the annotation
    		
    		# Look up annotation
    		
    		# Get result (jinja?) and add to PatAnnotations
    		

# ---------------------------- NEXT STEP: ReportMaker --------------------------------
>>>>>>> b00c19f362a80f284a99c1d4ace98f73dcd9e9ac
