import numpy as np

import os

from Bio.Align.Applications import ClustalwCommandline as c

from Bio import AlignIO

cline = c("clustalo", infile = "data/alignments/PA18_aln.fasta")

stdout, stderr = cline()

for group in groups:
					
					if group['value'] == 0.0:
						
						for hap in group['haps']:
							
							if "REF" in hap:
								
								print "Reference"
								
								break
							
							elif "Yes" in hap:
									
								print "Exact match:", hap
							
								return
									
							elif "No" in hap:
								
								print "Exact match without guideline:", hap
								
								return
					
					elif group['value'] != 0.0:
						
						# get next highest value with a true in there
						
						for hap in group['haps']:
							
							if "Yes" in hap:
									
									print "Closest match:", hap
									
									return
						
							elif "No" in hap:
								
									print "Closest match without guideline:", hap
							
							else:
								
								print "Test"
