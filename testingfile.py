import numpy as np

import os

from Bio.Align.Applications import ClustalwCommandline as c

from Bio import AlignIO

cline = c("clustalw2", infile = "data/alignments/PA18_aln.fasta")

stdout, stderr = cline()

