#!/usr/bin/python
# -*- coding: utf-8 -*-
from variant import Variant

rsids = []

with open('config/more_rsids.txt', 'r') as f:
    for line in f:
        rsids.append(line.strip())

with open('config/vcf_header.txt', 'r') as f:
    header = f.readlines()

header = "".join(header)

qual = '255'
filt = 'PASS'
info = 'DP=30'
form = 'GT:GQ'
test = '0/1:99'

with open('config/new_design.vcf', 'w') as f:

    recs = []

    f.write(header)

    for rsid in rsids:
        try:
            v = Variant(rsid, 'pharmgkb')
        except:
            print "error loading", rsid
            continue
        v.GetLocation()

        try:
            if '-' in v.ref or '-' in v.alt:
                continue

        except:
            continue

        print rsid

        rec = [
            str(v.chrom),
            str(v.begin),
            rsid,
            v.ref,
            v.alt,
            qual,
            filt,
            info,
            form,
            test
            ]
        rec_tab = '\t'.join(rec) + "\n"
        recs.append(rec_tab)

    recs.sort(key=lambda x: int(x.split('\t')[0]))
    f.writelines(recs)
