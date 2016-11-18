#!/usr/bin/python
# -*- coding: utf-8 -*-
from variant import Variant

rsids = []

with open('config/rsidlist', 'r') as f:
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


with open('config/new_design.vcf', 'w') as a, open('config/new_design_indels.vcf','w') as b:

    recs_a = []
    recs_b = []
    recs_x = []
    recs_y = []

    a.write(header)
    b.write(header)

    for rsid in rsids:
        try:
            v = Variant(rsid, 'pharmgkb')
        except:
            print "error loading", rsid
            continue
        v.GetLocation()

        print rsid
        try:
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
        except:
            continue

        rec_tab = '\t'.join(rec) + "\n"

        if '-' in v.ref or "(" in v.alt or'-' in v.alt:
                recs_b.append(rec_tab)
                continue
        if v.chrom == "X":
            recs_x.append(rec_tab)
        elif v.chrom == "Y":
            recs_y.append(rec_tab)
        else:
            recs_a.append(rec_tab)

    lsts = [recs_a, recs_x, recs_y]

    for lst in lsts:
        try:
            lst.sort(key=lambda x: (int(x.split('\t')[0]), int(x.split('\t')[1])))
        except ValueError:
            lst.sort(key=lambda x: (int(x.split('\t')[1])))

        a.writelines(lst)

    # make pretty later
        b.writelines(recs_b)