import json
import urllib2

gene = "PA145"
uri = "https://api.pharmgkb.org/v1/data/haplotype?gene.accessionId=%s&view=max" %gene
data = urllib2.urlopen(uri)
storage = json.load(data)
hapdict = {gene:[]}

for item in storage:
    rss = []
    star = item['name']
    for allele in item['alleles']:
        try:
            alt = allele['allele']
            rs = allele['location']['displayName']
            rss.append(rs)
        except:
            pass
    hgvs = item['hgvs']
    if "=" not in hgvs:
        hapdict[gene].append({"star":star, "alt":alt, "rs":rss, "hgvs":hgvs})

found_rs = [u'rs3918290', u'rs72549303', u'rs1801158', u'rs1801159', u'rs1801160']

for haplotype in hapdict[gene]:
    hap_rs = haplotype["rs"]
    comparison = set(found_rs).intersection(hap_rs)
    print len(found_rs), "/", len(comparison)

# highest match extractor etc
# highest likelyhood of matching? print star phen matching.