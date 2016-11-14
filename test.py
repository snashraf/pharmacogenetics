patienthaps = {"gene":haplotype}

test = ("gen", "hap", "id", "rs1&rs2&rs3")
patientvars = ["rs2", "rs3", "rs4"]
rsids = test[3].split("&")

id = test[2]
comparison = set(patientvars) & set(rsids)

print "match", len(comparison), "/", len(rsids)

if len(comparison) == len(rsids):
    print "whoehoe, match!"
    patienthaps.append(id)

