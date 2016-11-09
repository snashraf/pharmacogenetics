import pickle

with open("variants.pickle", "rb") as file:
    variants = pickle.load(file)

with open("genes.pickle", "rb") as file:
    genes = pickle.load(file)

with open("chemicals.pickle", "rb") as file:
    chemicals = pickle.load(file)

nohap = []

for gene in genes:
    haps = gene["alleles"]
    name = gene["symbol"]
    drugs = gene["drugs"]
    print drugs
    if len(haps) == 0:
        nohap.append(name)

print chemicals