from haplotyper import Haplotyper

if __name__ == "__main__":
    app = Haplotyper('data/hpc/test.vcf')
    #app.HapMaker("hgvs")
    app.HapMaker("hgvs")
    print app.haplotypes
else:
    pass