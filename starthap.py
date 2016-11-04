from haplotyper import Haplotyper

if __name__ == "__main__":
    app = Haplotyper('data/hpc/test.vcf')
    #app.HapMaker("hgvs")
    app.HapMaker("rs")
    print app.haplotypes
else:
    pass