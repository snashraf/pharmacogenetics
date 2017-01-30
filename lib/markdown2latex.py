import re

ex1 ='''
The Royal Dutch Pharmacists Association - Pharmacogenetics Working Group has evaluated therapeutic dose recommendations for azathioprine based on TPMT genotype [PMID:21412232].  They recommend selecting an alternative drug or reducing the initial dose for patients carrying inactive alleles.

| Phenotype (Genotype) | Therapeutic Dose Recommendation | Level of Evidence | Clinical Relevance |
| --- | --- | --- | --- |
| IM (one inactive allele: *2, *3, *4-*18) | Select alternative drug or reduce dose by 50%. Increase dose in response of hematologic monitoring and efficacy. | Published controlled studies of good quality* relating to phenotyped and/or genotyped patients or healthy volunteers, and having relevant pharmacokinetic or clinical endpoints. | Clinical effect (S): Failure of lifesaving therapy e.g. anticipated myelosuppression; prevention of breast cancer relapse; arrhythmia; neutropenia &lt; 0.5x10^9^/l; leucopenia &lt; 1.0x10^9^/l; thrombocytopenia &lt; 25x10^9^/l; life-threatening complications from diarrhea. |
| PM (two inactive alleles: *2, *3, *4-*18) | Select alternative drug or reduce dose by 90%. Increase dose in response of hematologic monitoring and efficacy. | Published controlled studies of good quality* relating to phenotyped and/or genotyped patients or healthy volunteers, and having relevant pharmacokinetic or clinical endpoints. | Clinical effect (S): death; arrhythmia; unanticipated myelosuppression. |

- *See [Methods ](http://www.pharmgkb.org/home/dutch_pharmacogenetics_working_group.jsp) or [PMID: 18253145] for definition of &#34;good quality.&#34;
- S: statistically significant difference.
'''

ex2 = '''
##  May 2014 Update on PharmGKB
- The CPIC authors recommend that the _DPYD*4_, _*5_, _*6_ and _*9A_ alleles be categorized as &#34;normal&#34; activity, in part based upon the recent publication [Comparative Functional Analysis of DPYD Variants of Potential Clinical Relevance to Dihydropyrimidine Dehydrogenase Activity](http://www.ncbi.nlm.nih.gov/pubmed/?term=24648345).

## December 2013 Publication

_Accepted article preview online August 2013; Advance online publication October 2013._
- Guidelines regarding the use of pharmacogenomic tests in dosing for fluoropyrimidines have been published in _Clinical Pharmacology and Therapeutics_ by the Clinical Pharmacogenetics Implementation Consortium ([CPIC](/contributors/consortia/cpic_profile.jsp)).
- These guidelines are applicable to:
  - at the time of this writing, there are no data available on the possible role of DPYD*2A, *13, or [variant:rs67376798] in 5-fluorouracil toxicities in pediatric patient populations; however, there is no reason to suspect that DPYD variant alleles would affect 5-fluorouracil metabolism differently in children compared to adults.
- Excerpt from the fluoropyrimidine dosing guideline based on DPYD genotype:
  - &#34;The strength of the dosing recommendations is based on the fact that some variants (DPYD*2A, *13, and [variant:rs67376798]) clearly affect DPD activity, and DPD activity is clearly related to 5-fluorouracil clearance, and 5-fluorouracil exposure is associated with its toxic effects. Therefore, reduction of fluoropyrimidine dosage in patients with these variants may prevent severe and possibly life-threatening toxicities. However, available evidence does not clearly indicate a degree of dose reduction needed to prevent fluoropyrimidine related toxicities...\[Based on literature review (see full manuscript),] our recommendation is to start with at least a 50% reduction of the starting dose followed by an increase in dose in patients experiencing no or clinically tolerable toxicity to maintain efficacy, a decrease in dose in patients who do not tolerate the starting dose to minimize toxicities or pharmacokinetic guided dose adjustments (if available). Patients who are homozygous for DPYD*2A, *13, or [variant:rs67376798] may demonstrate complete DPD deficiency and the use of 5-fluorouracil or capecitabine is not recommended in these patients.&#34;
- Download and read:
  - [Clinical Pharmacogenetics Implementation Consortium Guidelines for Dihydropyrimidine Dehydrogenase Genotype and Fluoropyrimidine Dosing](https://github.com/PharmGKB/cpic-guidelines/raw/master/fluoropyrimidines/2013/23988873.pdf)
  - [2013 supplement](https://github.com/PharmGKB/cpic-guidelines/raw/master/fluoropyrimidines/2013/23988873-supplement.pdf)

##  Table 1: Recommended dosing of fluoropyrimidines by genotype/phenotype.

_Adapted from Tables 1 and 2 of the 2013 guideline manuscript._

| Phenotype (genotype) | Examples of diplotypes | Implications for phenotypic measures | Dosing recommendations | Classification of recommendations ^a^ |
| --- | --- | --- | --- | --- |
| Homozygous wild-type or normal, high DPD activity (two or more functional *1 alleles) | *1/*1 | Normal DPD activity and &#34;normal&#34; risk for fluoropyrimidine toxicity | Use label-recommended dosage and administration | Moderate |
| Heterozygous or intermediate activity (~3-5% of patients), may have partial DPD deficiency, at risk for toxicity with drug exposure (one functional allele *1, plus one nonfunctional allele - *2A, *13 or rs67376798A ^c^) | *1/*2A; *1/*13; *1/ rs67376798A ^c^) | Decreased DPD activity (leukocyte DPD activity at 30% to 70% that of the normal population) and increased risk for severe or even fatal drug toxicity when treated with fluoropyrimidine drugs | Start with at least a 50% reduction in starting dose followed by titration of dose based on toxicity ^b^ or pharmacokinetic test (if available) | Moderate |
| Homozygous variant, DPD deficiency (~0.2% of patients), at risk for toxicity with drug exposure (2 nonfunctional alleles - *2A, *13 or rs67376798A ^c^) | *2A/*2A; *13/*13; rs67376798A ^c^ / rs67376798A ^c^ | Complete DPD deficiency and increased risk for severe or even fatal drug toxicity when treated with fluoropyrimidine drugs | Select alternate drug | Strong |

^a^ Rating scheme described in 2013 supplement.

^b^ Increase the dose in patients experiencing no or clinically tolerable toxicity to maintain efficacy; decrease the dose in patients who do not tolerate the starting dose to minimize toxicities.

^c^ Note that the rs67376798A allele refers to the allele on the positive chromosomal strand.  This is important because DPYD is on the minus chromosomal strand and [variant:rs67376798] is a T/A snp.  Therefore, the T allele on the gene confers the deficiency, while the complement on the positive chromosomal strand (A allele) is indicative of deficiency.
'''

replacedict = {
"\\":"",
"##":"",
"_":"",
"&#24":"",
"&#34":"",
"[":"",
"]":"",
"%":"\%",
"^/l":"/L",
"^":"",
";":", ",
"&lt;":":"
}

broken = ex2.split("\n")
intro = []
table = []
bullet = []
latex = ""

for line in broken:
    for a, b in replacedict.items():
        line = line.replace(a, b)
    line = line.lstrip("-, ")
    starter =  line[0:1]
    if "|" in starter:
        table.append(line)
        cols = line.count("|") - 1
    elif "-" in starter or "^" in starter:
        bullet.append(line)
    else:
        if not re.match(r'^\s*$', line):
            intro.append(line.strip("\n"))

for line in intro:
    latex += "\item {} \\newline\n".format(line)

latex += "\\vspace{1pt}\\newline"

latex += '''
\scriptsize
\\begin{center}
\\begin{tabularx}{0.9\\textwidth}{ b%s }
''' %("s" * cols)

for i, line in enumerate(table):
    items = [item for item in line.split("|") if item != "" and item != "---"]
    if i == 0:
        # headers
        items = ["\\textbf{%s}" %item for item in items]
    if i == 1:
        continue
    latex += "&".join(items) + "\\\\"
    if i < len(table) - 1:
        latex += '''
\\vspace{1pt}\\\\
\hline \\\\
\\vspace{1pt}\\\\
        '''
    elif i == len(table) - 1:
        latex += '''
\end{tabularx}
\end{center}
\\normalsize
\\vspace{10pt}
        '''
return latex
