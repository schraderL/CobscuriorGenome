# Genome paper metrics
## quast
*Version: 4.6.3*
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/genomePaper/
cd $base/quast
assembly=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.1.fa
cd /global/homes/jg/schradel/data/Pcal/minION/QUAST
nice quast.py $assembly -o CobsAlpha.v2.1.quast.out --eukaryote
```

```
Assembly                    Cobs.alpha.2.1
# contigs (>= 0 bp)         127
# contigs (>= 1000 bp)      127
# contigs (>= 5000 bp)      125
# contigs (>= 10000 bp)     121
# contigs (>= 25000 bp)     105
# contigs (>= 50000 bp)     91
Total length (>= 0 bp)      193051228
Total length (>= 1000 bp)   193051228
Total length (>= 5000 bp)   193047645
Total length (>= 10000 bp)  193025568
Total length (>= 25000 bp)  192755434
Total length (>= 50000 bp)  192237042
# contigs                   127
Largest contig              13148674
Total length                193051228
GC (%)                      41.02
N50                         6290588
N75                         4487289
L50                         11
L75                         21
# N's per 100 kbp           94.76
```


## TE islands
```bash

islands=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v.2.1.TEislands/te_islands.bed
# number of islands
wc -l $islands
#35
# number of scaffolds containing islands
cut -f 1 $islands|sort|uniq |wc -l
# 27
```

```bash
# size of islands
cat $islands |awk 'NF > 0 { print $1 "\t" ($3 - $2) }'|sort -nr -k2,2
```

```bash
scaffold1	2050000
scaffold2	1500000
scaffold15	1420358
scaffold17	1200000
scaffold20	1100000
scaffold14	1077151
scaffold5	1062622
scaffold4	1000000
scaffold19	1000000
scaffold16	1000000
scaffold25	906168
scaffold6	900000
scaffold9	890588
scaffold8	805125
scaffold7	800000
scaffold21	684715
scaffold24	631917
scaffold4	600000
scaffold27	600000
scaffold12	555492
scaffold30	500000
scaffold3	500000
scaffold25	500000
scaffold23	500000
scaffold22	500000
scaffold26	429807
scaffold2	360777
scaffold10	332196
scaffold21	300000
scaffold13	284943
scaffold18	236763
scaffold17	216616
scaffold9	90588
scaffold20	87289
scaffold22	41709
```

```bash
# sum of islands
cat $islands |awk 'NF > 0 {sum += $3-$2} END {print sum/1000000}'
#24.6648
```

```bash
# number of protein coding genes in TE islands
islands=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v.2.1.TEislands/te_islands.bed
gff=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.geneCombined/Cobs.alpha.v.2.1.geneannotation.1.5/Cobs.alpha.v.2.1.geneannotation.1.5.gff3

bedtools intersect -wao -a $islands -b $gff |awk '{if ($6=="gene") print $0}'|wc -l
#2118
bedtools intersect -v -b $islands -a $gff|awk '{if ($3=="gene") print $0}'|wc -l
#18855
```

```R

TEislandSize<-24.66
TEislandGene<-2118
geneCount<-20966
genomeSize<-193.05
LDRgenesObs<-geneCount-TEislandGene
TEgenesExp<-geneCount*(TEislandSize/genomeSize)
LDRgenesExp<-geneCount*(1-TEislandSize/genomeSize)
ft<-fisher.test(matrix(c(TEislandGene,LDRgenesObs , TEgenesExp, LDRgenesExp), 2, 2))

#Fisher's Exact Test for Count Data
#
#data:
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#0.7219377 0.8156277
#sample estimates:
#odds ratio
#0.7673988

```

## Retrieve Functional annotations of TE island genes
```bash
base=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/analyses/GOenrichment/TEislands
bedtools intersect -b $islands -a $gff |awk '{if ($3=="gene") print $0}'|perl -pe 's/.*ID=(.*?);.*/$1/g' > TEisland.genes.lst

grep -f TEisland.genes.lst /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.geneCombined/Cobs.alpha.v.2.1.geneannotation.1.5/Cobs.alpha.v.2.1.geneannotation.1.5.functionalAnnotation/Cobs.alpha.v.2.1.geneannotation.1.5.interproFormat.tsv
#in R
./analyses/GOenrichment/TEisland.GOenrichment.Rmd
```
