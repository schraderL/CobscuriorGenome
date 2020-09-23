
# calculate windows
```bash
rawbase=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/
base=$rawbase/analyses/visualizeTEislands
TEgff=$rawbase/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results/CobsTEsorted.csv.gff
genome=$rawbase/results/Cobs.alpha.2.1.fa
teisl=$rawbase/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v.2.1.TEislands/te_islands.bed
cd $base
mkdir windows
ln -s $genome .
samtools faidx Cobs.alpha.2.1.fa
cat Cobs.alpha.2.1.fa.fai |awk  -F $'\t' 'BEGIN {OFS = FS} {print $1,1,$2}'|sort -k3,3 -nr > CobsA2.1.bed
bedtools makewindows -b CobsA2.1.bed -w 100000 |bedtools sort > windows/100kwindows.bed
bedtools makewindows -b CobsA2.1.bed -w 10000 |bedtools sort > windows/10kwindows.bed
bedtools makewindows -b CobsA2.1.bed -w 20000 |bedtools sort > windows/20kwindows.bed


# get coverage of genes
gff=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.geneCombined/Cobs.alpha.v.2.1.geneannotation.1.5/Cobs.alpha.v.2.1.geneannotation.1.5.gff3
cat $gff| awk '{if ($3=="exon") print $0}'|bedtools sort -i - > tmpEXON.gff
sortBed -i tmpEXON.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > Genes.windows.coverage.bed

# get coverage of certain TEs
cat $TEgff|perl -pe 's/.*Name=(.*?);.*/$1/g' |cut -f 1 -d "/"|sort|uniq -c|perl -pe 's/^ +//g'|tr " " "\t" |sort -k1,1 -nr

#20756	DNA
#9310	LTR
#4946	LINE
#2696	Retro
#1994	Unknown
#1636	SINE

#LTRs
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=LTR/ {print $0}'|bedtools sort -i - > CobsA2.1.LTR.gff
sortBed -i CobsA2.1.LTR.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > LTR.windows.coverage.bed

#DNA
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=DNA/ {print $0}'|bedtools sort -i - > CobsA2.1.DNA.gff
sortBed -i CobsA2.1.DNA.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > DNA.windows.coverage.bed

#LINEs
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=LINE/ {print $0}'|bedtools sort -i - > CobsA2.1.LINE.gff
sortBed -i CobsA2.1.LINE.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > LINE.windows.coverage.bed

#Retros
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=Retro/ {print $0}'|bedtools sort -i - > CobsA2.1.Retro.gff
sortBed -i CobsA2.1.Retro.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > Retro.windows.coverage.bed

#Unknowns
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=Unknown/ {print $0}'|bedtools sort -i - > CobsA2.1.Unknown.gff
sortBed -i CobsA2.1.Unknown.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > Unknown.windows.coverage.bed

#SINEs
cat $TEgff| awk   -F $'\t' 'BEGIN {OFS = FS} $9 ~ /^Name=SINE/ {print $0}'|bedtools sort -i - > CobsA2.1.SINE.gff
sortBed -i CobsA2.1.SINE.gff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > SINE.windows.coverage.bed

#all TEs
sortBed -i $TEgff | gff2bed  |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > TEs.windows.coverage.bed

#TE islands
sortBed -i $teisl   |bedmap --echo --bases-uniq-f windows/10kwindows.bed - |tr "|" "\t" > TEislands.windows.coverage.bed
```
