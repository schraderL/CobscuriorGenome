# Calculate TE coverage across windows in Cobs.alpha.v2.1

```bash
base=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/TEanalysis/

# prepare karyotype file
export genomeFa=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa
export abbr=CobsA2.1

cat $genomeFa.fai |awk -F '\t' '{print "chr","-","Cobs."$1,$1,0,$2,"black",10000}'|tr " " "\t" > $base/$abbr.karyotype.txt
cat $genomeFa.fai |cut -f 1,2|sort -k2,2 -nr |awk -v OFS='\t' '{print $1,"0",$2}'> $abbr.bed
bedtools makewindows -b $abbr.bed -w 1000 |bedtools sort > $abbr.1kwindows.bed
bedtools makewindows -b $abbr.bed -w 200000 -s 1000 |bedtools sort > $abbr.200kslidingwindows.bed

# get TE annotation

# TE
rm -rf tmpTE.gff
export PATTERN="CobsAlpha1_REPET_SSRs"
cat CobsAlpha1.REPET.gff3 | awk '{if ($2!=ENVIRON["PATTERN"]) print $0}'|bedtools sort -i - > tmpTE.gff
rm -rf circos/data/plots/TE.windows.bed
bedtools map -a 10kwindows.bed -b tmpTE.gff -c 9 -o sum | sed -e s/'^'/"$abbr."/g > circos/data/plots/TE.windows.bed
bedtools intersect -a 10kwindows.bed -b tmpTE.gff -wao|awk -F '\t' '{a[$1"\t"$2"\t"$3]+=$13}END{for(i in a) print i,a[i],"."}'|tr " " "\t" |bedtools sort -i - > circos/data/plots/TE.windows.coverage.bed

# Calculate unique bases covered by TEs.
rm -f circos/data/plots/TE.windows.coverage.bed
sortBed -i tmpTE.gff | gff2bed  |bedmap --echo --bases-uniq-f 1kwindows.bed - |tr "|" "\t"|sed -e s/'^'/"$abbr."/g> circos/data/plots/TE.windows.coverage.bed

```bash
