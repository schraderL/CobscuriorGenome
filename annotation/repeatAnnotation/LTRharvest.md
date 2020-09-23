# LTRharvest annotation of Cobs.alpha.v2.1

We used LTRharvest to get a dedicated annotation of LTR transpososons in the genome of *C. obscurior* alpha.v2.1.

See http://avrilomics.blogspot.com/2015/09/ltrharvest.html for a description of the process.

LTRharvest is included in the genome tools (gt) package: http://genometools.org/tools.html

## Define environment
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/LTRharvest
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
out=Cobs.alpha.2.1.LTRharvest
gtpath=/global/homes/jg/schradel/software/gt-1.6.1-Linux_x86_64_x86_64-64bit/bin
cd $base
ln -s $genome .

# download gypsyDB collection
wget http://gydb.org/gydbModules/collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip

cd GyDB_collection/profiles/
ln -s /global/homes/jg/schradel/data/REPET/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm .
```

## File: filter_full_domain_set.lua
```bash
name        = "Protein Domain Filter"
author      = "Sascha Kastens"
version     = "1.0"
email       = "mail@skastens.de"
short_descr = "Filters out candidates without protein domains"
description = "Filters out a candidate if it does not contain at " ..
              "least one node of type 'protein_match'."

function filter(gn)
  gfi = gt.feature_node_iterator_new(gn)
  node = gfi:next()
  while not (node == nil) do
    if (node:get_type() == "protein_match") then
      return false
    end
    node = gfi:next()
  end
  return true
end
```


```bash
#http://avrilomics.blogspot.com/2015/09/ltrharvest.html

$gtpath/gt suffixerator -suf -lcp -tis -des -dna -ssp -db Cobs.alpha.2.1.fa -indexname Cobs.alpha.2.1.fa
$gtpath/gt ltrharvest -index Cobs.alpha.2.1.fa -seqids -tabout no > Cobs.alpha.2.1.ltrs.gff3
$gtpath/gt gff3 -sort Cobs.alpha.2.1.ltrs.gff3 > Cobs.alpha.2.1.ltrs_sorted.gff3
$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -aaout -outfileprefix ltrs_sorted -seqfile Cobs.alpha.2.1.fa -matchdescstart < Cobs.alpha.2.1.ltrs_sorted.gff3 > Cobs.alpha.2.1.ltrs.digested.gff3

# filter full domain set
$gtpath/gt select -rule_files filter_full_domain_set.lua -- < Cobs.alpha.2.1.ltrs.digested.gff3 > Cobs.alpha.2.1.ltrdigest.filtered.gff3
$gtpath/gt ltrharvest -index Cobs.alpha.2.1.fa -tabout no > Cobs.alpha.2.1.ltrs.legacy.gff3
# this only works on the legacy output without proper scaffold IDs
$gtpath/gt ltrclustering -psmall 80 -plarge 30 -o lologre_ltrclustering.out -seqfile Cobs.alpha.2.1.fa lologre_ltrdigest.gff3

# combine output from ltrclustering with ltrharvest
cat Cobs.alpha.2.1.ltrs.digested.gff3  |grep "##sequence-region"|sort -k 4,4 -nr > scf.seq.tsv
cat lologre_ltrclustering.out  |grep "##sequence-region"|sort -k 4,4 -nr > legacy.seq.tsv
paste scf.seq.tsv legacy.seq.tsv|perl -pe 's/ +/\t/g'|cut -f 2,6 > conversionTable.tsv
awk 'FNR==NR{a[$2]=$1;next}{print a[$1]"\t"$0}' conversionTable.tsv lologre_ltrclustering.out|cut -f 1,3-11|awk 'NF' > Cobs.alpha.2.1.ltrs.clustered.gff3

cat Cobs.alpha.2.1.ltrs.digested.gff3  |egrep "^#"|egrep "###" -v > Cobs.alpha.2.1.ltrs.digested.header
cat Cobs.alpha.2.1.ltrs.digested.header Cobs.alpha.2.1.ltrs.clustered.gff3 > Cobs.alpha.2.1.ltrs.clustered.header.gff3
$gtpath/gt select -rule_files filter_full_domain_set.lua -- < Cobs.alpha.2.1.ltrs.clustered.header.gff3 > Cobs.alpha.2.1.ltrs.clustered.filtered.gff3

```


## Processing LTR annotation
This was run locally, after downloading the LTRharvest results from the cluster.
```bash

base=/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.LTRharvest/
cd $base
# calculate frequencies of different LTR families
grep ltrfam Cobs.alpha.2.1.ltrs.clustered.filtered.gff3|perl -pe 's/.*ltrfam=(ltrfam_[0-9]+).*/$1/g'|sort |uniq -c|perl -pe 's/^ +([0-9]+) /$1\t/g'|sort -nr -k 1,1  > ltrfam.freq.tsv
```
`head ltrfam.freq.tsv`
>24	ltrfam_8
22	ltrfam_15
19	ltrfam_4
16	ltrfam_57
11	ltrfam_5
10	ltrfam_31
10	ltrfam_25
9	ltrfam_40
9	ltrfam_16
7	ltrfam_73

### Analyse individual families
```bash

cd $base/blast
# ltrfam_8 is the most abundant family
grep ltrfam_8$ ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3
grep "ltrfam=ltrfam_8$" ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3 |bedtools intersect -b - -a ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3|awk '{if ($3=="long_terminal_repeat") print $0}' > ltrfam8.LTRs.gff3

#grep "ID=repeat_region12;" ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3 |bedtools intersect -b - -a ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3|awk '{if ($3=="long_terminal_repeat") print $0}' > repeat_region12.LTRs.gff3

#grep ltrfam_15 ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3
#grep "ID=repeat_region35;" ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3 |bedtools intersect -b - -a ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3|awk '{if ($3=="long_terminal_repeat") print $0}' > repeat_region35.LTRs.gff3


# retrieve long_terminal_repeat sequences
bedtools getfasta -bed ltrfam8.LTRs.gff3 -fi ../../../results/Cobs.alpha.2.1.fa > ltrfam8.LTRs.fa

#bedtools getfasta -bed repeat_region12.LTRs.gff3 -fi ../../../results/Cobs.alpha.2.1.fa > repeat_region12.LTRs.fa
#bedtools getfasta -bed repeat_region35.LTRs.gff3 -fi ../../../results/Cobs.alpha.2.1.fa > repeat_region35.LTRs.fa

mafft --ep 0 --genafpair --maxiterate 1000  --clustalout ltrfam8.LTRs.fa > ltrfam8.LTRs.aln

#mafft --ep 0 --genafpair --maxiterate 1000  --clustalout repeat_region12.LTRs.fa
#mafft --ep 0 --genafpair --maxiterate 1000  --clustalout repeat_region35.LTRs.fa > repeat_region35.LTRs.aln
cat ltrfam8.LTRs.aln |perl -pe 's/(scaffold[0-9]+)\:([0-9]+)/$1.$2/g' > ltrfam8.LTRs.aln
trimal -in ltrfam8.LTRs.aln -st 1 -gt 0.6
"ccgtcatctgttgaggt"
"tgcccgcacctcgacg"

# best so far ()

##  bases 1-20 of LTR_retrotransposon12 long_terminal_repeat (scaffold9:6007484-6008227) see
>scaffold9:6007483-6007503
TGTGGCGGATCCCTCCGACG
##  bases 21-40 of LTR_retrotransposon12 long_terminal_repeat (scaffold9:6007484-6008227)
>scaffold9:6007504-6007524
CTACACGTCGGATGACGCCA

GCGAAGGAGGAAGCAAACTC

GTGTGGGTGCGGATTTGCCC

samtools faidx ltrfam8.LTRs.fa

# Divide LTR_retrotransposon12 long_terminal_repeat (scaffold9:6007484-6008227) into chunks of 20 bps

length=$((6008227-6007483-20))
for i in $(seq 1 $length);
do
  start=$((1+$i))
  end=$((20+$i))
  echo scaffold9:6007483-6008227:$start-$end
done > scaffold9.6007483.6008227.regions
# extract 20 bp chunks of this LTR
samtools faidx ltrfam8.LTRs.fa -r scaffold9.6007483.6008227.regions > scaffold9.6007483.6008227.regions.fa
# blast chunks of 20bp against the genome
blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -num_threads 4 -query scaffold9.6007483.6008227.regions.fa -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100 > scaffold9.6007483.6008227.regions.bls

# check how many perfect hits each chunk produces
cut -f 1 scaffold9.6007483.6008227.regions.bls |sort|uniq -c|perl -pe 's/^ +([0-9]+) /$1\t/g'|sort -nr -k 1,1 > scaffold9.6007483.6008227.regions.list

echo "TGTGGCGGATCCCTCCGACG"|blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query - -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100|wc -l
echo "CTACACGTCGGATGACGCCA"|blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query - -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100|wc -l
echo "GCGAAGGAGGAAGCAAACTC"|blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query - -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100|wc -l
echo "GTGTGGGTGCGGATTTGCCC"|blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query - -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100|wc -l


echo "tcactaatttgcggga"|blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query - -word_size 7 -dust no -evalue 5 -outfmt 6  -perc_identity 100 -qcov_hsp_perc 100|wc -l

# retrieved LTR of most abundant LTRfam from IGV and extracted first 20+20 bases for a primer

blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query testPrimer1.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6 > testPrimer1.bls
blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query testPrimer2.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6 > testPrimer2.bls
cat testPrimer1.bls |awk '{print $2"\t"$9"\t"$10"\t"$1";"$3";"$4";"$5";"$6";"$7";"$8";"$11";"$12}' > testPrimer1.bed


primersearch \
-seqall /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa \
-infile /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.LTRharvest/blast/primerPair1.fa \
-outfile test.primersearch \
-mismatchpercent 0
less test.primersearch
```

## Blast-based mcl clustering of ltrs
```bash
cd ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.LTRharvest/blast
ln -s ../fastas/ltrs_sorted_3ltr.fas .
makeblastdb -in ltrs_sorted_3ltr.fas -dbtype nucl
blastn -db  ltrs_sorted_3ltr.fas -query  ltrs_sorted_3ltr.fas -evalue 1e-10 -dust no -outfmt 6 > ltrs_sorted_3ltr.bls
cut -f 1,2,11 ltrs_sorted_3ltr.bls > seq.abc
mcxload -abc seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.mci -write-tab seq.tab
mcl seq.mci -I 1.4
mcxdump -icl out.seq.mci.I14 -tabr seq.tab -o dump.seq.mci.I14

 head -n 1 dump.seq.mci.I14 |tr "\t" "\n" > IDs.txt
 bioawk -cfastx 'BEGIN{while((getline k <"IDs.txt")>0)i[k]=1}{if(i[$name])print ">"$name"\n"$seq}' ltrs_sorted_3ltr.fas|mafft -


#>clusteredAln
#agatattgacaaa


blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query clusteredLTRs.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6
```

## Align entire LTRs of fam 8
```bash
grep "ltrfam_8$" ../Cobs.alpha.2.1.ltrs.clustered.filtered.gff3 > ltrfam_8.gff3
bedtools getfasta -bed ltrfam_8.gff3 -fi ../../../results/Cobs.alpha.2.1.fa > ltrfam_8.fa

mafft --ep 0 --genafpair --maxiterate 1000  --clustalout ltrfam_8.fa > ltrfam_8.aln
cat ltrfam_8.aln |perl -pe 's/(scaffold[0-9]+)\:([0-9]+)/$1.$2/g' > ltrfam_8.reformatted.aln
trimal -in ltrfam_8.reformatted.aln -st 1 -gt 0
# selct

blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query clusteredLTRs.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6
# does not work
```
https://academic.oup.com/nar/article/38/suppl_2/W313/1111009
https://link.springer.com/article/10.1007%2Fs11295-008-0182-9

<!----

```R
library(DECIPHER)
genome<-readDNAStringSet("/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa")
primers<-DNAStringSet(c("TGTGGCGGATCCCTCCGACG","CTCTCTCTCTCTCTCTCTCTA"))
test<-AmplifyDNA(primers,genome,annealingTemp=55, P=4e-7, maxProductSize=5000)
```


blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query repeat_region12.LTRs.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6

# best primer so far
TGTGGCGGATCCCTCCGACG
CTACACGTCGGATGACGCCA
GCGAAGGAGGAAGCAAACTC
see ```/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.LTRharvest/blast/bestPrimers.fa```
```
cd ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.LTRharvest/blast
blastn -db /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa -query bestPrimers.fa -word_size 7 -dust no -evalue 0.01 -outfmt 6

```




## Run LTRharvest
```bash
cd $base
$gtpath/gt suffixerator -db Cobs.alpha.2.1.fa -indexname Cobs.alpha.2.1.fa -tis -suf -lcp -des -ssp -sds -dna
$gtpath/gt ltrharvest -index Cobs.alpha.2.1.fa  -gff3 $out.gff -out $out.fa -seqids yes -tabout no > ltrharvest.out
$gtpath/gt ltrharvest -index Cobs.alpha.2.1.fa  -gff3 $out.2.gff -out $out.fa -seqids no -tabout no > ltrharvest.2.out


$gtpath/gt gff3 -sort ltrharvest.2.out  > ltrharvest.2.sorted.gff
$gtpath/gt gff3 -sort ltrharvest.out  > ltrharvest.sorted.gff
$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -outfileprefix Cobs.alpha.v2.1.ltrdigest ltrharvest.2.sorted.gff  Cobs.alpha.2.1.fa > Cobs.alpha.2.1.ltrdigest.gff
$gtpath/gt select -rule_files filter_full_domain_set.lua -- < Cobs.alpha.2.1.ltrdigest.gff > Cobs.alpha.2.1.ltrdigest.filtered.gff

$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -outfileprefix Cobs.alpha.v2.1.ltrdigest ltrharvest.sorted.gff  Cobs.alpha.2.1.fa > Cobs.alpha.2.1.ltrdigest.seqIDs.gff


#$gtpath/gt gff3 -sort $out.gff  > $out.sorted.gff
#$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -outfileprefix Cobs.alpha.v2.1.ltrdigest -encseq Cobs.alpha.2.1.fa -matchdescstart < $out.sorted.gff > Cobs.alpha.2.1.LTRdigest.gff
#$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -outfileprefix Cobs.alpha.v2.1.ltrdigest $out.sorted.gff Cobs.alpha.2.1.fa > Cobs.alpha.2.1.ltrdigest.gff

#$gtpath/gt ltrdigest -hmms GyDB_collection/profiles/*hmm -outfileprefix Cobs.alpha.v2.1.ltrdigest -encseq Cobs.alpha.2.1.fa -matchdescstart < ltrharvest.2.out > Cobs.alpha.2.1.LTRdigest.gff

mv Cobs.alpha.v2.1.ltrdigest_*.fas ltrdigest/
```
