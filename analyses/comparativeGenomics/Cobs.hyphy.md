# Cobs hyphy analysis


## setup environment
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/hyphy/
mkdir $base
cd $base
mkdir $base/data
mkdir $base/OGs
mkdir $base/fa
```
## retrieve orthofinder results
```bash

cd $base/data/

ln -s /global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/orthofinder/OrthoFinder/Results_Sep21/Species_Tree/SpeciesTree_rooted_node_labels.txt .
ln -s /global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/orthofinder/OrthoFinder/Results_Sep21/Orthogroups/Orthogroups_SingleCopyOrthologues.txt .
ln -s /global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/orthofinder/OrthoFinder/Results_Sep21/Orthogroups/Orthogroups.tsv .

#rerooted orthofinder tree to match phylogeny of myrmicines (Pbar as most basal species)
cat SpeciesTree_rerooted_node_labels.txt|perl -pe 's/(....).longestIsoform/$1/g'|sed s/\'//g > speciestree.tre
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1]' Orthogroups_SingleCopyOrthologues.txt Orthogroups.tsv > SCO.tsv

# extract CDS sequences from genbank files with gbseqextractor
## https://pypi.org/project/gbseqextractor/
# -rv reverse if gene on opposite strand

```

## retrieve data for other species
```bash

cd /global/homes/jg/schradel/data/genomes/NCBI/cds/
#search "myrmicinae" in ncbi assemblies. Click "Download assemblies", select "RefSeq" select "CDS from genomic FASTA (.fna)"
ls /global/homes/jg/schradel/data/genomes/NCBI/cds/ncbi-genomes-2020-09-22/
query=$(ls *cds|perl -pe 's/(.*).cds/$1/g'|grep GCF|tr "\n" "|")
esearch -db assembly -query $query | esummary | xtract -pattern DocumentSummary -element RefSeq,SpeciesName > genbank.table.tsv

# select longest isoform from all splice forms
cd $base/data

for file in /global/homes/jg/schradel/data/genomes/NCBI/cds/ncbi-genomes-2020-09-22/*
do
  acc=$(echo $file|perl -pe 's/.*\/(GCF_[0-9]{9}\.[0-9]).*/$1/g')
  short=$(grep $acc genbank.table.tsv|cut -f 2|perl -pe 's/(.).* (...).*/$1$2/g')
  zcat ${file}|perl -pe 's/\>(.*) \[gene=(.*?)\].*/\>$1\t$2/g'|  seqkit fx2tab -l -|tr ";" "\t"|\
  sort -k2,2 -k4,4nr |\
  sort -k2,2 -u -s |\
  awk  -v sho="$short" -v ac="$acc"  '{print ">"$2";species="sho";accession="ac";length="$4"\n"$3}'> ${short}.${acc}.longestIsoform.fa
done

```

## retrieve data for Cobs
```bash

seqkit fx2tab -l /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/Cobs.alpha.v.2.1.geneannotation.1.5.cds.fa|\
sort -k2,2 -k4,4nr |\
sort -k2,2 -u -s |\
perl -pe 's/gene\=//g' |\
awk '{print ">"$2";species=Cobs;accession=cobsalpha1.5;length="$4"\n"$3}'> Cobs.cds.longestIsoform.fa
```

## retrieve fasta for each SCO
```bash
## for each line, retrieve fasta files from all.transcript.fa
cat *.longestIsoform.fa|tr ";" "\t"  > all.cds.fa

cd $base
while read p; do
  OG=$(echo "$p"|cut -f 1)
  echo "$p"|cut -f 2-|perl -pe 's/\t/\n/g'|cut -f 1 -d "-" > $base/OGs/$OG.lst
done <$base/data/SCO.tsv

ls $base/OGs|parallel --nice 10 "seqtk subseq data/all.cds.fa OGs/{} > fa/{}.cds.fa"
```
## align with prank in codon mode
```bash
cd $base
mkdir $base/prank
mkdir $base/sets/

# remove anything but LOCUS names from fasta file
ls fa/*|parallel --nice 10 "perl -pe 's/\>(.*)species=(.*?)\t(.*)/\>"'$1'"/g' {} > sets/{/.}"

# replace species ids with gene LOC ids in tree
## Fancy AWK lookup pattern matching
ls fa/|parallel  --nice 10 "grep \">\" fa/{} |perl -pe 's/\>(.*)\sspecies=(....).*/\$2\t\$1/g' > sets/{=s:.lst.cds.fa:.species2locus.tsv:;=}"

# replace species ids with gene LOC ids in tree
## Fancy AWK lookup pattern matching
ls sets/*species2locus.tsv|parallel --nice 10 "awk -F \$'\t' 'FNR==NR { array[\$1]=\$2; next } { for (i in array) gsub(i, array[i]) }1' {} data/speciestree.tre > {=s:.species2locus.tsv:.tre:;=}"

# prank alignment
ls sets/*.cds|perl -pe 's/.*(OG[0-9]{7}).*/$1/g'|parallel --nice 10 prank -t=sets/{}.tre -d=sets/{}.lst.cds -o=sets/{}.aln -quiet -translate

# trim and remove terminal stops
ls sets/*best.pep.fas| parallel --nice 10 "pal2nal.pl {} {=s:.lst.aln.best.pep.fas:.lst.cds:;=} -nogap -nomismatch -output fasta > {=s:.lst.aln.best.pep.fas:.aln.trim:;=}"

```

## run Hyphy ABSREL
```bash
# run Hyphy ABSREL
ls sets/*.tre|parallel --nice 10 '/global/projects/programs/source/hyphy-2.5.6/hyphy absrel --alignment {=s:.tre:.aln.trim:;=} --tree {}'  
```
