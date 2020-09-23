---
title: 'Cobs Alpha minION assembly: remove contaminating scaffolds'
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 2
    collapsed: true
    number_sections: true
    df_print: paged
  pdf_document: default
---
<STYLE TYPE="text/css">
<!--
  td{
    font-family: Arial;
    font-size: 8pt;
    padding:0px;
    cellpadding="0";
    cellspacing="0"
  }
  th {
    font-family: Arial;
    font-size: 8pt;
    height: 20px;
    font-weight: bold;
    text-align: right;
    background-color: #ccccff;
  }
  table {
    border-spacing: 0px;
    border-collapse: collapse;
  }
--->
</STYLE>

# Screen alpha Cobs2.0 assembly for contaminations

This workflow will identify bacterial and mitochondrial scaffolds and remove them from the assembly v2.0 of *Cardiocondyla obscurior* alpha.

### INPUT
**Cobs.alpha version 2.0**
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa```**

### OUTPUT
**Cobs.alpha version 2.1**
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.1.fa```**

**Cobs.alpha mitochondrial scaffold**
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.0.mitochondrial.fa```**

**Cobs.alpha bacterial scaffolds *Westeberhardia* and *Wolbachia***
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.0.bacterial.fa```**

## Download set of genomes
I first download 2422 chromosome-level assemblies of prokaryotes from NCBI. These will be used as the bacterial database.
```bash
# Download list of complete bacterial genome sequences
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

# retrieve unique subgroup (field 7)(or group (field 6))
cat ./prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'|awk -F '\t' '!seen[$7]++'|wc -l
# % 2422 genomes
query=$(cat ./prokaryotes.txt|awk -F '\t' '{if ($16=="Chromosome") print $0}'|awk -F '\t' '!seen[$7]++'|cut -f 19|sed 1d|tr "\n" " "|perl -pe 's/ / OR /g'|perl -pe "s/ OR $//g")

esearch -db assembly  -query "$query" |elink -target nuccore|efetch -format fasta > db.2422genomes.fa

```
Format the 2422 genomes as a blast DB.
```bash
#format db
cd /global/scratch/schradel/BLAST/prokScreen/
makeblastdb -in db.2422genomes.fa -title db -dbtype nucl
```
#### Database prepared at
**```/global/scratch/schradel/BLAST/prokScreen/db.2422genomes.fa```**

## Blast assembly windows against bacterial DB
```bash
base=/global/scratch/schradel/BLAST/prokScreen.alphaCobs
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa
```
1. Setup environment
2. Create windows
3. Blast each window against the database using **megablast**
4. Delete empty files
```bash
#1. make folders
mkdir windows
mkdir bls

#prepare genome fasta
cut -f 1 -d "|" $genome > genome.fa
samtools faidx genome.fa

#2. calculate windows (store in windows folder)
cut -f 1,2 genome.fa.fai > genome.lengths.tsv
bedtools makewindows -g genome.lengths.tsv -w 2500  > genome.windows.tsv
#grep scaffold42 genome.windows.tsv > scf42.windows.tsv

#3. run blast for all windows
while read p; do
  scf=$(echo $p|cut -f 1 -d " ")
  start=$(echo $p|cut -f 2 -d " ")
  stop=$(echo $p|cut -f 3 -d " ")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  samtools faidx genome.fa $scf:$start-$stop > windows/$scf.$startL.$stopL.fa
  blastn -task megablast -query windows/$scf.$startL.$stopL.fa -db /global/scratch/schradel/BLAST/dbs/db.2422genomes.fa -num_threads 10 -evalue 1e-3 -num_alignments 1  -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore length pident sgi sacc stitle' > bls/$scf.$startL.$stopL.bls
done < genome.windows.tsv
# takes ~3 h

#4. delete empty files
find ./bls/ -empty -delete
```
# Blast assembly windows against ant DB
## Create ant database
Retrieve all available ant genome assemblies from NCBI (#48 in November 2019)

```bash
esearch -db assembly  -query "Formicidae"|elink -target nuccore|efetch -format fasta > db.antGenomes.fa
makeblastdb -in db.antGenomes.fa -title db -dbtype nucl
```

### Blast agains antDB

**```/global/scratch/schradel/BLAST/dbs/db.antGenomes.fa```**

Blast only those windows that had a significant hit in the prokaryotic DB
```bash
mkdir blsAnt

while read p; do
  window=$(echo $p|cut -f 1 -d " ")
  scf=$(echo $window|cut -f 1 -d ":")
  start=$(echo $window|cut -f 2 -d ":"|cut -f 1 -d "-")
  stop=$(echo $window|cut -f 2 -d ":"|cut -f 2 -d "-")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  echo $scf $start $stop
  blastn -task megablast -query windows/$scf.$startL.$stopL.fa -db /global/scratch/schradel/BLAST/dbs/db.antGenomes.fa -num_threads 10 -evalue 1e-3 -num_alignments 1  -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore length pident sgi sacc stitle' > blsAnt/$scf.$startL.$stopL.bls
done < genome.vs.2422genomes.bls

find ./blsAnt/ -empty -delete

```

### Retrieve all significant hits
```bash
head -n 1 -q blsAnt/*  > genome.vs.Antgenomes.bls
```
# Blast assembly windows against well assembled insect DB
**```/global/scratch/schradel/BLAST/dbs/db.insectGenomes.fa```**
### Prepare DB of well assembled insect genomes
I retrieved **insect** genomes that were assembled to **Chromosome** or **Complete Genome** level from NCBI.

```bash
# get list of all genomes from NCBI (~30.11.2019)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt
# get unique insect genera that were fully assembled
query=$(cat ./eukaryotes.txt|awk -F '\t' '{if ($17=="Chromosome" || $17=="Complete Genome") print $0}'|awk -F ' ' '!seen[$1]++'|awk -F '\t' '{if ($6=="Insects") print $0}'|cut -f 9|sed 1d|tr "\n" " "|perl -pe 's/ / OR /g'|perl -pe "s/ OR $//g")
# retrieve fasta of genomes
esearch -db assembly  -query "$query"|elink -target nuccore|efetch -format fasta > db.insectGenomes.fa
# make blast DB
makeblastdb -in db.insectGenomes.fa -title db -dbtype nucl
```
### Blast windows against insectDB
```bash

cd /global/scratch/schradel/BLAST/dbs
makeblastdb -in db.insectGenomes.fa -title db -dbtype nucl
cd ../prokScreen.alphaCobs/
mkdir blsInsect
while read p; do
  window=$(echo $p|cut -f 1 -d " ")
  scf=$(echo $window|cut -f 1 -d ":")
  start=$(echo $window|cut -f 2 -d ":"|cut -f 1 -d "-")
  stop=$(echo $window|cut -f 2 -d ":"|cut -f 2 -d "-")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  echo $scf $start $stop
  blastn -task megablast -query windows/$scf.$startL.$stopL.fa -db /global/scratch/schradel/BLAST/dbs/db.insectGenomes.fa -num_threads 10 -evalue 1e-3 -num_alignments 1  -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore length pident sgi sacc stitle' > blsInsect/$scf.$startL.$stopL.bls
done < genome.vs.2422genomes.bls

find ./blsInsect/ -empty -delete

```

# rRNA prediction
A lot of the hits against the bacterial genomes are in fact rRNAs. These are so conserved and eukaryotic rRNAs look very much like bacterial rRNAs, so that the blastn search against the prokaryotic DB produces a lot of significant hits. I use barrnap (https://github.com/tseemann/barrnap) to identify pro- and eukaryotic rRNAs.

###rRNA prediction with barrnap
```bash
base=/global/scratch/schradel/BLAST/prokScreen.alphaCobs
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa

mkdir rRNA
mkdir prok.rRNA

while read p; do
  window=$(echo $p|cut -f 1 -d " ")
  scf=$(echo $window|cut -f 1 -d ":")
  start=$(echo $window|cut -f 2 -d ":"|cut -f 1 -d "-")
  stop=$(echo $window|cut -f 2 -d ":"|cut -f 2 -d "-")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  echo $scf $start $stop
  barrnap -q -k euk windows/$scf.$startL.$stopL.fa > rRNA/$scf.$startL.$stopL.rRNA.gff3
done < genome.vs.2422genomes.bls

cat rRNA/*|egrep "^#" -v > euk.rRNA.gff3

while read p; do
  window=$(echo $p|cut -f 1 -d " ")
  scf=$(echo $window|cut -f 1 -d ":")
  start=$(echo $window|cut -f 2 -d ":"|cut -f 1 -d "-")
  stop=$(echo $window|cut -f 2 -d ":"|cut -f 2 -d "-")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  echo $scf $start $stop
  barrnap -q -k bac windows/$scf.$startL.$stopL.fa > prok.rRNA/$scf.$startL.$stopL.prok.rRNA.gff3
done < genome.vs.2422genomes.bls


cat prok.rRNA/*|egrep "^#" -v > pro.rRNA.gff3


join  -t $'\t' -j 1 <(sort -k1 euk.rRNA.gff3) <(sort -k1 pro.rRNA.gff3) > rRNA.gff3
```

# Find mitochondrial scaffolds
I downloaded the C. obscurior published mitochndrium and blasted it against the genome to identify mitochondrial scaffolds in the assembly.
```bash
base=/global/scratch/schradel/BLAST/prokScreen.alphaCobs
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa

# Cardiocondyla mitochondrium
esearch -db nucleotide -query "KX951753.1"| efetch -format fasta > ../dbs/db.KX951753.1.fa

# blast mitochndrium against genome assembly
mkdir blsMito
makeblastdb -in ../dbs/db.KX951753.1.fa -dbtype nucl

while read p; do
  scf=$(echo $p|cut -f 1 -d " ")
  start=$(echo $p|cut -f 2 -d " ")
  stop=$(echo $p|cut -f 3 -d " ")
  startL=$(printf "%08d" $start )
  stopL=$(printf "%08d" $stop )
  blastn -task megablast -query windows/$scf.$startL.$stopL.fa -db /global/scratch/schradel/BLAST/dbs/db.KX951753.1.fa -num_threads 6 -evalue 1e-10 -num_alignments 1  -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore length pident sgi sacc stitle' > blsMito/$scf.$startL.$stopL.bls
done < genome.windows.tsv

# blast mitochondtium against genome
makeblastdb -in /global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa -dbtype nucl
blastn -query /global/scratch/schradel/BLAST/dbs/db.KX951753.1.fa -db /global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa -outfmt '6 qseqid qstart qend sseqid sstart send evalue bitscore length pident sgi sacc stitle' > Cobs2.0.vs.CobsMito.bls

```


## Calculate GC content and length for each scaffold
```bash
# calculate gc and length
infoseq  -nocolumn -delimiter "\t" -auto -only  -name -length -pgc genome.fa >genome.GC.tsv
```



# Final results
The analysis shows that four scaffolds are bacterial:

* scaffold36 Westeberhardia
* scaffold37 Wolbachia
* scaffold51 Wolbachia
* scaffold83 Wolbachia

These will be removed from the assembly, thus generating CobsAlpha.v2.1


Scaffold106 appears to contain the mitochondrium twice. The coverage from short read data is >2400X according to QualiMap.
* scaffold106 mitochondrial

# Compare, annotate, and analyse in R
[R script to analyze these data](alphaCobs.screenAssembly.Rmd):
[**```/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/alphaCobs.analyse.screenAssembly.Rmd```**](/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/alphaCobs.analyse.screenAssembly.Rmd)


The analysis shows that four scaffolds are bacterial:

* **scaffold36  Westeberhardia**
* **scaffold37 Wolbachia**
* **scaffold51 Wolbachia**
* **scaffold83 Wolbachia**
* **scaffold106 Mitochondrial**

These will be removed from the assembly, thus generating CobsAlpha.v2.1

```{bash, eval=FALSE}

base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter

genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa

echo -e "scaffold36$\nscaffold37$\nscaffold51$\nscaffold83$" >contaminants.lst

seqkit grep -r -f contaminants.lst $genome > Cobs.alpha.2.0.bacterial.fa
seqkit grep -r -p "scaffold106$" $genome > Cobs.alpha.2.0.mitochondrial.fa
seqkit grep -rv -f contaminants.lst $genome > Cobs.alpha.2.05.fa
seqkit grep -rv -p "scaffold106$" Cobs.alpha.2.05.fa > Cobs.alpha.2.1.fa


```
![](assets/alphaCobs.screenAssembly-ed75db38.pdf)

# QC to get final numbers
## Run QUAST
```bash
cd /global/homes/jg/schradel/data/assemblies/Cobs.alpha/QC/QUAST
# named a11
a11=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.1.fa
nice quast.py $a11 -o a11.quast.out --eukaryote -s
```
[QUAST results for Cobs.alpha.v2.1](./results/quast/a11.quast.out/report.tsv)
```/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/quast/a11.quast.out/report.tsv```
| Assembly | Cobs.alpha.2.1 | Cobs.alpha.2.1_broken |
|--|--|--|
| # contigs (>= 0 bp) | 127 | - |
| # contigs (>= 1000 bp) | 127 | 188 |
| # contigs (>= 5000 bp) | 125 | 184 |
| # contigs (>= 10000 bp) | 121 | 177 |
| # contigs (>= 25000 bp) | 105 | 150 |
| # contigs (>= 50000 bp) | 91 | 112 |
| Total length (>= 0 bp) | 193051228 | - |
| Total length (>= 1000 bp) | 193051228 | 192868330 |
| Total length (>= 5000 bp) | 193047645 | 192858934 |
| Total length (>= 10000 bp) | 193025568 | 192812504 |
| Total length (>= 25000 bp) | 192755434 | 192313565 |
| Total length (>= 50000 bp) | 192237042 | 190916813 |
| # contigs | 127 | 188 |
| Largest contig | 13148674 | 12360777 |
| Total length | 193051228 | 192868330 |
| GC (%) | 41.02 | 41.02 |
| N50 | 6290588 | 5730056 |
| N75 | 4487289 | 4157591 |
| L50 | 11 | 13 |
| L75 | 21 | 23 |
| # N's per 100 kbp | 94.76 | 0.02 |


## Run BUSCO
```bash
cd /global/homes/jg/schradel/data/assemblies/Cobs.alpha/QC/BUSCO
export PATH=$PATH:"/global/projects/programs/source/augustus-3.3.1/bin/"
#cp -r /global/projects/programs/source/augustus-3.3.1/config/ .
export AUGUSTUS_CONFIG_PATH=$(readlink -f ./config/)
BUSCO_CONFIG_FILE=/global/projects/programs/source/busco/config/config.ini

a11=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/prokFilter/Cobs.alpha.2.1.fa
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a11 --out a11 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

```
