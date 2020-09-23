#*P. californicus* comparative genomics


# Orthofinder
## Prepare environment
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/orthofinder/
mkdir $base
cd $base
```

## Download data
I manually downloaded genbank files for all available myrmicines from NCBI.

```bash
GCF_000143395.1
GCF_000187915.1
GCF_000188075.2
GCF_000204515.1
GCF_000949405.1
GCF_000956235.1
GCF_001594045.1
GCF_001594055.1
GCF_001594065.1
GCF_001594075.1
GCF_001594115.1
GCF_003070985.1
GCF_003260585.2
```

These files are stored here:
`/global/homes/jg/schradel/data/genomes/NCBI/gbk/*/*.gbff`

## Convert gbk to protein fasta
Extracting the longest protein isoform from the gbk files with `funannotate util gbk2parts`.
```bash
cd $base
mkdir $base/data

for file in $(ls /global/homes/jg/schradel/data/genomes/NCBI/gbk/*/genomic.gbff)
do
  short=$(grep -m1 "ORGANISM" $file|perl -pe 's/.*ORGANISM\s+(.).* (...).*/$1$2/g')
  mkdir data/$short/
  funannotate util gbk2parts --gbk $file -o data/$short/$short

  # select longest isoform from all splice forms
  seqkit fx2tab -l data/$short/$short.proteins.fasta|\
  sort -k2,2 -k4,4nr |\
  sort -k2,2 -u -s |\
  awk '{print ">"$1";symbol="$2";length="$4"\n"$3}'> data/$short/$short.longestIsoform.fa
done


# Also get C. obscurior
mkdir data/Cobs/
# select longest isoform from all splice forms
seqkit fx2tab -l /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/Cobs.alpha.v.2.1.geneannotation.1.5.pep.fa|\
sort -k2,2 -k4,4nr |\
sort -k2,2 -u -s |\
awk '{print ">"$1";symbol="$2";length="$4"\n"$3}'> data/Cobs/Cobs.longestIsoform.fa
```

## Run orthofinder
We ran orthofinder (v.2.4.0) to retrieve orthogroups across all available myrmicines.
```bash

cd $base
mkdir $base/input
cd $base/input
ln -s $base/data/*/*.longestIsoform.fa .
/global/homes/jg/schradel/software/OrthoFinder/orthofinder -f $base/input/ -t 20
#Results: /global/homes/jg/schradel/data/Cardiocondyla/comparativeGenomics/orthofinder/OrthoFinder/Results_Sep18

# remove Pcal
cd $base
mkdir $base/input2
cd $base/input2
ln -s $base/data/*/*.longestIsoform.fa .
rm Pcal.longestIsoform.fa
cd $base
/global/homes/jg/schradel/software/OrthoFinder/orthofinder -f $base/input2/ -t 20
mv $base/input2/OrthoFinder/Results_Sep21/ $base/OrthoFinder

mv $base/input/OrthoFinder/ $base
cd $base/OrthoFinder/Results_*
tar -zcvf Gene_Trees.tar.gz Gene_Trees
```
