# TE phylogenies

Here is how I analysed TE sequences from the RepeatMasker+one-code-to-find-them-all C. obscurior alpha v.2.1 assembly annotations.

This script will
- extract TE sequences of particularly intersting TEs (mostly Gypsys, I guess).
- align these sequences
- Compute a phylogeny
- Create blast-based overviews of the coverage of the match relative to the reference.

## Extract repeats from RepeatMasker annotation

```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEphylogenies/
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
annotation=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/CobsTEsorted.csv

cd $base

# LTRs
TE="CobsR.176"
TE="CobsR.103"
TE="CobsR.1222"
TE="CobsR.23"
#DNA transposons
TE="CobsR.230"


cat $annotation|awk '{if ($8=="'"$TE"'") print $4"\t"$5"\t"$6"\t"$8"."++count[$8]"\t("$0")"}' > $TE.bed
bedtools getfasta -fi $genome  -bed $TE.bed -name > $TE.fa
```

## Align sequences
```bash
#with prank
#prank $TE.fa #killed
#with mafft
/global/homes/jg/schradel/conda/bin/mafft $TE.fa > $TE.mafft.aln
#cut -f 1 -d ":" $TE.mafft.aln > tmp.$TE.aln
```
## Phylogeny
```bash
#with FastTree
FastTreeMP $TE.mafft.aln > $TE.mafft.tre
```

## Compare TEs to reference
```bash

# get TE reference from library
TElib=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/LibraryPrep/FinalLibraryCobs.fasta

bioawk  -c fastx '$name ~ '/$TE#/' { print ">"$name"\n"$seq"\n"; }' $TElib > $TE.reference.fa
makeblastdb -in $TE.reference.fa -dbtype nucl
blastn -task megablast -query $TE.fa -db $TE.reference.fa -outfmt 6 > $TE.bls
#blastn -task megablast -query $TE.fa -db $TE.reference.fa -outfmt 7 > $TE.rawbls
#blastn -task megablast -query $TE.fa -db $TE.reference.fa -outfmt 17 > $TE.sam

```

## Annotate reference TE
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEphylogenies/
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
annotation=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/CobsTEsorted.csv
TElib=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/LibraryPrep/FinalLibraryCobs.fasta

cd /home/s/schradel/dbs/interpro/interproscan-5.44-79.0
mkdir $base/interpro/
./interproscan.sh \
--input $base/$TE.reference.fa \
--disable-precalc \
--output-dir $base/interpro/ \
--formats TSV,XML,GFF3 \
--goterms \
-t n

cd $base/interpro
```

## LTRharvest
```bash
cd $base
TE=CobsR.176
gt suffixerator -db $TE.reference.fa -indexname $TE.reference.fa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index $TE.reference.fa  -gff3 $TE.reference.ltrharvest.gff -out $TE.reference.ltrharvest.fa -seqids yes -tabout no > $TE.ltrharvest.out
gt ltrharvest -index $TE.reference.fa  -gff3 $TE.reference.ltrharvest2.gff -out $TE.reference.ltrharvest2.fa -seqids no -tabout no > $TE.ltrharvest.2.out
```

# Download data for plotting in R
```bash
cd /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/landscapes/phylogenies/
TE="CobsR.103"
TE="CobsR.1222"
mkdir $TE
mkdir $TE/interpro/

scp jgpogo:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEphylogenies/*$TE* ./$TE/
scp jgpogo:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEphylogenies/interpro/*$TE* ./$TE/interpro/

#Run /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotPhylo.Rmd
```
