


```bash
cd /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/analyses/GypsyAnalyses
grep "LTR/Gypsy" /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results/CobsTEsorted.csv.gff |perl -pe 's/.*Alias=(.*?)\;.*/$1\#/g'|sort|uniq -c|perl -pe 's/ +([0-9]+) (.*)/$1\t$2/g'|sort -k 1,1 -nr > Cobs.Gypsy.tsv
awk '{if ($1>10) print $2}' Cobs.Gypsy.tsv > Cobs.Gypsy.10.lst
seqkit grep -r -n -f Cobs.Gypsy.10.lst /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/LibraryPrep/FinalLibraryCobs.fasta > Cobs.Gypsy.10.fa
seqkit stat Cobs.Gypsy.10.fa

bioawk -c fastx '{ if(length($seq) > 1000) { print ">"$name; print $seq }}' Cobs.Gypsy.10.fa > Cobs.Gypsy.10.1kb.fa
bioawk -c fastx '{ if(length($seq) > 1000) { print $name"\t"length($seq) }}' Cobs.Gypsy.10.1kb.fa > Cobs.Gypsy.10.1kb.length.tsv

mafft Cobs.Gypsy.10.1kb.fa  > Cobs.Gypsy.10.1kb.mafft.aln
FastTree Cobs.Gypsy.10.1kb.mafft.aln > Cobs.Gypsy.10.1kb.mafft.FT.tre
cat Cobs.Gypsy.10.lst|perl -pe 's/#/;/g' > Cobs.Gypsy.10.lst2

grep -f Cobs.Gypsy.10.lst2  /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results/CobsTEsorted.csv.gff|perl -pe 's/(.*)\t.*?Alias=(.*?);.*/$1\t$2/g' > Cobs.Gypsy.10.gff
```

```bash
cd /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/analyses/GypsyAnalyses
grep "LTR/Gypsy" /Users/lukas/sciebo/Projects/TEannotation/RepBase25.04.fasta/invrep.ref
awk '{if ($1>10) print $2}' Cobs.Gypsy.tsv > Cobs.Gypsy.10.lst
seqkit grep -r -n -f Cobs.Gypsy.10.lst /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/LibraryPrep/FinalLibraryCobs.fasta > Cobs.Gypsy.10.fa
seqkit stat Cobs.Gypsy.10.fa

bioawk -c fastx '{ if(length($seq) > 1000) { print ">"$name; print $seq }}' Cobs.Gypsy.10.fa > Cobs.Gypsy.10.1kb.fa
bioawk -c fastx '{ if(length($seq) > 1000) { print $name"\t"length($seq) }}' Cobs.Gypsy.10.1kb.fa > Cobs.Gypsy.10.1kb.length.tsv

mafft Cobs.Gypsy.10.1kb.fa  > Cobs.Gypsy.10.1kb.mafft.aln
FastTree Cobs.Gypsy.10.1kb.mafft.aln > Cobs.Gypsy.10.1kb.mafft.FT.tre
cat Cobs.Gypsy.10.lst|perl -pe 's/#/;/g' > Cobs.Gypsy.10.lst2

grep -f Cobs.Gypsy.10.lst2  /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results/CobsTEsorted.csv.gff|perl -pe 's/(.*)\t.*?Alias=(.*?);.*/$1\t$2/g' > Cobs.Gypsy.10.gff
```


```bash
bioawk -c fastx '{ if(length($seq) > 10000) { print ">"$name; print $seq }}' Cobs.Gypsy.10.fa > Cobs.Gypsy.10.10kb.fa
bioawk -c fastx '{ if(length($seq) > 10000) { print $name"\t"length($seq) }}' Cobs.Gypsy.10.10kb.fa > Cobs.Gypsy.10.10kb.length.tsv

# get Gypsy-16_PBa-I
seqkit grep -r -n -p Gypsy-16_PBa-I  /Users/lukas/sciebo/Projects/TEannotation/RepBase25.04.fasta/invrep.ref > Gypsy-16_PBa-I.fa
seqkit grep -r -n -p "CobsR.176#"  /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/LibraryPrep/FinalLibraryCobs.fasta > CobsR176.fa
cat Gypsy-16_PBa-I.fa Cobs.Gypsy.10.1kb.fa  > gypsy.combined.fa
cat Gypsy-16_PBa-I.fa Cobs.Gypsy.10.10kb.fa  > gypsy.combined2.fa
mafft gypsy.combined.fa  > gypsy.combined.aln
mafft gypsy.combined2.fa  > gypsy.combined2.aln
FastTree -nt gypsy.combined.aln > gypsy.combined.aln.FT.tre
FastTree -nt gypsy.combined.aln > gypsy.combined2.aln.FT.tre

cat CobsR176.fa Gypsy-16_PBa-I.fa |mafft - > CobsR176.Gypsy-16_PBa-I.aln
```

```bash
## LTRharvest
```bash
cd $base

ln -s Gypsy-16_PBa-I.fa Gypsy-16_PBa-I.reference.fa
TE=Gypsy-16_PBa-I
gt suffixerator -db gypsy.combined.fa -indexname gypsy.combined.fa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index gypsy.combined.fa  -gff3 gypsy.combined.ltrharvest.gff -out gypsy.combined.ltrharvest.fa -seqids yes -tabout no > gypsy.combined.ltrharvest.out
gt ltrharvest -index $TE.reference.fa  -gff3 $TE.reference.ltrharvest2.gff -out $TE.reference.ltrharvest2.fa -seqids no -tabout no > $TE.ltrharvest.2.out
```

```
