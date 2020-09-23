# Process gene annotation for C. obscurior alpha v.2.1 assembly

Here, we will process and filter the gene annotations produced by BRAKER and GeMoMa.


For this, we will

1. Convert the raw gff file to proper gff3 format
2. Functionally annotate all proteins with interproscan
3. Subset the gff3 to only retain the longest isoform per gene
4. Flag putative TE-encoded genes by
   4.1. Blast proteins against RepBase
   4.2. TransposonPSI analysis
   4.3. Identify CDS overlapping annotated TEs
   4.4. Screening functional annotations for TE-related terms
5. Setup webapollo track with latest annotation


## 1. Convert raw gff file to proper gff3 format

This is run on pallas for simplicity.

####Upload data
```bash
scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.geneCombined/Cobs.alpha.v.2.1.combined_6ref_braker_v1.3.gff lschrader@pallas.bio.ku.dk:/usr/local/home/lschrader/data/IBE/Cardiocondyla/annotation/
```

#### Setup environment
```bash
export EVMutils=/usr/local/home/lschrader/software/EVidenceModeler-1.1.1/EvmUtils
export target_genome=/corefac/cse/lukas/IBE/Cardiocondyla/CardioMinION/results/Cobs.alpha.2.1.fa
export data=/corefac/cse/lukas/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.combined_6ref_braker_v1.3.gff
```
#### install tools
```bash
## GFFUtils (see https://github.com/fls-bioinformatics-core/GFFUtils and https://github.com/fls-bioinformatics-core/GFFUtils/issues/14)
## genometools (gt) http://genometools.org/
## gffread http://ccb.jhu.edu/software/stringtie/gff.shtml
```

```bash
# navigate
cd ~/data/IBE/Cardiocondyla/annotation

# Cleanup gff output files
GFFcleaner \
--clean-replace-attributes \
--add-missing-ids \
--add-exon-ids \
--report-duplicates \
--insert-missing= $data \
-o tmp.clean.gff 2>/dev/null

gt gff3 \
-sort \
-tidy  \
-checkids yes \
-retainids yes \
-fixregionboundaries tmp.clean.gff | \
sed s/GAF/GeMoMa/g|sed s/\+AnnotationEvidence//g> tmp.gff 2>/dev/null

#gt gff3 -sort -tidy  -addids -fixregionboundaries -checkids yes -retainids yes tmp.gff > tmp2.gff 2>/dev/null

awk '{if ( $3=="gene") print $0 }' tmp.gff > tmp.genes.gff 2>/dev/null

gffread tmp.gff -F -o tmp3.gff --force-exons -E 2>/dev/null

cat tmp.genes.gff tmp3.gff | gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -|egrep "^###" -v|sed s/geneID/Parent/g > tmp4.gff

GFFcleaner \
--clean-replace-attributes \
--add-missing-ids \
--add-exon-ids \
--report-duplicates \
--insert-missing= tmp4.gff \
-o tmp5.gff 2>/dev/null


# define output file name
export output=Cobs.alpha.v.2.1.geneannotation.1.3
gt gff3 -retainids -sort tmp5.gff > $output.gff3

# check that files are proper gff3
perl $EVMutils/gff3_gene_prediction_file_validator.pl $output.gff3

# export protein fasta sequences
perl $EVMutils/gff3_file_to_proteins.pl $output.gff3 $target_genome |sed $'s/[^[:print:]\t]/ /g' > $output.pep.fa 2>/dev/null &


# Cleanup gff output files: webapollo friendly
gffread tmp.gff -o tmp3.wa.gff --force-exons -E 2>/dev/null
cat tmp.genes.gff tmp3.wa.gff |gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -|egrep "^###" -v|sed s/geneID/Parent/g > tmp4.wa.gff
GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= tmp4.wa.gff -o tmp5.wa.gff 2>/dev/null
gt gff3 -retainids -sort tmp5.wa.gff > $output.webapolloFriendly.gff3

# cleanup
mv tmp* tmp/

```

####Download data
```bash
scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" lschrader@pallas.bio.ku.dk:/usr/local/home/lschrader/data/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.geneannotation.1.3.* /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.geneCombined/
```


# 2. interproScan for Cardiocondyla v2.1 gene annotation

```bash
#@jgpogo
cd /home/s/schradel/dbs/interpro/interproscan-5.44-79.0
cat /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.v.2.1.geneannotation.1.3.pep.fa|perl -pe 's/\*//g' > Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa

./interproscan.sh \
--input Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa \
--disable-precalc \
--output-dir /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro \
--formats TSV,XML,GFF3 \
--goterms

cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro

```

### Katharina Hoff's predictions
```bash
cd ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.functionalAnnotation
mv augustus.hints.aa.tsv augustus_Braker.aa.tsv

for f in *.tsv
do
  tag=$(echo $f|perl -pe 's/.*_(.*?)\..*/$1/g')
  #echo $tag
  cat $f|perl -pe "s/^/${tag}_/g"
done > interproscan.combined.aa

scp interproscan.combined.aa jgant2:/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/
```

### check functional annotations for TE-related terms
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro
cd $base
cat $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv|cut -f 5,6|sort |uniq -c|perl -pe 's/^ *([0-9]+) /$1\t/g'|sort -k 1,1 -nr > $base/Cobs.alpha.v.2.1.geneannotation.1.3.functionalAnnotations.frequency
cat $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv|cut -f 12,13|sort |uniq -c|perl -pe 's/^ *([0-9]+) /$1\t/g'|sort -k 1,1 -nr > $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproAnnotations.frequency
```
I manually screened all functional annotations and their frequency and compiled a list of 40 functional annotation terms that are TE associated. A list of all these functional annotation terms is stored in:
```/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.1.lst```

```
G3DSA:3.30.70.270	Reverse transcriptase/Diguanylate cyclase domain
SSF53098	"Ribonuclease	H-like superfamily"
PF00078	Reverse transcriptase (RNA-dependent DNA polymerase)
G3DSA:3.30.420.10	Ribonuclease H-like superfamily/Ribonuclease H
PS50878	Reverse transcriptase (RT) catalytic domain profile.
PS50994	Integrase catalytic domain profile.
G3DSA:3.10.10.10	HIV Type 1 Reverse Transcriptase, subunit A, domain 1
PF00665	Integrase core domain
G3DSA:3.10.20.370	retrotransposable element/transposon Tf2-type
PF17921	Integrase zinc binding domain
cd01647	RT_LTR
SSF57756	Retrovirus zinc finger-like domains superfamily
cd09274	RNase_HI_RT_Ty3
G3DSA:1.10.340.70	Ribonuclease H superfamily
PF17919	RNase H-like domain found in reverse transcriptase
cd00303	retropepsin_like
G3DSA:3.30.420.470	transposase type 1
cd01650	RT_nLTR_like
PF14529	Endonuclease-reverse transcriptase
PF05380	Pao retrotransposon peptidase
PF17917	RNase H-like domain found in reverse transcriptase
PF12259	Baculovirus F protein
PS50175	Aspartyl protease, retroviral-type family profile.
PF00077	Retroviral aspartyl protease
cd09077	R1-I-EN
PF07727	Reverse transcriptase (RNA-dependent DNA polymerase)
PF14223	gag-polypeptide of LTR copia-type
PF03732	Retrotransposon gag protein
PF04665	Poxvirus A32 protein
cd09272	RNase_HI_RT_Ty1
PF13359	DDE superfamily endonuclease
PF01498	Transposase
PF13976	GAG-pre-integrase domain
PF13975	gag-polyprotein putative aspartyl protease
PF14214	Helitron helicase-like domain at N-terminus
PF05699	hAT family C-terminal dimerisation region
cd09275	RNase_HI_RT_DIRS1
PF12017	Transposase protein
PF14787	GAG-polyprotein viral zinc-finger
cd09276	Rnase_HI_RT_non_LTR
PF01359	Transposase (partial DDE domain)

```

```bash
grep -f /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.1.lst $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv |cut -f 1 |sort|uniq|wc -l
```
In ```Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv```, there are 3868 proteins annotated with at least one of the above terms.

In addition, I screened interpro annotations.

```
IPR000477	Reverse transcriptase domain
IPR012337	Ribonuclease H-like superfamily
IPR036397	Ribonuclease H superfamily
IPR041577	Reverse transcriptase/retrotransposon-derived protein, RNase H-like domain
IPR008042	Retrotransposon, Pao
IPR041373	Reverse transcriptase, RNase H-like domain
IPR022048	Envelope fusion protein-like
IPR001995	Peptidase A2A, retrovirus, catalytic
IPR013103	Reverse transcriptase, RNA-dependent DNA polymerase
IPR005162	Retrotransposon gag domain
IPR038717	Tc1-like transposase, DDE domain
IPR027806	Harbinger transposase-derived nuclease domain
IPR002492	Transposase, Tc1-like
IPR002156	Ribonuclease H domain
IPR025724	GAG-pre-integrase domain
IPR025476	Helitron helicase-like domain
IPR008906	HAT, C-terminal dimerisation domain
IPR021896	Transposase protein
IPR004211	Recombination endonuclease VII
IPR010998	Integrase/recombinase, N-terminal
IPR041426	Mos1 transposase, HTH domain
IPR034132	Retropepsin Saci-like domain
IPR026103	Harbinger transposase-derived nuclease, animal
IPR024445	ISXO2-like transposase domain
IPR029526	PiggyBac transposable element-derived protein
IPR018289	MULE transposase domain
IPR029472	Retrotransposon Copia-like, N-terminal
IPR025898	Tc3 transposase, DNA binding domain
```

A list of all interpro annotation terms is stored in:
```/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.2.lst```

```bash
grep -f /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.2.lst $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv |cut -f 1 |sort|uniq|wc -l
```

In ```Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv```, there are 3004 proteins annotated with at least one of the above interpro terms.

### Combine both annotation term lists and screen
```bash
cat /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.1.lst /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.2.lst > /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.combined.lst

grep -f /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.combined.lst $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv|cut -f 1 |sort|uniq|wc -l
#3954

grep -f /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/TEterms.combined.lst $base/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv|cut -f 1 |sort|uniq > $base/Cobs.alpha.v.2.1.geneannotation.1.3.TE
```

Together, 3954 genes are identified.

There remain plenty of very likely TE-derived genes with functional annotations that do occur in TEs but could also characterize proper ant genes.

Best choice will be to use the following logic:
- proteins that have a TE functional annotation (based on the above lists), overlap predicted TEs and have blast hits against RepBase proteins are **TEderived=3**

- proteins that have a TE functional annotation (based on the above lists) and have blast hits against RepBase proteins are **TEderived=2**
- proteins that have no functional annotation (except MobiDBLite) and a blast hit against RepBase proteins are **TEderived=2**
- proteins that have a blast hit against RepBase proteins are **TEderived=1**
- proteins that have a TE functional annotation (based on the above lists) are **TEderived=1**
- proteins that overlap an annotated TE **TEderived=1**


#3. Subset the gff3 to only retain the longest isoform per gene
work on jgant2, as here conda and agat are installed

#### Upload data
```bash
ssh jgant2
scp -o ProxyCommand="ssh -W %h:%p zcm795@ssh-bio-labnet.science.ku.dk" lschrader@pallas.bio.ku.dk:/usr/local/home/lschrader/data/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.geneannotation.1.3.* /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/
```

#### Select longest isoform with agat
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/
# alternative: agat! https://github.com/NBISweden/AGAT/wiki
agat_sp_keep_longest_isoform.pl -gff Cobs.alpha.v.2.1.geneannotation.1.3.gff3 > Cobs.alpha.v.2.1.geneannotation.1.3.longestIsoform.gff3
```
-------------------------
#4. Flag putative TE-encoded genes by

4.1. Blast proteins against RepBase
4.2. TransposonPSI analysis
4.3. Identify CDS overlapping annotated TEs
4.4. Screening functional annotations for TE-related terms

-------------------------

## 4.1. blast against RepBase library

### Setup variables
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/
repbaseDB=/global/homes/jg/schradel/data/repbase20.05_aaSeq_cleaned_TE.fa
fastaData=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.v.2.1.geneannotation.1.3.pep.fa
outname=Cobs.alpha.v.2.1.geneannotation.1.3
```
### Setup environment
```bash

mkdir $base/data
mkdir $base/RepBaseBlast
makeblastdb -in $repbaseDB -dbtype prot
```
### Run blast vs RepBase20.05
```bash
#blastp
cd $base/RepBaseBlast
cat $fastaData | \
parallel \
--block 10k \
--recstart '>' \
--pipe blastp \
-evalue 1e-10 \
-outfmt 6 \
-db $repbaseDB \
-query - > $outname.repbase.bls

cat $outname.repbase.bls|sort -gk11,11 |sort -buk 1,1  > $outname.repbase.bls.bls.bh

```
8156 proteins are identified as homologous to RepBase proteins. These proteins will be flagged.

### 4.2. TransposonPSI analysis

TransposonPSI is run on pallas, as it was installed there already.

### Run transposon psi to identify potential TE-derived or TE-contaminated genes.
```bash
#edited the TransposonPSI.pl  to use 20 threads in blastal˘l (-a 20). Uses blastgpg though unless nuc database
base=/usr/local/home/lschrader/data/IBE/Cardiocondyla/TransposonPSI/
```

```bash
# run TransposonPSI
cd $base

# Cobs cleaned up v2.1 gene annotation v1.2 (same as v1.3 but gff formatting is slightly different)
data=/corefac/cse/lukas/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.geneannotation.1.2.fa
perl ~/software/TransposonPSI_08222010/transposonPSI.pl $data prot 1>Cobs.alpha.v.2.1.geneannotation.1.2.PSI.out 2> Cobs.alpha.v.2.1.geneannotation.1.2.TPSI.err

# Cobs cleaned up v2.1 gene annotation v1.3
data=/corefac/cse/lukas/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.geneannotation.1.3.pep.fa
perl ~/software/TransposonPSI_08222010/transposonPSI.pl $data prot 1>Cobs.alpha.v.2.1.geneannotation.1.3.PSI.out 2> Cobs.alpha.v.2.1.geneannotation.1.3.TPSI.err

```

### 4.3. Identify CDS overlapping annotated TEs

Find overlap between CDS and TE annotations

#### upload data
```bash
scp /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v2.1.TE.gff3 jgant2:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/
```

#### setup environment
```bash
# run on jgant2
base=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/overlap
gene=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.v.2.1.geneannotation.1.3.gff3
repeat=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Cobs.alpha.v2.1.TE.gff3
```

```bash
cd $base
awk '{if ($3=="CDS") print $0}' $gene > CDS.gff3

# bedtools intersect with minimal overlap of CDS by 50 %
bedtools intersect -wa -a CDS.gff3 -b $repeat -f 0.5 > CDS.repeat.overlap.gff
bedtools intersect -wo -a CDS.gff3 -b $repeat -f 0.5 > CDS.repeat.overlap.wo.gff
cat CDS.repeat.overlap.gff |perl -pe 's/.*Parent=(.*)/$1/g'|perl -pe 's/;.*//g'|sort|uniq > mRNA.repeat.overlap.lst

# index gff
gffutils-cli create $gene

# split mRNA list into chunks of 2000 and retrieve all parent IDs using gffutils-cli
mkdir splitFiles
mkdir parents
cd splitFiles
split --lines=2000 ../mRNA.repeat.overlap.lst
cd $base
for file in $base/splitFiles/*
do
  set=$(echo $file|perl -pe 's/.*\/(.*?)/$1/g')
  echo $set
  ids=$(cat $file |tr "\n" ","|sed 's/.$//')
  gffutils-cli parents --exclude-self $gene.db $ids |perl -pe 's/.*ID=(.*?);.*|$/$1/g' > parents/$set.parents.lst
done

# file containing all genes that are overlapping TEs
cat parents/* |sort|uniq> gene.repeat.overlap.lst

```

### 4.4 Combine results and update gff
```bash
base=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/combine
gene=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.v.2.1.geneannotation.1.3.gff3
overlap=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/overlap/mRNA.repeat.overlap.lst
repbase=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/RepBaseBlast/Cobs.alpha.v.2.1.geneannotation.1.3.repbase.bls.bls.bh
funan=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/Cobs.alpha.v.2.1.geneannotation.1.3.TEcandidates.lst

# OVERLAP
#######################################################
#add annotation tag to each gene, that it's CDS overlap with an annotated TE
mkdir tmp

# use simple grep -f to lookup all genes that should be flagged and append a "TEoverlap=true" to the corresponding gff entry
cat $overlap|perl -pe 's/^(.*)$/ID=$1;/g' >  tmp/mRNA.repeat.overlap.tmp.lst
grep -f tmp/mRNA.repeat.overlap.tmp.lst $gene > tmp/tmp.overlapT.out
grep -v -f tmp/mRNA.repeat.overlap.tmp.lst $gene > tmp/tmp.overlapF.out
awk '{if ($3 == "mRNA") print $0}' tmp/tmp.overlapF.out > tmp/tmp.overlapF.mRNA.out
awk '{if ($3 != "mRNA") print $0}' tmp/tmp.overlapF.out > tmp/tmp.overlapF.other.out

# add true/false flag for repeatoverlap
cat tmp/tmp.overlapT.out|perl -pe 's/^(.*)$/$1;repeatoverlap=true/g' > tmp/tmp2.overlapT.out
cat tmp/tmp.overlapF.mRNA.out|perl -pe 's/^(.*)$/$1;repeatoverlap=false/g' > tmp/tmp2.overlapF.mRNA.out
# combine altered and unaltered gff entries
cat tmp/tmp.overlapF.other.out tmp/tmp2.overlapT.out tmp/tmp2.overlapF.mRNA.out > tmp/tmp.overlap.final.out

#######################################################


# REPBASE HIT
#######################################################
target=tmp/tmp.overlap.final.out
cut -f 1 $repbase|perl -pe 's/^(.*)$/ID=$1;/g' >  tmp/mRNA.repeat.repbase.tmp.lst
grep -f tmp/mRNA.repeat.repbase.tmp.lst $target > tmp/tmp.repbaseT.out
grep -v -f tmp/mRNA.repeat.repbase.tmp.lst $target > tmp/tmp.repbaseF.out
awk '{if ($3 == "mRNA") print $0}' tmp/tmp.repbaseF.out > tmp/tmp.repbaseF.mRNA.out
awk '{if ($3 != "mRNA") print $0}' tmp/tmp.repbaseF.out > tmp/tmp.repbaseF.other.out

# add true/false flag for repbasematch
cat tmp/tmp.repbaseT.out|perl -pe 's/^(.*)$/$1;repbasematch=true/g' > tmp/tmp2.repbaseT.out
cat tmp/tmp.repbaseF.mRNA.out|perl -pe 's/^(.*)$/$1;repbasematch=false/g' > tmp/tmp2.repbaseF.mRNA.out
# combine altered and unaltered gff entries
cat tmp/tmp.repbaseF.other.out tmp/tmp2.repbaseT.out tmp/tmp2.repbaseF.mRNA.out > tmp/tmp.repbase.final.out

#######################################################

# FUNAN FLAG
#######################################################
target=tmp/tmp.repbase.final.out
cut -f 1 $funan|perl -pe 's/^(.*)$/ID=$1;/g' >  tmp/mRNA.repeat.funan.tmp.lst
grep -f tmp/mRNA.repeat.funan.tmp.lst $target > tmp/tmp.funanT.out
grep -v -f tmp/mRNA.repeat.funan.tmp.lst $target > tmp/tmp.funanF.out
awk '{if ($3 == "mRNA") print $0}' tmp/tmp.funanF.out > tmp/tmp.funanF.mRNA.out
awk '{if ($3 != "mRNA") print $0}' tmp/tmp.funanF.out > tmp/tmp.funanF.other.out

# add true/false flag for funanmatch
cat tmp/tmp.funanT.out|perl -pe 's/^(.*)$/$1;functionalAnnotationTE=true/g' > tmp/tmp2.funanT.out
cat tmp/tmp.funanF.mRNA.out|perl -pe 's/^(.*)$/$1;functionalAnnotationTE=false/g' > tmp/tmp2.funanF.mRNA.out
# combine altered and unaltered gff entries
cat tmp/tmp.funanF.other.out tmp/tmp2.funanT.out tmp/tmp2.funanF.mRNA.out > tmp/tmp.funan.final.out

#######################################################

# transposonPSI hit

# final tmp gff for cleanup:
gff2clean=$base/tmp/tmp.funan.final.out

# cleanup gff3

agat_convert_sp_gxf2gxf.pl -g $gff2clean -o tmp/tmp.clean.gff3
gt gff3 -retainids -tidy -checkids -sort tmp/tmp.clean.gff3 > Cobs.alpha.v.2.1.geneannotation.1.3.TEflags.gff3
readlink -f Cobs.alpha.v.2.1.geneannotation.1.3.TEflags.gff3
```

# Filter gff based on different flags
see [alphaCobs.v2.1.TEproteinScreen.Rmd](alphaCobs.v2.1.TEproteinScreen.html) for processing GFF.
A file containing all genes to be retained in the final gene set is found at:
[./TEprotein.screen/genes2keep.tsv](./TEprotein.screen/genes2keep.tsv)

```bash
scp Cobs.alpha.v.2.1.geneannotation.1.3.TEflags2.gff3 jgant2:/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/filter/

# navigate to proper folder
base=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/filter
keeper=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/filter/genes2keep.tsv
gff2filter=$base/Cobs.alpha.v.2.1.geneannotation.1.3.TEflags2.gff3

cd $base

# filter out those genes that are meant to be kept
agat_sp_filter_feature_by_attribute_value.pl --gff $gff2filter --attribute keep --value false -v -p gene -o ./Cobs.alpha.v.2.1.geneannotation.1.4.keep.gff3

#sanity check
awk '{if ($3=="gene") print $0}' Cobs.alpha.v.2.1.geneannotation.1.4.keep.gff3 |perl -pe 's/.*?ID=(.*?);.*/$1/g' > Cobs.alpha.v.2.1.geneannotation.1.4.keep.lst
cat Cobs.alpha.v.2.1.geneannotation.1.4.keep.lst genes2keep.tsv |sort|uniq -c|egrep "1 "
```

# Finalize
```bash
mkdir /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results
cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results

cp /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/filter/Cobs.alpha.v.2.1.geneannotation.1.4.keep.gff3 Cobs.alpha.v.2.1.geneannotation.1.4.gff3
```
# Rename gene and mRNA ids
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results
agat_sp_manage_IDs.pl --gff Cobs.alpha.v.2.1.geneannotation.1.4.gff3 --prefix COBS -o tmp.gff3  -p gene --ensembl T
# agat renames outfile to tmp3.gff
awk '{if ($3=="mRNA") print $9}' tmp3.gff |cut -f 1,2 -d ";"|perl -pe 's/(.*);Parent=(.*)/$1\tParent=$2\tID=$2-mRNA/g' > tmp.tsv
# incremental counter of occurence of a Parent ID
awk '{ if (word == $2) { counter++ } else { counter = 1; word = $2 }; print $0 "\t"$3"-" counter }' tmp.tsv |cut -f 1,4|perl -pe 's/ID=//g' > conversionTable.mRNA.tsv
# search and replace from lookup table
time awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1'  "conversionTable.mRNA.tsv" tmp3.gff > tmp.renamed.gff3

# sanity check
awk '{if ($3=="mRNA") print $0}' tmp.renamed.renamed3.gff |perl -pe 's/.*ID=(.*?);Parent=(.*?);.*/$1\t$2/g'|perl -pe 's/-mRNA-.*?\t/\t/g'|awk '{if ($1!=$2) print $0}'
cat tmp.renamed.renamed3.gff |perl -pe 's/COBSG000000/COBS/g' > Cobs.alpha.v.2.1.geneannotation.1.5.gff3

#head -n 20 tmp3.gff > tmp.20.gff3
#time awk  'NR==FNR{a[$1]=$2;next} {for (i in a) { gsub(i,a[i],$9)}}1' "conversionTable.mRNA.tsv" tmp.20.gff3

```
# Rename removed genes to CobsTE
```bash
# better: rename genes in gff version 1.3 then filter?
# OR: Also rename genes in the gff containing only the removed genes

cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/removed
removed=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/filter/Cobs.alpha.v.2.1.geneannotation.1.4.removed.gff3
agat_sp_manage_IDs.pl --gff $removed --prefix COBS -o tmp.removed.gff  -p gene --ensembl T
cat tmp.removed.gff |perl -pe 's/COBSG000000/CobsTE/g' > tmp.removed2.gff
awk '{if ($3=="mRNA") print $9}' tmp.removed2.gff |cut -f 1,2 -d ";"|perl -pe 's/(.*);Parent=(.*)/$1\tParent=$2\tID=$2-mRNA/g' > tmp.removed.tsv
awk '{ if (word == $2) { counter++ } else { counter = 1; word = $2 }; print $0 "\t"$3"-" counter }' tmp.removed.tsv |cut -f 1,4|perl -pe 's/ID=//g' > conversionTable.removed.mRNA.tsv
# search and replace from lookup table
time awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1'  "conversionTable.removed.mRNA.tsv" tmp.removed2.gff > tmp.removed.renamed.gff3

# sanity check
awk '{if ($3=="mRNA") print $0}' tmp.removed.renamed.gff3 |perl -pe 's/.*ID=(.*?);Parent=(.*?);.*/$1\t$2/g'|perl -pe 's/-mRNA-.*?\t/\t/g'|awk '{if ($1!=$2) print $0}'
cat tmp.removed.renamed.gff3 > Cobs.alpha.v.2.1.geneannotation.1.5.TEs.gff3

```

# Extract fasta sequences
```bash
# get protein sequences
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa

gffread -g $genome -y tmp.COBS.aa -w Cobs.alpha.v.2.1.geneannotation.1.5.exons.fa -x Cobs.alpha.v.2.1.geneannotation.1.5.cds.fa Cobs.alpha.v.2.1.geneannotation.1.5.gff3
gffread -g $genome -y tmp.TEs.aa -w Cobs.alpha.v.2.1.geneannotation.1.5.TEs.exons.fa -x Cobs.alpha.v.2.1.geneannotation.1.5.TEs.cds.fa Cobs.alpha.v.2.1.geneannotation.1.5.TEs.gff3

cat tmp.COBS.aa|perl -pe 's/\.$/\*/g'> Cobs.alpha.v.2.1.geneannotation.1.5.pep.fa
cat tmp.TEs.aa|perl -pe 's/\.$/\*/g'> Cobs.alpha.v.2.1.geneannotation.1.5.TEs.pep.fa

tar cvzf Cobs.alpha.v.2.1.geneannotation.1.5.tar.gz Cobs.alpha.v.2.1.geneannotation.1.5.*


```

# create clean functional annotation files for v1.5
```bash
ln -s /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/interpro/Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv .
cd /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/functionalAnnotation.v1.5
cat ../tmp/conversionTable.mRNA.tsv ../removed/conversionTable.removed.mRNA.tsv |perl -pe 's/COBSG000000/COBS/g' > conversionTable.all.tsv

# search and replace from lookup table
time awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1'  "conversionTable.removed.mRNA.tsv" tmp.removed2.gff > tmp.removed.renamed.gff3

awk 'FNR==NR{a[$1]=$2;next}{print a[$1]"\t"$0}' conversionTable.all.tsv Cobs.alpha.v.2.1.geneannotation.1.3.interproFormat.fa.tsv |cut -f 1,3- > Cobs.alpha.v.2.1.geneannotation.1.5.interproFormat.all.tsv

egrep "^COBS" Cobs.alpha.v.2.1.geneannotation.1.5.interproFormat.all.tsv > Cobs.alpha.v.2.1.geneannotation.1.5.interproFormat.tsv
egrep "^CobsTE" Cobs.alpha.v.2.1.geneannotation.1.5.interproFormat.all.tsv > Cobs.alpha.v.2.1.geneannotation.1.5.TEs.interproFormat.tsv


```

Make webapollo friendly gffs
```bash

cat Cobs.alpha.v.2.1.geneannotation.1.5.gff3 |perl -pe 's/(.*Parent=.*?);.*/$1/g' > Cobs.alpha.v.2.1.geneannotation.1.5.webapolloFriendly.gff3
```

<!----

# update gff database
db=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/combine/Cobs.alpha.v.2.1.geneannotation.1.3.TEflags.gff3
rm $db.db
gffutils-cli create $db

# split gene list into chunks of 3000
mkdir splitFiles
mkdir records
cd splitFiles
split --lines=3000 ../genes2keep.tsv
cd $base

for file in $base/splitFiles/*
do
  set=$(echo $file|perl -pe 's/.*\/(.*?)/$1/g')
  echo $set
  ids=$(cat $file |tr "\n" ","|sed 's/.$//')
  gffutils-cli children $db.db $ids  > records/$set.records.gff
done



# rough gff file containing all genes that are not overlapping TEs and not having a repbase match
cat records/* > tmp/tmp.goodGenes.gff3

# cleanup gff3
gffread tmp/tmp.goodGenes.gff3 -o tmp/tmp.goodGenes.clean.gff3 -F -O
sed '1d' tmp/tmp.goodGenes.clean.gff3 > tmp/tmp.goodGenes.clean2.gff3
gt gff3 -retainids -sort tmp/tmp.goodGenes.clean2.gff3 > Cobs.alpha.v.2.1.geneannotation.1.4.gff3

agat_sp_webApollo_compliant.pl -g Cobs.alpha.v.2.1.geneannotation.1.4.gff3 > Cobs.alpha.v.2.1.geneannotation.1.4.wa.gff3
cat Cobs.alpha.v.2.1.geneannotation.1.4.gff3 |perl -pe 's/(.*Parent=.*?);.*/$1/g' > Cobs.alpha.v.2.1.geneannotation.1.4.webapolloFriendly.gff3


# gff only containing those genes that are repbasematch=false and repeatoverlap=false
/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/combine/Cobs.alpha.v.2.1.geneannotation.1.4.wa.gff3
/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/combine/Cobs.alpha.v.2.1.geneannotation.1.4.gff3


```
