


# Global approach
1. Run RepeatModeler2 on genome
2. Run EDTA on genome
3. Combine all with CD-hit
4. Run TEsorter on remaining Unknowns
5. Run PASTEC on remaining Unknowns
6. Run deepTE
8. Blast unknowns against arthropod Uniprot and remove those that hit a protein
9. Download RepBase arthropod
10. Run RepeatMasker with Ant Library + RepBase

------------------------------------------
# De novo repeat prediction
De novo repeat annotation for Cobs.v2.1.alpha assembly
## RepeatModeler2
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/RepeatModeler

PATH=$(getconf PATH)
#reconfigure RM
cd ~/software/RepeatModeler-2.0.1
perl ./configure

cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/RepeatModeler
~/software/RepeatModeler-2.0.1/BuildDatabase -name Cobs.alpha.2.1 Cobs.alpha.2.1.fa
nohup ~/software/RepeatModeler-2.0.1/RepeatModeler -database Cobs.alpha.2.1 -pa 20 -LTRStruct >& run.out &

# final RepMod Library
Cobs.alpha.2.1-families.fa
```

## EDTA
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/EDTA
ln -s /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa .
EDTA.pl --threads 20 --genome Cobs.alpha.2.1.fa --sensitive 1

# final EDTA library
Cobs.alpha.2.1.fa.mod.EDTA.TElib.fa
```

------------------------------------------

#Remove redundancies
## CD-HIT
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/CDHIT
 ln -s ../RepeatModeler/Cobs.alpha.2.1-families.fa  .
 ln -s ../EDTA/Cobs.alpha.2.1.fa.mod.EDTA.TElib.fa .
 cat Cobs.alpha.2.1-families.fa Cobs.alpha.2.1.fa.mod.EDTA.TElib.fa > Cobs.alpha.2.1.RepMod.EDTA.fa
 cd-hit-est -c 0.8 -n 5 -i Cobs.alpha.2.1.RepMod.EDTA.fa -o Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa

cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa|cut -f 1 -d "#" > Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa

```

------------------------------------------
# Classify repeats in de novo Library

## PASTEclassifier
link: https://urgi.versailles.inra.fr/Tools/PASTEClassifier
paper: https://pubmed.ncbi.nlm.nih.gov/24786468/
The template config file is contained below: [PASTEclassifier.template.cfg](#PASTEclassifier-cfg)

### Prepare run
#### Define environment

```bash
PATH=$(getconf PATH)
projectID=CobsApastec
# add REPET_PATH to PATH
export REPET_PATH=/global/projects/programs/source/REPET_linux-x64-2.5
export PATH=$REPET_PATH/bin:$PATH
# ADD libraries to $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/global/projects/programs/source/libs/
# add REPET to PYTHONPATH
export PYTHONPATH=$REPET_PATH:$PYTHONPATH
# add trf
export PATH=$PATH:/global/homes/jg/schradel/data/REPET/software/trf409
# add seqkit
export PATH=$PATH:/global/homes/jg/schradel/data/REPET/software/seqkit_0.9.0

# Define MYSQL settings
export REPET_HOST=ebbsrv05
export REPET_USER=schradel
export REPET_PW=schr4d3l
export REPET_DB=REPET_schradel2
#Data folder
export datafolder="/global/homes/jg/schradel/data/REPET/"

base=/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/PASTEC

```

#### Retrieve preformatted libraries and config
```bash
cd $base
ln -s $datafolder/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa .

ln -s $datafolder/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm.* .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa_* .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa_* .

cat $datafolder/PASTEclassifier.template.cfg|sed "s|PROJECTDIR|$base|g"|sed "s|PROJECTNAME|$projectID|g" > PASTEclassifier.cfg

```

#### Retrieve fasta file to classify
```bash
cd $base
#ln -s ../RepeatModeler/Cobs.alpha.2.1-families.fa .
ln -s ../CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa .

```

#### Clean MySQL database
Remove all mysql tables for a give projectID
```bash
mysql -u $REPET_USER -D $REPET_DB -e "show tables" -s --host $REPET_HOST -p'*****' | egrep "^$projectID" |xargs -I "@@" echo 'mysql -u $REPET_USER -D $REPET_DB -p"schr4d3l" -e "DROP TABLE @@" --host $REPET_HOST'
```

### Run PASTEClassifier
```bash
##get unclassified FASTA
#cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |seqkit grep -p "#Unknown" -nr |cut -f 1 -d " "|sed "s/#/./g"|sed "s/\//-/g"> $projectID.unknown.fa
#PASTEClassifier.py -i $projectID.unknown.fa -C PASTEclassifier.cfg -p

##get all de novo repeats
cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |cut -f 1 -d " "|sed "s/#/./g"|sed "s/\//-/g"> $projectID.fa
PASTEClassifier.py -i $projectID.fa -C PASTEclassifier.cfg -p


mkdir tmp
mv * tmp
mv tmp/*negStrandReversed* .
# tmp_negStrandReversed.classif contains the annotation from PASTEC
# It only recognizes MITES that were previously not annotated
## PasteClassifier output
```
The file ```*_negStrandReversed.classif``` contains the annotations in tsv.

#deepTE
Link: https://github.com/LiLabAtVT/DeepTE
Paper: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa519/5838183
DeepTE is aimed to classify transposons with unknown classification via Convolutional Neural Network.

check ```/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/```

```bash

conda create -n py36 python=3.6
conda activate py36
conda install tensorflow-gpu=1.14.0
conda install biopython
conda install keras=2.2.4
conda install numpy=1.16.0

cd ~/software/
git clone https://github.com/LiLabAtVT/DeepTE.git

mkdir /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/
# classify domains first
cd ~/software/DeepTE
python DeepTE_domain.py \
-d /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/ \
-o /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/ \
-i /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa \
-s /global/homes/jg/schradel/software/DeepTE/supfile_dir/ \
--hmmscan /global/projects/programs/bin/hmmscan

# classify TEs informed by the domain search
cd ~/software/DeepTE
python DeepTE.py \
-d /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/ \
-o /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/ \
-i /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa \
-m_dir /global/homes/jg/schradel/dbs/deepTE/download_M_model_dir/Metazoans_model \
-modify /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/opt_te_domain_pattern.txt \
-sp M



```

## Run deepTE to differentiate CDS and TEs
```bash
cd ~/software/DeepTE
python DeepTE.py \
-d /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/UNS \
-o /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/UNS \
-i /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/UNS/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa \
-m_dir /global/homes/jg/schradel/dbs/deepTE/download_U_model_dir/UNS_model \
-UNS yes

```

# TEsorter
link: https://github.com/zhangrengang/TEsorter
TEsorter can be used to classify any TE sequence, including Class I and Class II elements which are covered by the [REXdb](http://repeatexplorer.org/?page_id=918) database. It is installed along with EDTA.

```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEsorter
ln -s ../CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.no .

TEsorter Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa -p 20 -cov 10 -eval 1e-3
```

------------------------------------------
# Gather classifications
```bash
# deepTE
/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/opt_DeepTE.txt
# TEsorter
/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEsorter/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa.rexdb.cls.tsv
# PASTEC
/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/PASTEC/CobsApastec_negStrandReversed.classif
```

```bash
grep ">" ../TEsorter/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa.rexdb.cls.lib |tr " " "\t"> TEsorter.tsv
grep ">" ../PASTEC/CobsApastec_negStrandReversed_WickerH.fa |perl -pe 's/(.*?)_(.*)/$1\t$2/g' > PASTEC.tsv
grep ">" ../deepTE/opt_DeepTE.fasta |perl -pe "s/__/\t/g" > deepTE.tsv

cat ../CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |grep  "#Unknown" |perl -pe 's/\>(.*?)\#.*/$1#/g' > unknown.lst
cat ../CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |grep  "#Unknown" |perl -pe 's/\>(.*?)\#.*/$1\./g' > unknown.pastec.lst

grep -f unknown.lst TEsorter.tsv|grep -v Unknown
grep -f unknown.lst deepTE.tsv|grep -v Unknown
grep -f unknown.pastec.lst PASTEC.tsv|grep -v noCat
```
------------------------------------------
#Remove host proteins
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/tmp/
ln -s ../CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa .


blastx -soft_masking false -evalue 1e-10 -outfmt 8 -db /global/databases/uniprot_arthropoda/2018_04_27/uniprot-arthropoda+6656.fasta -query Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa -out Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.bls -num_threads 20

cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.bls|cut -f 1 |sort|uniq > Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.hits

#cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |seqkit grep -p "#Unknown" -nr > unknown.fa
#blastx -soft_masking false -evalue 1e-10 -outfmt 8 -db /global/databases/uniprot_arthropoda/2018_04_27/uniprot-arthropoda+6656.fasta -query unknown.fa -out Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.unknown.bls -num_threads 20

```

------------------------------------------
#Combine everything
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/RepeatMasker/input
# retrieve repbase arthropods
perl /global/homes/jg/schradel/conda/share/RepeatMasker/util/queryRepeatDatabase.pl -species "arthropod" > RepBase.arthropod.fa
perl /global/homes/jg/schradel/conda/share/RepeatMasker/util/queryRepeatDatabase.pl > RepBase.fa
grep ">" RepBase.fa |perl -pe 's/.*(\#.+?)\ .*/$1/g'|sort|uniq > RepBase.lst

cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/RepeatMasker/input
# retrieve de novo Library
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/CDHIT/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa .

# retrieve host gene list
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/tmp/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.hits .

cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.fa |seqkit grep -f Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.hits -nrv >  Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.fa

# retrieve classifications
# deepTE
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/opt_DeepTE.txt .
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/deepTE/opt_DeepTE.fasta .
# TEsorter
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/TEsorter/Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa.rexdb.cls.lib .
# PASTEC
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/PASTEC/CobsApastec_negStrandReversed.classif .

# split fasta in known and unknown
cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.fa |seqkit grep -p "#Unknown" -nr > Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.fa
cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.fa |seqkit grep -p "#Unknown" -nrv> Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.known.fa

# create list of unknown repeats
grep ">" Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.fa|perl -pe 's/\>(.+?\#).*/$1/g' > Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.lst

# get classifications from TEsorter for these
#Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.lst
cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa.rexdb.cls.lib |seqkit grep -f Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.lst -nr|seqkit grep -p "#Unknown" -nrv |cut -f 1 -d " ">  TEsorter.classified.fa

cat Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noClassification.fa.rexdb.cls.lib |seqkit grep -f Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.unknown.lst -nr|seqkit grep -p "#Unknown" -nr |grep ">"|perl -pe 's/\>(.+?\#).*/$1/g' >  TEsorter.unclassified.lst

# get classifications from PASTEclassifier for the repeats unclassified by TEsorter
ln -s /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/PASTEC/CobsApastec_negStrandReversed_WickerH_noCatBestHit.fa .

cat CobsApastec_negStrandReversed_WickerH_noCatBestHit.fa |sed 's/\./\#/g'|perl -pe 's/(?:(#.*))\-/$1\//g'|seqkit grep -f TEsorter.unclassified.lst -nr|seqkit grep -p "noCat" -nrv > PASTEC.classified.fa
cat PASTEC.classified.fa|perl -pe 's/>(.+?)\_(.+)\#.*/\>$2\t$1/g'|sed 's/-incomp//g'|awk -F $'\t' 'BEGIN {OFS = FS} FNR==NR { a[$1]=$2; next } $2 in a { $2=a[$2] }1' PASTE2RepBase.tsv - |tr "\t" "#"> PASTEC.repBaseClassified.fa

cat CobsApastec_negStrandReversed_WickerH_noCatBestHit.fa |sed 's/\./\#/g'|perl -pe 's/(?:(#.*))\-/$1\//g'|seqkit grep -f TEsorter.unclassified.lst -nr|seqkit grep -p "noCat" -nr|grep ">"|perl -pe 's/\>.+?_(.*)/$1/g' > PASTEC.unclassified.lst

# get classifications from deepTE for the repeats unclassified by TEsorter & PASTEC

cat opt_DeepTE.fasta |seqkit grep -f PASTEC.unclassified.lst -nr |seqkit grep -p "Class" -nr > deepTE.classified.fa
cat opt_DeepTE.fasta |seqkit grep -f PASTEC.unclassified.lst -nr |seqkit grep -p "Class" -nrv > deepTE.unclassified.fa

cat deepTE.classified.fa|perl -pe 's/(.+)\#.+\_\_(.*)/$1\t$2/g' |awk -F $'\t' 'BEGIN {OFS = FS} FNR==NR { a[$1]=$2; next } $2 in a { $2=a[$2] }1' deepTE2RepBase.tsv - |tr "\t" "#" > deepTE.repBaseClassified.fa


# combine all
cat \
Cobs.alpha.2.1.RepMod.EDTA.nonRedundant.noHostGene.known.fa \
TEsorter.classified.fa \
PASTEC.repBaseClassified.fa \
deepTE.repBaseClassified.fa \
> Cobs.alpha.2.1.deNovoLibrary.fa

######################################
Cobs.alpha.2.1.deNovoLibrary.fa
######################################

perl ~/Master/RepeatMasker/RepeatMasker/util/buildRMLibFromEMBL.pl

cat RepBase.arthropod.fa Cobs.alpha.2.1.deNovoLibrary.fa > library.fa

```
# RepeatMasker
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/RepeatMasker

ln -s /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa .
RepeatMasker -s -gff -a -inv -small -excln -cutoff 250 -nolow -lib input/library.fa -pa 20 -html Cobs.alpha.2.1.fa

perl /global/homes/jg/schradel/conda/share/RepeatMasker/util/queryRepeatDatabase.pl
perl /global/homes/jg/schradel/conda/share/RepeatMasker/util/buildSummary.pl -species Cobs Cobs.alpha.2.1.fa.out > Cobs.alpha.2.1.summary
```
# ltr_finder
```bash
ltr_finder Cobs.alpha.2.1.fa
```


# Supplement

##PASTEclassifier config template {#PASTEclassifier-cfg}

```PASTEclassifier.template.cfg``` stored at ```/global/homes/jg/schradel/data/REPET/```

```
[repet_env]
repet_version: 2.5
repet_host: ebbsrv05
repet_user: schradel
repet_pw: schr4d3l
repet_db: REPET_schradel2
repet_port: 3306
repet_job_manager: SGE

[project]
project_name: PROJECTNAME
project_dir: PROJECTDIR

[detect_features]
resources: h_vmem=5G
tmpDir: /global/scratch/schradel/tmp
term_rep: yes
polyA: yes
tand_rep: yes
orf: yes
blast: blastplus
TE_BLRn: no
TE_BLRtx: yes
TE_nucl_bank: repbase20.05_ntSeq_cleaned_TE.fa
TE_BLRx: yes
TE_prot_bank: repbase20.05_aaSeq_cleaned_TE.fa
TE_HMMER: yes
TE_HMM_profiles: ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm
TE_HMMER_evalue: 10
HG_BLRn: no
HG_nucl_bank: <bank_of_host_genes>
rDNA_BLRn: no
rDNA_bank: <bank_of_rDNA_sequences_from_eukaryota>
tRNA_scan: no
TRFmaxPeriod: 15
RepScout: no
RepScout_bank:
clean: no

[classif_consensus]
resources: h_vmem=10G
tmpDir: /global/scratch/schradel/tmp
limit_job_nb: 0
max_profiles_evalue: 1e-3
min_TE_profiles_coverage: 20
min_HG_profiles_coverage: 75
max_helitron_extremities_evalue: 1e-3
min_TE_bank_coverage: 5
min_HG_bank_coverage: 95
min_HG_bank_identity: 90
min_rDNA_bank_coverage: 95
min_rDNA_bank_identity: 90
min_SSR_coverage: 0.75
max_SSR_size: 100
remove_redundancy: no
min_redundancy_identity: 95
min_redundancy_coverage: 98
rev_complement: yes
add_wicker_code: yes
add_noCat_bestHitClassif: yes
clean: no
```
<!---
#TE_class
## installation
TEclass installation runs a download script that downloads dependencies. Some of these are outdated and need to be manually corrected.
```bash
curl -o 'librf.tar.gz' 'https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/librf/librf.0.1.tar.gz'
curl -o 'blast.tar.gz' 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz'

# manually install libsvm with make
# embl REPBASE setup can be skipped
# download the classifiers
PATH=$PATH:"/global/homes/jg/schradel/software/TEclass-2.1.3/libsvm-3.24"
./TEclassTest.pl -o test -c ~/software/TEclass-2.1.3 testfile.fa
# can't install lvq_pak-3.1
```


https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Helitron_Scanner.html#gsc.tab=0
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/tmp
ln -s /global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa .
java -jar /global/homes/jg/schradel/software/HelitronScanner/HelitronScanner.jar scanTail -bs 1000000 -g Cobs.alpha.2.1.fa  -th 16 -o tail.helitronscanner.out
