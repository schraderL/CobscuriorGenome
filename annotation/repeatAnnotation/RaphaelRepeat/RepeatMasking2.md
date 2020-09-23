#TE annotation pipeline

Used versions:
>RepeatMasker version open-4.0.8 - Local, other versions problematic with version of TRF (analysis broke down near the end)
>Search Engine: NCBI/RMBLAST [ 2.9.0+ ] -
>RebBaseLib-Sources (used because newest versions available)
>Artefacts RELEASE 20190301;
>Dfam RELEASE Dfam_3.0;
>RepBase RELEASE 20181026;


#Library prep

##Adding ArtTEdb (see above)
Download all fasta files from http://db.cbi.pku.edu.cn/arte/download.html to
```bash
cd ~/Master/Library
#Combine these files into one single file.
cat *.fasta > ArtTEdb.fasta
```

* **.fasta includes all fasta files downloaded from http://db.cbi.pku.edu.cn/arte/download.html
-> Delete all downloaded fasta files after cat*

```bash
#Change the Format of the ArTEdb.fasta file so that RepeatMasker is able to read it
sed -i s@\\t@#@g ArtTEdb.fasta
sed -i s@\\tSINE@#SINE@g ArtTEdb.fasta
sed -i s@\\tRetro@#Retro@g ArtTEdb.fasta
sed -i s@\\tDNA@#DNA@g ArtTEdb.fasta
sed -i s@\\tRC@#RC@g ArtTEdb.fasta
sed -i s@\\tLTR@#LTR@g ArtTEdb.fasta
sed -i s@\\tUnknown@#Unknown@g ArtTEdb.fasta
```

For the additional RepeatScout/PASTEClassifier sequences:
-> copy the R-skript resulting out-files into one directory, then do:
```bash
for f in *.txt ; do sed 's/^\([^acgt]\)/>\1/' $f > $f.fa; done
for f in *.fa ; do tr "\t" "\n" < $f > $f.fasta ; done
```
## PASTEClassifier and RepeatScout-run to find and classify de novo TE sequences

Create [PASTEClassifier.cfg](#PASTEClassifier.cfg) file with the following content in your working directory (change project name to $1 and project_dir to the appropiate path):

##Run the following bash script (example species: Cobs)
```bash
#$ -S /bin/bash
#$ -l h_vmem=2G
#$ -pe smp 8
#$ -o /global/homes/jg/r_stef04/Master/TEpipeline/RS/Cobs/CobsOut
#$ -e /global/homes/jg/r_stef04/Master/TEpipeline/RS/Cobs/CobsErr
#$ -wd /global/homes/jg/r_stef04/Master/TEpipeline/RS/Cobs
#$ -V

#RepeatScout
export PATH=$PATH:"/global/projects/programs/source/censor-4.2.29/src/lcfilter/nseg/"
# add REPET_PATH to PATH
export REPET_PATH=/global/projects/programs/source/REPET_linux-x64-2.5
export PATH=$REPET_PATH/bin:$PATH
# ADD libraries to $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/global/projects/programs/source/libs/
# add REPET to PYTHONPATH
export PYTHONPATH=$REPET_PATH:$PYTHONPATH
# add RepeatScout
export PATH=$PATH:/global/projects/programs/source/RepeatScout-1
# Define MYSQL settings
export REPET_HOST=ebbsrv05
export REPET_USER=schradel
export REPET_PW=schr4d3l
export REPET_DB=REPET_schradel2
#Data folder
export datafolder="/global/homes/jg/r_stef04/REPET"

base=/global/homes/jg/r_stef04/Master/TEpipeline/RS/Cobs
cd $base
cat /global/homes/jg/r_stef04/Master/GenomesTE/NewGenomes/Cobs/Cobs.fa | cut -f 1 -d "|" > genome.fa
genome=$(readlink -f genome.fa)

build_lmer_table -sequence $genome -freq $genome.lmer.frequency
RepeatScout -sequence $genome -output $genome.repeats.fas  -freq $genome.lmer.frequency
filter-stage-1.prl $genome.repeats.fas > $genome.repeats.fas.filtered_1

#PasteClassifier
cd $base
ln -s $datafolder/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_ntSeq_cleaned_TE.fa .
ln -s $datafolder/RepBase20.05_REPET.embl/repbase20.05_aaSeq_cleaned_TE.fa .

cut -f 1 -d " " $genome.repeats.fas.filtered_1|sed 's/=/./g'  > $genome.RSlib.fa
PASTEClassifier.py -i $genome.RSlib.fa -C PASTEClassifier.cfg -p
```

Step 2: REPETclassification
 - Combine the resulting PASTEClassifier Classif and $genome.RSlib.fa-files into fasta file with RepeatMasker-compatible headers.
 - First: Transform the fasta out-file of PASTEClassifier into a data frame using the the following bash command, then continue in R:
```bash
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <  <$genome.RSlib.fa> > $1.out.fa
```

Run R-Skript REPETclassification.R over each respective $1.out.fa file. This skript contains the following steps:
 - For-loop to get the classification data (only) from the Classif-file.
 - Creation of new columns containing only actual classes, families and subfamilies (merging other data as well).
 - Changing all entries into the appropiate RepeatMasker-compatible classifications.
 - Read in the $1.out.fa.
 - Create new headers for the final Fasta-file.
 - Combine new headers with their respective sequences.

End of R-Skript. Delete first superfluous line with rownames from final outfile, $1.txt.

If necessary:
Adding header ">" to future fasta
```sed 's/^\([^acgt]\)/>\1/' $1.txt > $1.fa```
Exchanging "tab"-Delimiter with New Line characters in HymRep3.txt
```tr "\t" "\n" < $1.fa > $1.fasta```

Repeat for each genome. Concenate the resulting out-files into one big file.

## Filter out all non-arthropod seqs out of RepBase

```bash
cd /global/homes/jg/r_stef04/Master/RepeatMasker/RepeatMasker/util
perl queryRepeatDatabase.pl -species "arthropods" > lib.lib
```

with RepBase RELEASE 20181026 only *(in /global/homes/jg/r_stef04/Master/RepeatMasker/RepeatMasker/Libraries/TheRest/RepeatMaskerNoLin.embl, swap with normal RepeatMaskerLib.embl in Libraries and rename it to RepeatMaskerLib.embl since it doesn't contain the "linear;" identifier which has proven problematic because it replaces the actual identifier of the sequences)*
-> copy to ~/Master/Library
-> delete the command lines at the start from lib.lib
-> rename to RepBaseArthro.fasta

### Concatenating and cleaning library

- Concenate with custom TE/Repeat-library consisting of arthropod sequences of RebBase, a de novo sequence RepeatScout library consisting of sequences of all analysed species and the hymenoptera sequences of the ArTEdb data base

```bash
cd ~/Master/Library
cat *.fasta > CustomArTElib.fasta
```

Change towards our classification termini to avoid needless fragmentation of results
```bash
sed -i s@#RC@#DNA/RC@g CustomArTElib.fasta
sed -i s@#DNA/Helitron@#DNA/RC/Helitron@g CustomArTElib.fasta
```
Manual control of whether classification is correct for all sequences, especially in regards to the ArTEdb.

## Filtering Output file by removing protein-gene-similar sequences:

First: Extract and blast unknown sequences against annotated OR genes:
```bash
awk '/^>/ {P=index($0,"#Unknown")!=0} {if(P) print} ' CustomArTELib.fa > UnclassSeq

makeblastdb -in /home/r/r_stef04/Master/Library/$1/$1.OR.fa -dbtype prot

blastx -soft_masking false -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 7 -db Obir.OR.fa -query UnclassSeq -out ORObir
```
Second: Blast against protein annotation of genome
```bash
makeblastdb -in /home/r/r_stef04/Master/Library/Obir/Obir.OR.fa -dbtype prot

blastx -soft_masking false -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -outfmt 7 -db ObirProtein -query UnclassSeq -out ORObir
```

```bash
#Removing everything but name of reference seq from blast output
egrep "^#" -v ResultMarinerAllSimpleRep |cut -f 1|sort|uniq
```

Add all unique names of sequences with hits to a file called SeqNames, one line per sequence name and remove the classifier (everything after and including the '#') from each line.

```bash
sed -i s@#Unknown/Unknown@#@g SeqNames
sed -i s@#Unknown@#@g SeqNames
#Remove all protein similar unknown sequences from the fasta library
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' CustomArTELib.fa| grep -Fvf SeqNames - | tr "\t" "\n" > CustomArTELib.fa
```

- Remove redundancies in new Library with cd-hit-est
```bash
./cd-hit-est -T 0 -c 0.8 -n 5 -i CustomArTElibx.fasta -o CustomArTElib.fa
```

- Extract and remove all Unknown (i.e. not classified) Sequences from the composite library (rename it FinalCustomLibrary). Put the extracted sequences into an out.fa-file as input for the TEclass suit http://www.compgen.uni-muenster.de/tools/teclass/generate/index.pl?lang=en (09.0
3.2020). The file may have to be seperated into multiple smaller files before entering them into TEclass.
- Combine the resulting out-files into one single file, then change its format using the following steps:
```bash
cat UnknownClass* > UnknownCl
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <  /global/homes/jg/r_stef04/Master/L
ibrary/FinalLib/UnknownClass/UnknownCl > UC.fa
sed -i s@\|forward@@g UC.fa
sed -i s@\|reverse@@g UC.fa
sed -i s@\|ORF*@@g UC.fa
sed -i s@nonLTR@Retro@g UC.fa
sed -i s@DNAs@DNA' '@g UC.fa
sed -i s@Retros@Retro' '@g UC.fa
sed -i s@LINEs@LINE' '@g UC.fa
sed -i s@LTRs@LTR' '@g UC.fa
sed -i s@SINEs@SINE' '@g UC.fa
sed -i s@unclear@Unknown@g UC.fa
sed -i s@Unknowns@Unknown' '@g UC.fa
sed -i s@' 'complemented@@g UC.fa
```
- Delete everything in the headers of UC.fa after the classification:
 LINE[^\t]* -> LINE
 LTR[^\t]* -> LTR
 DNA[^\t]* -> DNA
 Retro[^\t]* -> Retro
 Unknown[^\t]* -> Unknown

- Extract header of sequences with the first part of the Rskript TEclassClass -> Output file is ListHeader.csv
- Seperate RepeatScout-Headers (Current format: [Species abbreviation]R.[Number]) from those of other sources, which are put into a seperate file (here AddUnCl), keeping the RepeatScout-Headers in ListHeader.csv.
- Use ListHeader.csv as input for the second part of TEclassClass -> output is ListHeader1.csv
- Extract Headers of the original sequence file (here out.fa) and sort AddUnCl according to it manually (TEclass resorts these sequences and they have to be in the previous order for the next step to work)
- Make two adjustments to the header lines before proceeding:
	- Delete the '|TEclass' string in all headers of AddUnCl
	- Replace 'Unknown/Unknown' with 'Unknown'
- Use the now sorted AddUnCl as input for the third part of TEclassClass -> Output is ListHeader2.csv
- Add the ListHeader2.csv headers according to the position of their sequences in out.fa to the RepeatScout headers. Since the latter are alphabetically sorted according to species name, there should be a single block between the sequences of two of those species in which the Unknown sequences of non-RepeatScout sequences can be found. The structure of the combined header out-file needs to be the same as that of the original out.fa-file. Name the resulting tab-delimited outfile SeqNames
- Run UnknownToHeader.py against out.fa, which will replace its header names with the ones from SeqNames
```bash
python UnknownToHeader.py -i out2.fa -r SeqNames -o FinalUnknownSeq
```
-> Add the FinalUnknownSeq-file to FinalCustomLibrary (position irrelevant).

## Complete Genome Annotation with concenated RepeatMaskerLib.fasta (RML.lib + REPET-Annotationen + ArTEdb)
```bash
cd ~/Master/RepeatMasker/RepeatMasker
perl RepeatMasker -s -gff -a -excln -cutoff 250 -nolow -lib /home/r/r_stef04/Master/Library/CustomArTElib.fasta -pa 1 $species
```

#### Further Filtering through Onecodetofindthemall:

```bash
mkdir $genomebase/$1/OneCode$1
cd /global/homes/jg/r_stef04/Master/RepeatMasker/Tutorial/Onecodetofindthemall
./build_dictionary.pl --rm $species.out --unknown > $genomebase/$1/OneCode$1/LTR$1 #unknown because certain sequences aren't classified, not fuzzy since it can lead to loss of information
#Add genome file and RM output file to $species/OneCode$1-folder
./one_code_to_find_them_all.pl --rm $species.out --ltr $genomebase/$1/OneCode$1/LTR$1 --strict --fasta --unknown $genomebase/$1/OneCode$1/$1.fa
```

Change .elem_sorted.csv output files into single gff3:
```bash
cd $species/OneCode$1
cat *.elem_sorted.csv|egrep "^###"|perl -pe 's/ +/\t/g'|cut -f 2-7,9-11,15,17 > $1TEsorted.csv

#if error 'bash: /bin/cat: Argument list too long' appears, do the following steps instead:

find ./ -name "*.elem_sorted.csv" -exec cat {} \; > /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted.csv
egrep "^###" /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted.csv > /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted2.csv
perl -pe 's/ +/\t/g' /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted1.csv > /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted2.csv
cut -f 2-7,9-11,15,17 /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted2.csv > /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted3.csv

#Rename the last output file to .../$1TEsorted.csv.

#Change every instance of LTR/Unknown or Unknown/ to Unknown (wrong class).

sed -i s@LTR/Unknown@Unknown@g /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted.csv
sed -i s@Unknown/@Unknown@g /global/homes/jg/r_stef04/Master/GenomesTE/$1/$1TEsorted.csv
```
Run OneCodeToGFF.R script (change the file paths in the script accordingly).
Replace first line of each gff-outfile with "##gff-version 3".

## bed tools analysis
Preparation:
```bash
#Create directories including all $1TEsorted.csv, all $1TEsorted.csv.bed, all $1TEsorted.csv.gff and all $1.fa genome files.
#cd to this directory
#Genome files to bed
for f in *.fa; do samtools faidx $f > /global/homes/jg/r_stef04/Master/GenomesTE/TEannotation/Fasta/$f_1.fa; done
for f in *.fai; do awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $f > /global/homes/jg/r_stef04/Master/GenomesTE/TEannotation/Fasta/$f.bed; done
```
Calculate total coverage of TEs, classified TEs and seperate subfamilies using the TEcontent.R script. To do that prepare multiple files.

```bash
cd /home/r/r_stef04/Master/ResultFiles/TEbed
#"Merging" the coverage of loci with multiple, nested TE sequences, so that bases are not included into the total sum multiple times.
for f in *.csv.bed ; do cut -f 1 -d ";" $f|bedtools merge -i - -c 4 -o distinct > $f.MergedTE.bed; done
#Remove unclassified sequences from bed files
for f in *.csv.bed ; do sed -n '/Unknown/!p' $f > $f.TEClass.bed;done
for f in *.TEClass.bed ; do cut -f 1 -d ";" $f|bedtools merge -i - -c 4 -o distinct > $f.MergedClassTE.bed; done
#Rename files to $110KB.bed
```

## Utilities:

```bash
#Genome quality
/home/r/r_stef04/Master/Software/assembly-stats-master/assembly-stats -t /home/r/r_stef04/Master/GenomesTE/$1/$1.fa

#Embl into FASTA
perl ~/Master/RepeatMasker/RepeatMasker/util/buildRMLibFromEMBL.pl <file>

#Reworking Summary Output
perl /global/homes/jg/r_stef04/Master/RepeatMasker/RepeatMasker/util/buildSummary.pl -species <> <input> > <output>

#Transform fasta into table file:
 ./seqkit fx2tab /home/r/r_stef04/Master/denovoHymenopteraRepeats.fa -o DeHymRep

#Transform table into fasta
./seqkit tab2fx /home/r/r_stef04/Master/HymRep.txt -o DeHymRep.fasta

#Potential Host Genes Blast against nt database
blastx -db /global/databases/ncbi-nr/2018_04_27/nr -query PotentialHostObirREPET.fa -out PotHosObiREP -outfmt 6
while read line<&3; do      proteinacc=$(echo $line | cut -d " " -f2);      ~/bin/edirect/elink -db protein -id "$proteinacc" -target gene | ~/bin/edirect/efetch -format uilist | sed '2!d'; done 3< ~/Master/GenomesTE/Obir/PHOR >> ~/Master/GenomesTE/Obir/PotHosObiREPTab

#For transformation of fasta into data.frame (\n -> \t) (R skript to transform REPET annotations into RepeatMasker-readable library needs this)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <  /global/homes/jg/r_stef04/Master/REPET-master/annotations/Obir/ObirDenovo/Obir_Blaster_GrpRecPil_Map_TEclassif_Filtered/Obir_sim_denovoLibTEs_filtered.fa > out.fa

#Extracting specific sequence via header from Fasta
#1. Linearize the fasta file:
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < /home/r/r_stef04/Master/GenomesTE/Obir/Obir.fa > Obir1.fa
tail -n +2 Obir1.fa > Obir22.fa
#grep the sequence
grep -A1 --no-group-separator ">NC_039506.1 Ooceraea biroi isolate clonal line C1 chromosome 1" Obir22.fa > ObirChr2.fa

#Extracting fasta sequences using list-file with parts of header names (example: multifasta file here: /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/Obir.OR.fa; listfile here: /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/ORGeneListArrays.txt)
IFS=$'\n'; for i in $(cat /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/ORGeneListArrays.txt);do line=$(grep -nr "$i" /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/Obir.OR.fa); if [[ ! -z $line ]];then for j in $line;do lineNumber=$(echo $j | cut -d':' -f1); sed -n "$lineNumber p" /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/Obir.OR.fa; awk -v nb=$lineNumber 'NR > nb {if ($0 ~ ">") exit; else print $0 }' /global/homes/jg/r_stef04/Master/ORgenes/Obir/final/finalSet/Obir.OR.fa; done;fi;done > ObirArray.fasta
```

------------------------------------

#PASTEClassifier.cfg {#PASTEClassifier.cfg}
```
[repet_env]
repet_version: 1.0
repet_host: ebbsrv05
repet_user: schradel
repet_pw: schr4d3l
repet_db: REPET_schradel2
repet_port: 3306
repet_job_manager: SGE

[project]
project_name: Cobs
project_dir: /global/homes/jg/r_stef04/Master/TEpipeline/RS/Cobs

[detect_features]
term_rep: yes
polyA: yes
tand_rep: yes
orf: yes
blast: blastplus
TE_BLRn: yes
TE_BLRtx: yes
TE_nucl_bank: repbase20.05_ntSeq_cleaned_TE.fa
TE_BLRx: yes
TE_prot_bank: repbase20.05_aaSeq_cleaned_TE.fa
HG_BLRn: no
HG_nucl_bank: <bank_of_host_genes>
TE_HMMER: no
TE_HMM_profiles: ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm
TE_HMMER_evalue: 10
rDNA_BLRn: no
rDNA_bank: <bank_of_rDNA_sequences_from_eukaryota>
tRNA_scan: no
TRFmaxPeriod: 15
clean: yes

[classif_consensus]
max_profiles_evalue: 1e-3
min_TE_profiles_coverage: 20
min_HG_profiles_coverage: 75
max_helitron_extremities_evalue: 1e-3
min_TE_bank_coverage: 5
min_HG_bank_coverage: 95
min_rDNA_bank_coverage: 95
min_HG_bank_identity: 90
min_rDNA_bank_identity: 90
min_SSR_coverage: 0.75
max_SSR_size: 100
remove_redundancy: no
min_redundancy_identity: 95
min_redundancy_coverage: 98
rev_complement: no
add_wicker_code: yes
add_noCat_bestHitClassif: yes
clean: yes
```
