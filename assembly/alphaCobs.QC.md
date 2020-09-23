# Nanopore raw read QC


## Porechop to trim adapters

QC and porechop was run on two minIon libraries generated for Cobs.alpha
```bash
# /global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019
library=Alpha_29-10-2019
# /global/homes/jg/schradel/data/minIonData/alphaDNA_07-08-2019
library=alphaDNA_07-08-2019
```

#### Raw data QC with minIONQC
```bash
# minIONQC scans for the sequencing_summary.txt files produced by guppy
cd /global/homes/jg/schradel/data/minIonData/$library/QC/minIONQC
 Rscript /global/homes/jg/schradel/bin/MinIONQC.R -p 6 -i /global/homes/jg/schradel/data/minIonData/$library/fq/ --outputdirectory /global/homes/jg/schradel/data/minIonData/$library/QC/minIONQC/raw/ --combined-only=TRUE

# final output
ls /global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/QC/minIONQC/raw/combinedQC
ls /global/homes/jg/schradel/data/minIonData/alphaDNA_07-08-2019/QC/minIONQC/raw/combinedQC
```

#### Prepare files & folders for porechop
We used ```porechop``` to remove adapter contaminations from minION reads.

```bash
mkdir /global/homes/jg/schradel/data/minIonData/$library/QC/porechop/
cd /global/homes/jg/schradel/data/minIonData/$library/QC/porechop/
mkdir ./porechopped/
datafolder=/global/homes/jg/schradel/data/minIonData/$library/fq/*/*.fastq.gz
readlink -f $datafolder> list.of.files.txt
```

Batch job script for porechop ```batchPorechop.sh```
```bash
#!/bin/bash
#$ -cwd         #run in current working dir
#$ -pe 1    # how many CPUs?
#$ -l h_vmem=2G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N porechop     # job name
#$ -S /bin/bash    # tell the SGE to use bash

# actual command
echo "Task id is $SGE_TASK_ID"
processing=$(sed "${SGE_TASK_ID}q;d" ${file})
porechop -i ${processing} -o ./porechopped/${SGE_TASK_ID}.porechop.fastq.gz
```

#### Submit as array job
```
qsub -v file="list.of.files.txt" -t 1-$(wc -l < list.of.files.txt) batchPorechop.sh
```
#### cleanup
```
mkdir porechopOut/
mv porechop.e* porechopOut
mv porechop.o* porechopOut

```
## Filtering reads with FiltLong
* see https://github.com/rrwick/Filtlong

```bash
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/FiltLong
base1=/global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/
base2=/global/homes/jg/schradel/data/minIonData/alphaDNA_07-08-2019/
```

1. Collect Illumina data for minION QC ***Cardiocondyla*** alpha **Illumina**
   ```bash
   illuminaData=/global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq
   seqtk sample $illuminaData 0.1 > $base/Cobs.SRR1564444_1.sample.fq
   ```
2. Collect porechopped fastq files from both libraries and concatenate
    ```bash
    cat $base1/QC/porechop/porechopped/*.fastq.gz $base2/QC/porechop/porechopped/*.fastq.gz > $base/alphaDNA.porechopped.fq.gz
    ```
3. Filter reads with FiltLong
    ```bash
    nanoporeReads=$base/alphaDNA.porechopped.fq.gz
    illSample=$base/Cobs.SRR1564444_1.sample.fq
    /global/homes/jg/schradel/software/Filtlong/bin/filtlong --split 500 --min_length 3000 --keep_percent 90 --target_bases 6000000000 --trim --illumina_1 $illSample $nanoporeReads | gzip > $base/alphaDNA.porechopped.filtlong.fq.gz
    ```

Final set of filtered and trimmed reads are found at
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/FiltLong/alphaDNA.porechopped.filtlong.fq.gz```**
