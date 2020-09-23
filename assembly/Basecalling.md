---
title: "Cobs Genome repeat landscapes: TEislands"
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

# Guppy basecalling
**GUPPY v 3.3.0** was used for base-calling minION data generated for **Cobs alpha** in two flow cells. We used default settings for base-calling.

#### Download most recent release of guppy
https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.3.0_linux64.tar.gz

#### Check resources
https://community.nanoporetech.com/posts/guppy-on-high-computer-clu
https://medium.com/@kepler_00/nanopore-gpu-basecalling-using-guppy-on-ubuntu-18-04-and-nvidia-docker-v2-with-a-rtx-2080-d875945e5c8d
https://github.com/rrwick/August-2019-consensus-accuracy-update

### Flow cell and kit
For **HMW DNA** of the ***C. obscurior*** **alpha** colony, we used
* **R9.4.1 FLO-MIN106 (FAK34237)**
* **SQK-LSK109**



## Running GUPPY 3.3.0 on a batch system

Basecalling was run on two minIon libraries generated for Cobs.alpha
```bash
library=Alpha_29-10-2019/
library=alphaDNA_07-08-2019
```


#### Prepare files & folders
```bash
cd /global/homes/jg/schradel/data/minIonData/$library/
datafolder=/global/homes/jg/schradel/data/minIonData/$library/softlinks/
mkdir $datafolder

## library1
#for i in $(ls /global/homes/jg/schradel/data/minIonData/alphaDNA_07-08-2019/alphaDNA_07-08-2019/20190807_1546_MN29376_FAK34237_3959fdcc/fast5/*.fast5)

##library2
for i in $(ls /global/homes/jg/schradel/data/minIonData/$library/alphaDNA_07-08-2019/20190807_1546_MN29376_FAK34237_3959fdcc/fast5/*.fast5)
do
  #echo $i
  file=$(echo $i|perl -pe "s/.*\///g")
  mkdir $datafolder/$file/
  cd $datafolder/$file/
  ln -s $i .
done

```


```
cd /global/homes/jg/schradel/data/minIonData/
readlink -f softlinks/*> list.of.folders.txt
mkdir /global/homes/jg/schradel/data/minIonData/$library/fq/
```
####Batch script:
files saved as ```batchGuppy.sh```

```bash
#!/bin/bash
#$ -S /bin/bash # tell the SGE to use bash
#$ -cwd         #run in current working dir
#$ -pe smp 2    # how many CPUs?
#$ -l h_vmem=3G # how much RAM?
#$ -w e         # reject jobs with error
#$ -V           # export all environment variables
#$ -N guppy     # job name

echo "Task id is $SGE_TASK_ID"
processing=$(sed "${SGE_TASK_ID}q;d" ${file})

# actually run the command
guppy_basecaller \
--config /global/homes/jg/schradel/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg \
--compress_fastq \
-i ${processing} \
-s ${fqPath}/${SGE_TASK_ID} \
–cpu_threads_per_caller 4 \
–num_callers 1
```

#### Run qsub on all fast5 files
```bash
qsub -v file="list.of.folders.txt" -v fqPath="/global/homes/jg/schradel/data/minIonData/$library/fq/" -t 1-$(wc -l < list.of.folders.txt) batchGuppy.sh
# cleanup
mkdir guppy.out
```
### Final fq files
```
#library1
/global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/fq/*/*.fastq.gz

#library2
/global/homes/jg/schradel/data/minIonData/alphaDNA_07-08-2019/fq/*/*.fastq.gz
```



--------------------------------


<!--
## Running GUPPY 3.3.0 on individual fast5 files

#### Basecalling with GPU
```{bash}
guppy_basecaller \
--config /global/homes/jg/schradel/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg \
--compress_fastq \
-i /global/homes/jg/schradel/data/minIonData/test/ \
-s /global/homes/jg/schradel/data/minIonData/test/fq2/ \
--gpu_runners_per_device 1
```


#### options with config
```{bash}
guppy_basecaller \
--config /global/homes/jg/schradel/software/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg \
--compress_fastq \
-i /global/homes/jg/schradel/data/minIonData/test/ \
-s /global/homes/jg/schradel/data/minIonData/test/fq/ \
–cpu_threads_per_caller 2 \
–num_callers 1
```
#### options with flow cell and kit
```{bash}
guppy_basecaller \
--compress_fastq \
-i /global/homes/jg/schradel/data/minIonData/test/ \
-s /global/homes/jg/schradel/data/minIonData/test/fq/ \
–cpu_threads_per_caller 2 \
–num_callers 1
--flowcell FLO-MIN106 \
--kit SQK-LSK109
``` -->
