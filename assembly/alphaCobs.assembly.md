
# Nanopore Assembly


## Assembly with CANU v1.9
Canu was run on porechopped minION data produced in two runs ```Alpha_29-10-2019``` and ```alphaDNA_07-08-2019```.

```bash
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/assembly
fastq=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/FiltLong/alphaDNA.porechopped.filtlong.fq.gz

/global/homes/jg/schradel/software/canu-1.9/Linux-amd64/bin/canu \
  -s canu.cfg \
  -p CobsAlpha  \
  -correctedErrorRate=0.085  \
  genomeSize=190m \
  -nanopore-raw $fastq\
  -d $base/
```


##### canu.cfg
  * Only setting RAM and not specifying cores to request (no ```-pe smp 10```) in ```canu.cfg```.
    ```bash
    useGrid=remote
    gridEngine=sge
    gridOptions = "-l h_vmem=30G"
    ```

### CANU contig assembly {# a1}
The contig assembly produced by CANU is found at
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/assembly/CobsAlpha.contigs.fasta```**



![](alphaCobs.assembly.assets/alphaCobs.assembly-3a003cc3.pdf)
*Figure 1: Nx plot for CANU contig assembly (produced with [QUAST](# QUAST))*

## Scaffolding using MinIonData with SSPACE
We use the porechopped MinIon reads to scaffold the assembly.
```bash
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding
cd $base
porechopped=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/FiltLong/alphaDNA.porechopped.filtlong.fq.gz
seqtk seq -a $porechopped > $base/alphaDNA.porechopped.fa

genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/assembly/CobsAlpha.contigs.fasta
perl /global/homes/jg/schradel/software/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl -c $genome -p $base/alphaDNA.porechopped.fa -b $base/CobsAlpha.scf1.fasta -t 10
```

### minION-scaffolded assembly {# a2}
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding/CobsAlpha.scf1.fasta/scaffolds.fasta```**


![](assets/alphaCobs.assembly-9da067d3.pdf)
*Figure 2: Nx plot for [minION-scaffolded CANU assembly](# a2) (produced with [QUAST](# QUAST))*

## Scaffolding using 454 data with SSPACE

#### Retrieve data
We have two 454 long-insert libraries for the alpha colony. We can use these in addition to the MinIon Data for scaffolding.
* **20 kb C. obs alpha 454** data:  https://www.ncbi.nlm.nih.gov/sra/SRX692533[accn]
  ```https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-5/SRR1565732/SRR1565732.1```
* **8 kb C. obs alpha 454** data: https://www.ncbi.nlm.nih.gov/sra/SRX692534[accn]
  ```https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos/sra-pub-run-5/SRR1565733/SRR1565733.1```

#### 1. Convert SRA format to split FASTQ file
We use the SRA toolkit to convert the SRA file to 3 fastq files.
* *_1.fastq = 1st mate
* *_2.fastq = 2nd mate
* *_3.fastq = unpaired

see
* https://www.biostars.org/p/222122/

```
cd /global/homes/jg/schradel/data/Cardiocondyla/alpha.data/454/
fastq-dump --split-files SRR1565732.1 ## not working properly
fastq-dump --split-3 SRR1565732.1
fastq-dump --split-3 SRR1565732.1
```
#### Prepare ```library.txt``` for ```SSPACE_Standard_v3```

Library file ```library.txt```

According to ```/Volumes/CobsData/Genome/MWG 20 kb 8 kp lpe lib/MWG 20kb 8kb pel/data/ReadMe.pdf``` the libraries have the following characteristics:
* **"20kb library": 30kb inserts with ~10% variation**
* **"8kb library": 8kb inserts with ~10% variation**
```
Lib20kb bwasw /global/homes/jg/schradel/data/Cardiocondyla/alpha.data/454/SRR1565732.1_1.fastq /global/homes/jg/schradel/data/Cardiocondyla/alpha.data/454/SRR1565732.1_2.fastq 30000 0.1 FF
Lib08kb bwasw /global/homes/jg/schradel/data/Cardiocondyla/alpha.data/454/SRR1565733.1_1.fastq /global/homes/jg/schradel/data/Cardiocondyla/alpha.data/454/SRR1565733.1_2.fastq 8000 0.1 FF
```

#### Scaffolding with SSPACE
We use ```SSPACE_Standard_v3.0.pl``` to further scaffold the MinIon-scaffolded assembly.
see:
* https://sites.google.com/site/104kananongstest/home/sta_command/sta_result/sspace
* https://www.biostars.org/p/17104/
```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/454.scaffolding
cd $base/
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding/CobsAlpha.scf1.fasta/scaffolds.fasta
perl /global/homes/jg/schradel/software/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l library.txt -s $genome -T 10 -b CobsAlpha.scf2 -p 1
```

### 454-scaffolded assembly {# a3}
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/454.scaffolding/CobsAlpha.scf2/CobsAlpha.scf2.final.scaffolds.fasta```**

![](assets/alphaCobs.assembly-cea274ad.pdf)
*Figure 3: Nx plot for [454-scaffolded CANU assembly](# a3) (produced with [QUAST](# QUAST))*

# Gap filling with LR_closer
see https://github.com/CAFS-bioinformatics/LR_Gapcloser

We used **LR_closer** to fill N-gaps in the scaffolds. This is achieved by including the corrected reads generated by Canu

```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/gapfilling/
cd $base/

export PATH=$PATH:/global/homes/jg/schradel/software/LR_Gapcloser/src/
# chmod u+x /global/homes/jg/schradel/software/LR_Gapcloser/src/*
# gunzip -c ../../assembly/Cobs.alpha/CobsAlpha.correctedReads.fasta.gz > reads.corrected.fasta
ln -s /global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding/alphaDNA.porechopped.fa .
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/454.scaffolding/CobsAlpha.scf2/CobsAlpha.scf2.final.scaffolds.fasta
bash LR_Gapcloser.sh -i $genome -l alphaDNA.porechopped.fa -s n
```

### gap-filled assembly {#4}
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/gapfilling/Output/iteration-3/gapclosed.fasta```**

## Scaffolding using ILLUMINA data with SSPACE

#### Retrieve data
We have one large short insert Illumina library for the alpha colony. We can use this in addition to the MinIon and 454 Data for scaffolding.

#### Prepare ```library.txt``` for ```SSPACE_Standard_v3```

Library file ```library.txt```

According to ```/Volumes/CobsData/Genome/MWG 20 kb 8 kp lpe lib/MWG 20kb 8kb pel/data/ReadMe.pdf``` the libraries have the following characteristics:
* **"Illumina library": 200bp inserts **

```
Lib200bp bwasw /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq 200 0.1 FR
```


#### Scaffolding with SSPACE
We use ```SSPACE_Standard_v3.0.pl``` to further scaffold the MinIon-scaffolded assembly.
see:
* https://sites.google.com/site/104kananongstest/home/sta_command/sta_result/sspace
* https://www.biostars.org/p/17104/
```bash
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/illumina.scaffolding
cd $base/
genome=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/gapfilling/Output/iteration-3/gapclosed.fasta
perl /global/homes/jg/schradel/software/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l library.txt -s $genome -T 10 -b CobsAlpha.scf3 -p 1
```

### illumina-scaffolded assembly {# a5}
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/illumina.scaffolding/CobsAlpha.scf3/CobsAlpha.scf3.final.scaffolds.fasta```**

![Illumina](assets/alphaCobs.assembly-89effcd0.pdf)
*Figure 4: Nx plot for [illumina-scaffolded CANU assembly](# a5) (produced with [QUAST](# QUAST))*


## Polishing with ntedit
see https://github.com/bcgsc/ntedit
```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/ntedit
cp /global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/illumina.scaffolding/CobsAlpha.scf3/CobsAlpha.scf3.final.scaffolds.fasta ./genome.fa
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq .
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq .

cd $base
```
### Run ntedit-make
```
/global/homes/jg/schradel/software/ntEdit/ntedit-make ntedit draft=genome.fa reads=SRR1564444 k=64 cutoff=2 t=10
```

### ntedit-polished assembly {# a6}
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/ntedit/genome_k64_edited.fa```**



## Polishing with pilon
### Pilon round 1 {# a7}
see https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/pilon/mapping.html
```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon
cd $base
cp /global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/ntedit/genome_k64_edited.fa ./genome.fa
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq .
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq .

bwa index genome.fa
bwa mem -t 6 genome.fa SRR1564444_1.fastq SRR1564444_2.fastq | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam

samtools index mapping.sorted.bam

# pilon
java -Xmx20G -jar /global/homes/jg/schradel/software/pilon/pilon-1.23.jar --frags mapping.sorted.bam --genome genome.fa

```
#### Pilon round 1
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon/pilon.fasta```**


### Pilon round 2 {# a8}
see https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/pilon/mapping.html
```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2
cd $base
cp /global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon/pilon.fasta ./genome.fa
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq .
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq .

bwa index genome.fa
bwa mem -t 6 genome.fa SRR1564444_1.fastq SRR1564444_2.fastq | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam

samtools index mapping.sorted.bam

# pilon
java -Xmx20G -jar /global/homes/jg/schradel/software/pilon/pilon-1.23.jar --frags mapping.sorted.bam --genome genome.fa

```
#### Pilon round 2
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2/pilon.fasta```**

### Pilon round 3 {# a9}

```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon3
cd $base
cp /global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2/pilon.fasta ./genome.fa
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq .
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq .

bwa index genome.fa
bwa mem -t 6 genome.fa SRR1564444_1.fastq SRR1564444_2.fastq | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam

samtools index mapping.sorted.bam

# pilon
java -Xmx20G -jar /global/homes/jg/schradel/software/pilon/pilon-1.23.jar --frags mapping.sorted.bam --genome genome.fa

```

#### Pilon round 3
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon3/pilon.fasta```**


### Pilon round 4 {# a10}

```
base=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon4
cd $base
cp /global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2/pilon.fasta ./genome.fa
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_1.fastq .
ln -s /global/homes/jg/merrbii/data/alpha_NGS/SRR1564444/SRR1564444_2.fastq .

bwa index genome.fa
bwa mem -t 6 genome.fa SRR1564444_1.fastq SRR1564444_2.fastq | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam

samtools index mapping.sorted.bam

# pilon
java -Xmx20G -jar /global/homes/jg/schradel/software/pilon/pilon-1.23.jar --frags mapping.sorted.bam --genome genome.fa

```
#### Pilon round 4
**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon4/pilon.fasta```**


## FINAL ASSEMBLY
```
cd /global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly
cat /global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon4/pilon.fasta|cut -f 1 -d "|" > Cobs.alpha.2.0.fa
```

**```/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa```**
![](assets/alphaCobs.assembly-2cbcbffa.pdf)


### Continue here:
**[/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/alphaCobs.processingAssembly.md](alphaCobs.processingAssembly.md)**



# Quality Assessment
see https://github.com/NBISweden/workshop-genome_assembly/wiki/Tools

|Abbr|Assembly|type|N50|scaffolds|largest|total|N|[BUSCO n:4415](# BUSCO)|file|
|-|-|-|-|-|-|-|-|-|-|
|a1|canu|contigs|4.62 MB|291|9.92 MB|192.650 MB|0|C:72.9%[S:72.7%,D:0.2%],F:19.7%,M:7.4%|[/global/homes/jg/schradel/data/assemblies/Cobs.alpha/assembly/CobsAlpha.contigs.fasta](# a1)|
|a2|SSPACE-long|scf|5.45 MB|192|13.08 MB|193.179 MB|529104|C:73.0%[S:72.8%,D:0.2%],F:19.6%,M:7.4%|[/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding/CobsAlpha.scf1.fasta/scaffolds.fasta](# a2)|
|a3|SSPACE-454|scf|5.72 MB|156|13.08 MB|193.375 MB|724288|C:72.9%[S:72.7%,D:0.2%],F:19.7%,M:7.4%|[/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/454.scaffolding/CobsAlpha.scf2/CobsAlpha.scf2.final.scaffolds.fasta](# a3)|
|a4   |gapfilled   |scf   | 5.72 MB |156   |13.08 MB   | 193.376 MB  | 190005 |  C:73.0%[S:72.8%,D:0.2%],F:19.7%,M:7.3% |[/global/homes/jg/schradel/data/assemblies/Cobs.alpha/gapfilling/Output/iteration-3/gapclosed.fasta](#4)|
|a5   |SSPACE-Illumina   |scf   |  5.99 MB |132   |  13.08 MB| 193.376 MB  | 190119  |  C:72.8%[S:72.6%,D:0.2%],F:19.8%,M:7.4% |   [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/illumina.scaffolding/CobsAlpha.scf3/CobsAlpha.scf3.final.scaffolds.fasta](#5)|
|a6   |ntEdit-polished   | scf  | 6.02 MB  | 132	 | 13.13 MB   |  194.146 MB | 190119  | C:97.3%[S:96.8%,D:0.5%],F:2.0%,M:0.7%  |  [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/ntedit/genome_k64_edited.fa](#6)|
|a7   |pilon1-polished   | scf  | 6.03 MB  | 132 | 13.15 MB   |194.482 MB  | 189436  |  C:98.3%[S:97.8%,D:0.5%],F:1.1%,M:0.6% |  [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon/pilon.fasta](#7)|
|a8   |pilon2-polished   | scf  | 6.03 MB | 132 |   13.15 MB| 194.468 MB | 188842 | C:98.3%[S:97.8%,D:0.5%],F:1.1%,M:0.6% |  [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2/pilon.fasta](#8)|
|a9   |pilon3-polished   | scf  |6.03 MB |132|13.15 MB|194.467 MB| 188493| C:98.3%[S:97.8%,D:0.5%],F:1.1%,M:0.6%|  [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon3/pilon.fasta](#9)|
|a10   |pilon3-polished   | scf  |6.03 MB |132|13.15 MB| 194.467 MB|188493 |  C:98.3%[S:97.8%,D:0.5%],F:1.1%,M:0.6%|  [/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon4/pilon.fasta](#10)|**



### Files
```bash
a1=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/assembly/CobsAlpha.contigs.fasta
a2=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/minION.scaffolding/CobsAlpha.scf1.fasta/scaffolds.fasta
a3=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/454.scaffolding/CobsAlpha.scf2/CobsAlpha.scf2.final.scaffolds.fasta
a4=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/gapfilling/Output/iteration-3/gapclosed.fasta
a5=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/scaffolding/illumina.scaffolding/CobsAlpha.scf3/CobsAlpha.scf3.final.scaffolds.fasta
a6=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/ntedit/genome_k64_edited.fa
a7=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon/pilon.fasta
a8=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon2/pilon.fasta
a9=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon3/pilon.fasta
a10=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/polishing/pilon4/pilon.fasta
```
### Busco
#### Prepare BUSCO
```
cd /global/homes/jg/schradel/data/assemblies/Cobs.alpha/QC/BUSCO
export PATH=$PATH:"/global/projects/programs/source/augustus-3.3.1/bin/"
# cp -r /global/projects/programs/source/augustus-3.3.1/config/ .
export AUGUSTUS_CONFIG_PATH=$(readlink -f ./config/)
BUSCO_CONFIG_FILE=/global/projects/programs/source/busco/config/config.ini
```

#### Run BUSCO {# BUSCO}
```
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a1 --out a1 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a2 --out a2 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a3 --out a3 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a4 --out a4 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a5 --out a5 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a6 --out a6 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a7 --out a7 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a8 --out a8 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a9 --out a9 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a10 --out a10 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

```


#### QUAST {# QUAST}
```
cd /global/homes/jg/schradel/data/assemblies/Cobs.alpha/QC/QUAST/
nice quast.py $a1 -o a1.quast.out --eukaryote
nice quast.py $a2 -o a2.quast.out --eukaryote
nice quast.py $a2 -o a2.s.quast.out --eukaryote -s
nice quast.py $a3 -o a3.quast.out --eukaryote -s
nice quast.py $a4 -o a4.quast.out --eukaryote -s
nice quast.py $a5 -o a5.quast.out --eukaryote -s
nice quast.py $a6 -o a6.quast.out --eukaryote -s
nice quast.py $a7 -o a7.quast.out --eukaryote -s
nice quast.py $a8 -o a8.quast.out --eukaryote -s
nice quast.py $a9 -o a9.quast.out --eukaryote -s
nice quast.py $a10 -o a10.quast.out --eukaryote -s

# Quast compare two
nice quast.py -t 3 -o a1.vs.a8 -R $a1 $a8 --eukaryote -s

```

<!--


# Quality Assessment
see https://github.com/NBISweden/workshop-genome_assembly/wiki/Tools

```bash
# n50
a1=/global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/assembly/Cobs.alpha/CobsAlpha.contigs.fasta



n50simple.pl $a1
n50simple.pl $a2
n50simple.pl $a3
n50simple.pl $a5
n50simple.pl $a6
n50simple.pl $a7
n50simple.pl $a8

# number of scf
grep ">" -c $a1
grep ">" -c $a2
grep ">" -c $a3
grep ">" -c $a4
grep ">" -c $a5
grep ">" -c $a6
grep ">" -c $a7
grep ">" -c $a8

# total length
grep -v ">" $a1 | tr -d '\n' | wc -c
grep -v ">" $a2| tr -d '\n' | wc -c
grep -v ">" $a3| tr -d '\n' | wc -c
grep -v ">" $a4| tr -d '\n' | wc -c
grep -v ">" $a5| tr -d '\n' | wc -c
grep -v ">" $a6| tr -d '\n' | wc -c
grep -v ">" $a7| tr -d '\n' | wc -c
grep -v ">" $a8| tr -d '\n' | wc -c

# longest scf
samtools faidx $a1
cat $a1.fai|sort -k2,2 -nr|head

samtools faidx $a2
cat $a2.fai|sort -k2,2 -nr|head

samtools faidx $a3
cat $a3.fai|sort -k2,2 -nr|head

samtools faidx $a5
cat $a5.fai|sort -k2,2 -nr|head

samtools faidx $a6
cat $a6.fai|sort -k2,2 -nr|head

samtools faidx $a7
cat $a7.fai|sort -k2,2 -nr|head

samtools faidx $a8
cat $a8.fai|sort -k2,2 -nr|head

# Ns
cat $a1|grep ">" -v|egrep "N|n" -o|wc -l
cat $a2|grep ">" -v|egrep "N|n" -o|wc -l
cat $a3|grep ">" -v|egrep "N|n" -o|wc -l
cat $a4|grep ">" -v|egrep "N|n" -o|wc -l
cat $a5|grep ">" -v|egrep "N|n" -o|wc -l
cat $a6|grep ">" -v|egrep "N|n" -o|wc -l
cat $a7|grep ">" -v|egrep "N|n" -o|wc -l
cat $a8|grep ">" -v|egrep "N|n" -o|wc -l


# Quast
base=/global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/assembly/QC/QUAST
cd $base
cd polished
quast.py -t 3 -o $base -R $a1 $a5
cd pilon
quast.py -t 3 -o $base -R $a1 $a6
```
|Assembly|type|N50|scaffolds|largest|total|N|BUSCO n:4415|
|-|-|-|-|-|-|-|-|
|canu|contigs|2.43 MB|972|9.8 MB|186.625 MB|0|C:58.3%[S:58.2%,D:0.1%],F:24.2%,M:17.5%
|SSPACE-long|MinIon scaffolding|4.04 MB|466|12.1 MB|187.693 MB|996302|
|SSPACE-Standard|454 scaffolding|4.04 MB|384|12.1 MB|188.305 MB|1608852|
|SSPACE-Standard|gapFilled|4.04 MB|384|12.1 MB|188.305 MB|530551|C:59.2%[S:59.1%,D:0.1%],F:23.9%,M:16.9%|
|SSPACE-Standard|polished|4.05 MB|384|12.15 MB|189.236 MB|530551|C:93.7%[S:93.3%,D:0.4%],F:4.6%,M:1.7%|
|PILON 1st | pilon round 1   |4.07 MB   |384|12.19 MB | 189.847 MB|528250   | C:97.6%[S:97.2%,D:0.4%],F:1.7%,M:0.7%  |
|PILON 2nd | pilon round 2   |4.07 MB   |384|12.19 MB | 189.808 MB|526516   | C:98.2%[S:97.7%,D:0.5%],F:1.3%,M:0.5%  |
|PILON 3rd | pilon round 3   |4.07 MB   |384|12.19 MB | 189.815 MB|524764   | C:98.2%[S:97.7%,D:0.5%],F:1.3%,M:0.5%  |

### Busco
#### Prepare BUSCO
```
cd /global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/assembly/QC
export PATH=$PATH:"/global/projects/programs/source/augustus-3.3.1/bin/"
# cp -r /global/projects/programs/source/augustus-3.3.1/config/ .
export AUGUSTUS_CONFIG_PATH="/global/homes/jg/schradel/data/minIonData/Alpha_29-10-2019/assembly/QC/config"
BUSCO_CONFIG_FILE=/global/projects/programs/source/busco/config/config.ini
```

#### Run BUSCO
```
python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a1 --out canu --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a5 --out polished --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a6 --out pilon1 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a7 --out pilon2 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

python /global/projects/programs/source/busco/scripts/run_BUSCO.py --in $a8 --out pilon3 --lineage /global/projects/programs/source/busco/hymenoptera_odb9/ --mode genome --cpu 1 -f

```
