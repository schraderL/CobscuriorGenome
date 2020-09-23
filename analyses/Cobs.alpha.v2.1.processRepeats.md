# Processing RepeatMasker output for Cobs.alpha.v2.1 assembly

```
# bioPerl installed with conda on jgpogo
```

```bash
#upload data
scp /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael.zip jgpogo:~/data/Cardiocondyla/repeatAnnotation/Raphael/
ssh jgpogo
base=~/data/Cardiocondyla/repeatAnnotation/Raphael
scripts=/global/homes/jg/schradel/software/Parsing-RepeatMasker-Outputs
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
```


## Create tsv file for Kimura divergence and landscape plits
```bash
cd $base
perl /global/projects/programs/source/RepeatMasker/util/calcDivergenceFromAlign.pl Cobs.fa.align -a Cobs.fa.align.new -s Cobs.fa.align.summary
```


## Create file for landscape plits and time estimates
**https://github.com/4ureliek/Parsing-RepeatMasker-Outputs**

**From help text of parseRM.pl**
<font size="1">
>There are 3 non exclusive parsings types (they can be set together):
 -p To get a summary of the masking, as well as amount or DNA,
      counts of fragments, + several other details for each repeat name
      (all-repeats file), family, class and total amount (summary file)
      To deal well with the .align positions instead of segments are considered,
      so it is slow (several hours on a bird genome)
 -a To determine the amounts of DNA in a genome that is masked by repeats
      of different lineages / %divergence categories
 -l To split the amount of DNA by bins of %div or My, allowing to generate
      landscape graphs for each repeat name, family or class (one output for each)

> -a,--age (STRING)
         Option1: load a file
            File should contain TE age data, tabular delimited. Columns should be:
            Rname \t Rclass \t Rfam \t Lineage
            The Rname needs to be an exact match to the repeat names in the RM output;
            class and family can be changed
            If -d is chosen, the file must contain info for ALL repeats of all files
         Option2: 1 or 2 values (INT or FLOAT). Use a comma to sepearate 2 values;
            their order does not matter (-a 13.5,15 and -a 15,13.5 are equivalent)
            Any DNA masked with %div < the smallest number given (here 13.5) will be
            placed as lineage specific. Any DNA masked with %div > the largest number
            given (here 15) will be placed as ancient, and the rest will be placed as "nd".
          **If -m is set, then numbers here HAVE TO BE in My instead of %divergence**
            If -d is chosen and different values should be used for each input files
               you can provide a file containing 2 columns: filename \t X
               filename = name of the file to be parsed
               X = the value (in nt)

</font>

## Run parseRM.pl
```bash
cd $base/
# m according to  mutation rate estimate recently published for bumblebees (3.6e-9 per nt per MY, Liu et al. 2017)
$scripts/parseRM.pl -i Cobs/Results/Cobs.fa.align -l 50,0.2 --parse --glen $genome -m 0.0036 --age 0.1,1

$scripts/parseRM.pl -i Cobs/Results/Cobs.fa.align -l 50,0.5 --parse --glen $genome -m 0.0036 --age 0.1,1


```

## Output
The two scripts above create the following output:
```
Cobs.fa.align.landscape.Div.Rclass.tab
Cobs.fa.align.landscape.Div.Rfam.tab
Cobs.fa.align.landscape.Div.Rname.landscape.gypsy.pdf
Cobs.fa.align.landscape.Div.Rname.tab
Cobs.fa.align.landscape.My.Rclass.tab
Cobs.fa.align.landscape.My.Rfam.tab
Cobs.fa.align.landscape.My.Rname.landscape.LTR.0-50.pdf
Cobs.fa.align.landscape.My.Rname.landscape.LTR.0.10.pdf
Cobs.fa.align.landscape.My.Rname.tab
Cobs.fa.align.new
Cobs.fa.align.parseRM.all-repeats.tab
Cobs.fa.align.parseRM.summary.tab
Cobs.fa.align.splitage.tab
Cobs.fa.align.summary
```

## Split Cobs.fa.align.new in TE-island and LDR
```bash
cd ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/landscapes/
perl ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/scripts/splitAlign.define.pl Cobs.fa.align.new > Cobs.fa.align.new.bed
ln -s /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v.2.1.TEislands/te_islands.bed .
bedtools intersect -a te_islands.bed -b Cobs.fa.align.new.bed -wb > Cobs.fa.align.new.TEisland.bed
perl ~/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/scripts/splitAlign.split.pl Cobs.fa.align.new Cobs.fa.align.new.TEisland.bed
mv nonoverlap.align Cobs.fa.LDR.align
mv overlap.align Cobs.fa.TEisl.align

# calculate TEisland total size
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' te_islands.bed
##24664824
# upload data
gzip Cobs.fa.TEisl.align
gzip Cobs.fa.LDR.align
scp Cobs.fa.*.align.gz jgpogo:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/landscape


#jgpogo
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/landscape
scripts=/global/homes/jg/schradel/software/Parsing-RepeatMasker-Outputs
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
seqkit stat $genome
# length=
# genome minus TE islands
# 193051228-24664824=168386404

gunzip Cobs.fa.*.align.gz

$scripts/parseRM.pl -i Cobs.fa.TEisl.align -l 50,0.2 --parse --glen 24664824 -m 0.0036 --age 0.1,1
$scripts/parseRM.pl -i Cobs.fa.LDR.align -l 50,0.2 --parse --glen 168386404 -m 0.0036 --age 0.1,1

# download data
scp jgant2:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/landscape/Cobs.fa.TEisl.align.* /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/landscapes/TEislands/

scp jgant2:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Raphael/Cobs/Results/landscape/Cobs.fa.LDR.align.* /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/landscapes/LDRs/

```

# Plotting
Continue in R with
```bash
/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.Rmd
/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.LDRs.Rmd
/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.TEislands.Rmd
```



# Process DNApipeTE output
```bash
cd /global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/DNApipeTE/
genome=/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/data/Cobs.alpha.2.1.fa
scripts=/global/homes/jg/schradel/software/Parsing-RepeatMasker-Outputs
ln -s Itabuna.Trinity.fasta.cat Itabuna.Trinity.fasta.align
ln -s Tenerife.Trinity.fasta.cat Tenerife.Trinity.fasta.align
$scripts/parseRM.pl -i Itabuna.Trinity.fasta.align -l 50,0.2 --parse --glen 193051228 -m 0.0036 --age 0.1,1
$scripts/parseRM.pl -i Tenerife.Trinity.fasta.align -l 50,0.2 --parse --glen 193051228 -m 0.0036 --age 0.1,1

$scripts/parseRM.pl -i Itabuna.Trinity.fasta.out -l 50,0.2

perl /global/projects/programs/source/RepeatMasker/util/calcDivergenceFromAlign.pl Itabuna.Trinity.fasta.align -a Itabuna.Trinity.fasta.align.new -s Itabuna.Trinity.fasta.align.summary
perl /global/projects/programs/source/RepeatMasker/util/calcDivergenceFromAlign.pl Tenerife.Trinity.fasta.align -a Tenerife.Trinity.fasta.align.new -s Tenerife.Trinity.fasta.align.summary


```

## Plotting
Continue in R with
```bash
# download data
cd /Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/analyses/DNApipeTE/

scp jgpogo:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/DNApipeTE/*.summary .
scp jgpogo:/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/DNApipeTE/*.tab .

/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.Rmd
/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.LDRs.Rmd
/Users/lukas/sciebo/Projects/CardiocondylaGenome/CardioMinION/code/analyses/RepeatLandscapes/plotLandscapes.TEislands.Rmd
```
