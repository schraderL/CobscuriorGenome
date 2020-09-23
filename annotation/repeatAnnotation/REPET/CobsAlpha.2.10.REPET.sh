#Cardiocondyla obscurior alpha strain
# Cobs2.10 assembly
species=CobsAlpha2
fasta=/global/homes/jg/schradel/data/assemblies/Cobs.alpha/finalAssembly/Cobs.alpha.2.0.fa

# prepare environment
source /global/homes/jg/schradel/projects/REPET/TEdenovoPipeline/sourceTEdenovo.sh $species $fasta
bash /global/homes/jg/schradel/projects/REPET/TEdenovoPipeline/prepareTEdenovo.sh

# /global/scratch/schradel/REPET/CobsAlpha1
cd $base
#Run TEdenovo
bash /global/homes/jg/schradel/projects/REPET/TEdenovoPipeline/TEdenovo.sh


# Run TEverify
source /global/homes/jg/schradel/sciebo/Projects/repet/TEannotVerifyPipeline/sourceTEannotVerify.sh $species
bash  /global/homes/jg/schradel/sciebo/Projects/repet/TEannotVerifyPipeline/prepareTEannotVerify.sh
cd $base
bash /global/homes/jg/schradel/projects/REPET/TEannotVerifyPipeline/TEannot1.sh


# Run TEannot
source /global/homes/jg/schradel/sciebo/Projects/repet/TEannotGenomePipeline/sourceTEannotGenome.sh $species
mkdir $base
cd $base
bash /global/homes/jg/schradel/sciebo/Projects/repet/TEannotGenomePipeline/prepareTEannotGenome.sh
bash /global/homes/jg/schradel/sciebo/Projects/repet/TEannotGenomePipeline/TEannot2.sh
