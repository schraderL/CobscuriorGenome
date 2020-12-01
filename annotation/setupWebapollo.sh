# Working progress!!
# on jgant2
cd /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha

# Cobs2.1.alpha
species="Cobs.a.2.1"
genome="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa"
geneGFF="/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/Cobs.alpha.v.2.1.geneannotation.1.5.webapolloFriendly.gff3"
TEgeneGFF="/global/homes/jg/schradel/data/Cardiocondyla/geneAnnotation/results/Cobs.alpha.v.2.1.geneannotation.1.5.TEs.webapolloFriendly.gff3"
orGFF="/global/homes/jg/schradel/data/Cardiocondyla/Cobs.alpha.v2.1.OR.gff3"
repeatsGFF="/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/Cobs.alpha.v2.1.TE.RS2.gff3"
repetGFF="/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/CobsAlpha2.REPET.gff3"
ltrGFF="/global/homes/jg/schradel/data/Cardiocondyla/repeatAnnotation/LTRharvest/Cobs.alpha.2.1.ltrs.clustered.filtered.gff3"

software=/home/s/schradel/software/JBrowse-1.15.0/bin

awk '{if ($2=="CobsAlpha2_REPET_TEs") print $0}' $repetGFF > repet.te.gff
awk '{if ($2=="CobsAlpha2_REPET_SSRs") print $0}' $repetGFF > repet.ssr.gff
awk '{if ($2=="CobsAlpha2_REPET_tblastx") print $0}' $repetGFF > repet.blx.gff
awk '{if ($2=="CobsAlpha2_REPET_blastx") print $0}' $repetGFF > repet.blx2.gff



perl $software/prepare-refseqs.pl --fasta $genome --out ./$species
perl $software/flatfile-to-json.pl --gff $geneGFF --type mRNA --trackLabel CobsOGS.1.5 --out ./$species &
perl $software/flatfile-to-json.pl --gff $TEgeneGFF --type mRNA --trackLabel CobsOGS.removed.1.5 --out ./$species &
perl $software/flatfile-to-json.pl --gff $orGFF --type mRNA --trackLabel OR --out ./$species &
perl $software/flatfile-to-json.pl --gff $repeatsGFF --type "transposable element" --trackLabel repeat --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.te.gff --trackLabel repetTE --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.ssr.gff --trackLabel repetSSR --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.blx.gff --trackLabel repetBLX --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.blx2.gff --trackLabel repetBLX2 --out ./$species &
perl $software/flatfile-to-json.pl --gff $ltrGFF --type "repeat_region" --trackLabel LTRharvest --out ./$species &
perl $software/generate-names.pl -v --out ./$species

# Add bigwig file
cut -f 1-2 $genome.fai > chrom.size

sort -k 1,1 -k2,2n /global/homes/jg/schradel/sciebo/Cobs2.1_annotation/GeMoMa_out/all_samples_coverage.bedgraph > bedSort.bedgraph
bedtools slop -i bedSort.bedgraph -g chrom.size -b 0 | bedClip stdin chrom.size bedSort.bedgraph.clip
bedGraphToBigWig bedSort.bedgraph.clip chrom.size coverage.bw

cd /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha/Cobs.a.2.1
mv ../coverage.bw /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha/Cobs.a.2.1
$software/add-bw-track.pl --in /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha/Cobs.a.2.1/trackList.json --bw_url coverage.bw --label RNAseqCoverage --key "RNAseq coverage"

# tar folder and copy
cd /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha
tar -zcvf $species.tar.gz $species 

scp $species.tar.gz schradel@iebapollo:/opt/apollo
ssh iebapollo

cd /opt/apollo
species="Cobs.a.2.1"
tar xvf $species.tar.gz
chmod -R 755 $species


# setup genome at https://iebapollo.uni-muenster.de:8443/apollo-2.1.0/annotator/index
#
