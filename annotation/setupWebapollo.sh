# Working progress!!
# on jgant2
cd /global/homes/jg/schradel/data/webapollo/Cobs2.1.alpha

# Cobs2.1.alpha
species="Cobs.a.2.1"
genome="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/results/Cobs.alpha.2.1.fa"
geneGFF="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.GeMoMa/Cobs.alpha.v.2.1.GeMoMa.6ref.gff"
orGFF="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v2.1.ORannotation/Cobs.alpha.v2.1.OR.gff3"
repeatsGFF="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/Cobs.alpha.v2.1.TE.gff3"
repetGFF="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v2.1.REPET/CobsAlpha2.REPET.gff3"
brakerGFF="/global/homes/jg/schradel/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.BRAKER/augustus.hints.gtf"
software=/home/s/schradel/software/JBrowse-1.15.0/bin

awk '{if ($2=="CobsAlpha2_REPET_TEs") print $0}' $repetGFF > repet.te.gff
awk '{if ($2=="CobsAlpha2_REPET_SSRs") print $0}' $repetGFF > repet.ssr.gff
awk '{if ($2=="CobsAlpha2_REPET_tblastx") print $0}' $repetGFF > repet.blx.gff
awk '{if ($2=="CobsAlpha2_REPET_blastx") print $0}' $repetGFF > repet.blx2.gff



perl $software/prepare-refseqs.pl --fasta $genome --out ./$species
perl $software/flatfile-to-json.pl --gff $geneGFF --type gene --trackLabel GeMoMa --out ./$species &
perl $software/flatfile-to-json.pl --gff $geneGFF --type prediction --trackLabel GeMoMa_prediction --out ./$species &
perl $software/flatfile-to-json.pl --gff $orGFF --type mRNA --trackLabel OR --out ./$species &
perl $software/flatfile-to-json.pl --gff $repeatsGFF --type "transposable element" --trackLabel repeat --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.te.gff --trackLabel repetTE --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.ssr.gff --trackLabel repetSSR --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.blx.gff --trackLabel repetBLX --out ./$species &
perl $software/flatfile-to-json.pl --trackType CanvasFeatures --gff repet.blx2.gff --trackLabel repetBLX2 --out ./$species &
perl $software/flatfile-to-json.pl --gff $brakerGFF --trackLabel BRAKER --out ./$species
perl $software/flatfile-to-json.pl --gff $brakerGFF --type transcript --trackLabel BRAKER_transcript --out ./$species



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
tar -zcvf $species.tar.gz $species &


#perl $software/flatfile-to-json.pl --gff ./Aech_data/Aech.OR.exons.bls.sorted.gff --type HSP --trackLabel exonblastn --out ./$species
#perl $software/flatfile-to-json.pl --gff $repeatGFF --type mRNA --trackLabel REPBASE --out ./$species


# perl $software/flatfile-to-json.pl --gff ./Ains_data/Ains.fullORannotation.gff3 --type mRNA --trackLabel ORautomatic --out ./$species
# perl $software/flatfile-to-json.pl --gff ./Ains_data/Ains.OR.gff3 --type mRNA --trackLabel ORfinal --out ./$species
# perl $software/flatfile-to-json.pl --gff ./Ains_data/Ains.vs.Aech.GeMoMa.final.gff --type mRNA --trackLabel ORgemoma --out ./$species
# perl $software/flatfile-to-json.pl --gff ./Ains_data/Annotations.gff3 --type mRNA --trackLabel ORmanual_Raphael --out ./$species
scp $species.tar.gz schradel@iebapollo:/opt/apollo
ssh iebapollo

cd /opt/apollo
species="Cobs.a.2.1"
tar xvf $species.tar.gz
chmod -R 755 $species


# setup genome at https://iebapollo.uni-muenster.de:8443/apollo-2.1.0/annotator/index
#
