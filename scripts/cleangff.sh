export EVMutils=/usr/local/home/lschrader/software/EVidenceModeler-1.1.1/EvmUtils
export target_genome=/corefac/cse/lukas/IBE/Cardiocondyla/CardioMinION/results/Cobs.alpha.2.1.fa

# install tools
## GFFUtils (see https://github.com/fls-bioinformatics-core/GFFUtils and https://github.com/fls-bioinformatics-core/GFFUtils/issues/14)
## genometools (gt) http://genometools.org/
## gffread http://ccb.jhu.edu/software/stringtie/gff.shtml


# navigate
cd ~/data/IBE/Cardiocondyla/annotation

# Cleanup gff output files
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= Cobs.alpha.v.2.1.combined_6ref_braker_v1.2.gff 2>/dev/null
  gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries Cobs.alpha.v.2.1.combined_6ref_braker_v1.2_clean.gff| sed s/prediction/mRNA/g|sed s/GAF/GeMoMa/g|sed s/\+AnnotationEvidence//g> tmp.gff 2>/dev/null
  #gt gff3 -sort -tidy  -addids -fixregionboundaries -checkids yes -retainids yes tmp.gff > tmp2.gff 2>/dev/null
  awk '{if ( $3=="gene") print $0 }' tmp.gff > tmp.genes.gff 2>/dev/null
  gffread tmp.gff -F -o tmp3.gff --force-exons -E 2>/dev/null
  cat tmp.genes.gff tmp3.gff |gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -|egrep "^###" -v|sed s/geneID/Parent/g > tmp4.gff
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= tmp4.gff -o tmp5.gff 2>/dev/null
  gt gff3 -retainids -sort tmp5.gff > Cobs.alpha.v.2.1.geneannotation.1.2.gff3
  perl $EVMutils/gff3_gene_prediction_file_validator.pl Cobs.alpha.v.2.1.geneannotation.1.2.gff3
  perl $EVMutils/gff3_file_to_proteins.pl   Cobs.alpha.v.2.1.geneannotation.1.2.gff3 $target_genome |sed $'s/[^[:print:]\t]/ /g' > Cobs.alpha.v.2.1.geneannotation.1.2.fa 2>/dev/null


# Cleanup gff output files: webapollo friendly
  gffread tmp.gff -o tmp3.wa.gff --force-exons -E 2>/dev/null
  cat tmp.genes.gff tmp3.wa.gff |gt gff3 -sort -tidy  -checkids yes -retainids yes -fixregionboundaries -|egrep "^###" -v|sed s/geneID/Parent/g > tmp4.wa.gff
  GFFcleaner --clean-replace-attributes --add-missing-ids --add-exon-ids --report-duplicates --insert-missing= tmp4.wa.gff -o tmp5.wa.gff 2>/dev/null
  gt gff3 -retainids -sort tmp5.wa.gff > Cobs.alpha.v.2.1.geneannotation.1.2.webapolloFriendly.gff3
  mv tmp* tmp/
