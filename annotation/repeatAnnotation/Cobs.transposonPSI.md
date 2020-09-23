## transposonPSI on pallas

### Run transposon psi to identify potential TE-derived or TE-contaminated genes.
```bash
#edited the TransposonPSI.pl  to use 20 threads in blastall (-a 20). Uses blastgpg though unless nuc database
base=/usr/local/home/lschrader/data/IBE/Cardiocondyla/TransposonPSI/

```

```bash
#############################
# run TransposonPSI on each genome in parallel
cd $base
# GEMOMA
data=/corefac/cse/lukas/IBE/Cardiocondyla/CardioMinION/annotation/Cobs.alpha.v.2.1.GeMoMa/Cobs.alpha.v.2.1.GeMoMa.6ref.proteins.fa
perl ~/software/TransposonPSI_08222010/transposonPSI.pl $data prot 1>gemoma.PSI.out 2> gemoma.blast.err
# BRAKER
data=/corefac/cse/lukas/IBE/Cardiocondyla/CardioMinION/annotation/Cobs.alpha.v.2.1.BRAKER/augustus.hints.aa
perl ~/software/TransposonPSI_08222010/transposonPSI.pl $data prot 1>BRAKERv.1.0.PSI.out 2> BRAKERv.1.0.TPSI.err
# Cobs cleaned up v2.1 gene annotation
data=/corefac/cse/lukas/IBE/Cardiocondyla/annotation/Cobs.alpha.v.2.1.geneannotation.1.2.fa
perl ~/software/TransposonPSI_08222010/transposonPSI.pl $data prot 1>Cobs.alpha.v.2.1.geneannotation.1.2.PSI.out 2> Cobs.alpha.v.2.1.geneannotation.1.2.TPSI.err

############################
