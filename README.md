# Cardiocondyla alpha v.2.1 assembly and analysis

Documentation for the assembly, annotation and processing of C. obscurior genome project "alpha version 2". We sequenced the genome of the alpha colony of the C. obcurior lab population from Brazil. The genome is assembled from minIon data.

## Overview

1. Base calling raw nanopore data:
[./assembly/Basecalling.md](./assembly/Basecalling.md)
2. Quality control of raw reads:
[./assembly/alphaCobs.QC.md](./assembly/alphaCobs.QC.md)
3. Assembly with Canu and processing of assembly:
[./assembly/alphaCobs.assembly.md](./assembly/alphaCobs.assembly.md)
4. Screen assembly for bacterial contaminations:
[./assembly/alphaCobs.screenAssembly.md](./assembly/alphaCobs.screenAssembly.md)
[./assembly/alphaCobs.analyse.screenAssembly.Rmd](./assembly/alphaCobs.analyse.screenAssembly.Rmd)
5. Repeat Annotation:
   - REPET [./annotation/repeatAnnotation/REPET/CobsAlpha.2.10.REPET.sh](./annotation/repeatAnnotation/REPET/CobsAlpha.2.10.REPET.sh) *run on version 2.0*
   - RSannotation [./annotation/repeatAnnotation/RSannotation/RepeatMasking2.md](./annotation/repeatAnnotation/RSannotation/RepeatMasking2.md) *(Raphael Steffen)*
     - Results: ```~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results```
   - LTRharvest [./annotation/repeatAnnotation/LTRharvest.md](./annotation/repeatAnnotation/LTRharvest.md)
   - transposonPSI [./annotation/repeatAnnotation/Cobs.transposonPSI.md](./annotation/repeatAnnotation/Cobs.transposonPSI.md)
   - EDTA+RepMod [./annotation/repeatAnnotation/EDTA.RepMod.TEannotation.md](./annotation/repeatAnnotation/EDTA.RepMod.TEannotation.md)


6. Process gene annotation (culminating in annotation v1.5):
[./annotation/geneAnnotation/alphaCobs.processGeneAnnotation.md](./annotation/geneAnnotation/alphaCobs.processGeneAnnotation.md)

7. Setup webapollo:
[./annotation/setupWebapollo.sh](./annotation/setupWebapollo.sh)

8. Different analyses of genome structure, etc:
[./analysis/](./analysis/)
