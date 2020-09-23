# Cardiocondyla alpha v.2.1 assembly

Documentation for the assembly, annotation and processing of C. obscurior genome project "alpha version 2". We sequenced the genome of the alpha colony of the C. obcurior lab population from Brazil. The genome is assembled from minIon data.

## Overview

1. Base calling raw nanopore data
[./assembly/Basecalling.md](./assembly/Basecalling.html)
2. Quality control of raw reads
[./assembly/alphaCobs.QC.md](./assembly/alphaCobs.QC.html)
3. Assembly with Canu and processing of assembly
[./assembly/alphaCobs.assembly.md](./assembly/alphaCobs.assembly.html)
4. Screen assembly for bacterial contaminations
[./assembly/alphaCobs.screenAssembly.md](./assembly/alphaCobs.screenAssembly.html)
[./assembly/alphaCobs.analyse.screenAssembly.Rmd](./assembly/alphaCobs.analyse.screenAssembly.html)
5. Repeat Annotation
   - REPET [./annotation/repeatAnnotation/REPET/CobsAlpha.2.10.REPET.sh](./annotation/repeatAnnotation/REPET/CobsAlpha.2.10.REPET.sh) *run on version 2.0*
   - Raphael [./annotation/repeatAnnotation/RaphaelRepeat/RepeatMasking2.md](./annotation/repeatAnnotation/RaphaelRepeat/RepeatMasking2.html) *(Raphael Steffen)*
     - Results: ```~/sciebo/Projects/CardiocondylaGenome/CardioMinION/annotation/Cobs.alpha.v.2.1.repeats/TE_newlibrary_19.03.20/Cobs.RepeatDataRaphael/Results```
   - LTRharvest [./annotation/repeatAnnotation/LTRharvest.md](./annotation/repeatAnnotation/LTRharvest.html)
   - transposonPSI [./annotation/repeatAnnotation/Cobs.transposonPSI.md](./annotation/repeatAnnotation/Cobs.transposonPSI.html)
6. Gene annotation (Jens Keilwagen & Katharina Hoff)

7. Process gene annotation (culminating in annotation v1.5)
[./annotation/geneAnnotation/alphaCobs.processGeneAnnotation.md](./annotation/geneAnnotation/alphaCobs.processGeneAnnotation.html)

8. Setup webapollo
[./annotation/setupWebapollo.sh](./annotation/setupWebapollo.sh)
