# FRONTIER

### meth/preprocess.R
Data was processed using the `minfi` packages in R (R Foundation for Statistical Computing, Vienna). Data from the 450k (`IlluminaHumanMethylation450k.ilmn12.hg19`) and EPIC platforms (`IlluminaHumanMethylationEPICanno.ilm10b2.hg19`) were processed seperately. Detection P-values were calculated for each probe and sample, and samples with an average detection P-value > 0.01 were removed from follow-up analysis. Data was normalized using XXXX. Probes on sex chromosomes and known cross-reactive probes were removed, as were probes mapping to known SNPs and probes with a detection P-value > 0.01. Finally, data from different platforms was merged using `minfi::combineArrays`.

### purity/pames.R
Tumor purity was calculated using the `PAMES` package in R. Normal central nervous system samples from the German Cancer Research Center (DKFZ) were used as a control. PAMES operates in three steps. First, AUCs are calculated for each probe discriminating between tumor and normal. Second, a selection of the most informative probes is made. Third, tumor purity is calculated on input samples using these probes.

### meth/train-methylationclasses.R
Glioma methylation subtype classification was performed using L2-regularized logistic regression using the R package `LiblineaR`. Classifiers were trained and evaluated on a set of XXXX common probes from XXX TCGA glioma samples with known methylation subtypes. The classes `LGm6-GBM` and `PA-like` were merged into a single class `LGm6-PA` as the seperation between these classes was based on phenotype. To improve classification accuracy of samples with low tumor purity, DKFZ controls were added to the classifier as seperate classes. Overall accuracy was XX and AUCs for each class ranged from XX to XX. 

### expr/train-transcriptomeclasses.R
Transcriptome classes were predicted using methylation B-values and L2-regularized logistic regression from the `LiblineaR` R package. 