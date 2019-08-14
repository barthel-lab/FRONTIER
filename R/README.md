# FRONTIER

### cnv/conumee-cnv.R
The R package `conumee` was used to estimate copy number from 450k and EPIC data.

### cnv/cnv-matrix-cluster.R
Sample-only clustering of CNV, across chromosomes.

### expr/runssgsea-toronto.R
We used `ssGSEA` to determine transcriptional clusters (Verhaak 2010, 2013 revision or Wang 2017) for a set of Toronto samples with matching 450k methylation and U133 RNA microarray. 

### expr/toronto-expression.R
See above (consider merging scripts)

### expr/train-transcriptome-classes.R
Transcriptome classes were predicted using methylation B-values and L2-regularized logistic regression from the `LiblineaR` R package. 

### meta/make-master-table.R
Data from various analysis modules (copy number, purity, subtype prediction, etc) were merged into a per-sample table.

### meta/make-sample-sheet.R
Raw Illumina idat files from several sources (VUmc cohort, Toronto cohort, UCSF cohort and DKFZ cohort) were merged and processed individually. 

### meta/plot-metadata-overview.R
An overview of all data was drawn using `ggplot`.

### meth/preprocess.R
Data was processed using the `minfi` packages in R (R Foundation for Statistical Computing, Vienna). Data from the 450k (`IlluminaHumanMethylation450k.ilmn12.hg19`) and EPIC platforms (`IlluminaHumanMethylationEPICanno.ilm10b2.hg19`) were processed seperately. Detection P-values were calculated for each probe and sample, and samples with an average detection P-value > 0.01 were removed from follow-up analysis. Data was normalized using XXXX. Probes on sex chromosomes and known cross-reactive probes were removed, as were probes mapping to known SNPs and probes with a detection P-value > 0.01. Finally, data from different platforms was merged using `minfi::combineArrays`.

### meth/train-methylation-classes.R
Glioma methylation subtype classification was performed using L2-regularized logistic regression using the R package `LiblineaR`. Classifiers were trained and evaluated on a set of XXXX common probes from XXX TCGA glioma samples with known methylation subtypes. The classes `LGm6-GBM` and `PA-like` were merged into a single class `LGm6-PA` as the seperation between these classes was based on phenotype. To improve classification accuracy of samples with low tumor purity, DKFZ controls were added to the classifier as seperate classes. Overall accuracy was XX and AUCs for each class ranged from XX to XX. 

### purity/pames.R
Tumor purity was calculated using the `PAMES` package in R. Normal central nervous system samples from the German Cancer Research Center (DKFZ) were used as a control. PAMES operates in three steps. First, AUCs are calculated for each probe discriminating between tumor and normal. Second, a selection of the most informative probes is made. Third, tumor purity is calculated on input samples using these probes.