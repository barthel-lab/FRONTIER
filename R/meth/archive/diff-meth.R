## Script in progress
## Conduct differential methylation analysis

# Get M-values
vumc_m = getM(vumc_g_filt)

## Get Beta values
vumc_b = getBeta(vumc_g_filt)

# Set variables
vumc_targ = vumc_targ %>% mutate(PatientID = factor(sprintf("PT%s", Patient)),
                                 Type = factor(ifelse(M.PA == "Normal", "Normal", ifelse(M.IDH == "IDH mut", "IDHmut", ifelse(M.IDH == "IDH wt", "IDHwt", NA)))))

# use the above to create a design matrix
design = model.matrix(~ 0 + Type + PatientID, data = vumc_targ)

fit = lmFit(vumc_m, design)

contrasts = makeContrasts((TypeIDHwt+TypeIDHmut)-TypeNormal,
                          TypeIDHmut-TypeNormal,
                          TypeIDHwt-TypeNormal,
                          TypeIDHwt-TypeIDHmut,
                          levels=design)
contrasts

# fit the contrasts
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)

annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
ann = annEPIC[match(rownames(vumc_m),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

DMPs = topTable(fit2, num=Inf, coef=2, genelist = ann)

par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(vumc_b, cpg = cpg, pheno = vumc_targ$Type, ylab = "Beta values")
})

library(DMRcate)

my_anno = cpg.annotate(object = vumc_m, datatype = "array", what = "M", 
                       analysis.type = "differential", design = design, 
                       contrasts = TRUE, cont.matrix = contrasts, 
                       coef = "TypeIDHmut - TypeNormal", arraytype = "EPIC")

DMRs = dmrcate(my_anno, lambda=1000, C=2)

# convert the regions to annotated genomic ranges
data(dmrcatedata)
results.ranges <- extractRanges(DMRs, genome = "hg38")

# set up the grouping variables and colours
groups <- pal[1:length(unique(vumc_targ$Type))]
names(groups) <- levels(factor(vumc_targ$Type))
cols <- groups[as.character(factor(vumc_targ$Type))]

pal = brewer.pal(8,"Dark2")

samps <- 1:nrow(vumc_targ)

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges=results.ranges, dmr=1, CpGs=vumc_b, phen.col=cols, what = "Beta",
         arraytype = "EPIC", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg38", samps=samps)
