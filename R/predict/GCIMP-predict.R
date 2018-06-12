
setwd("~/projects/MSG/")

load('results/MSG.master_meta.Rdata')

## 7-probe G-CIMP low predictive signature
## Samples wher
## prognostication value and FDR were assigned at n ≥ 5 probes
probes = c("cg09732711" = 0.7, "cg09326832" = 0.28, "cg24665265" = 0.67, "cg06220958" = 0.17, "cg10245915" = 0.12, "cg11689625" = 0.31, "cg11799650" = 0.4)

load('results/MSG.QC.filtered.normalized.anno.filt_only.Rdata')

vumc_targ = vumc_targ %>% mutate(M.PA = ifelse(PC1 > -0.5 | PC2 > 0.5, "Tumor", "Normal"), 
                                 M.PAI = ifelse(PC1 > 1 | PC2 > 0.5, "Tumor", ifelse(PC1 > -0.5, "Intermediate", "Normal")),
                                 M.IDH = ifelse(PC1 > -0.5 & M.PA == "Tumor", "IDH mut", ifelse(PC2 > 0.5 & M.PA == "Tumor", "IDH wt", NA)))

tmp = pData(vumc_g_filt) %>% as.data.frame() %>% left_join(vumc_targ) %>% DataFrame()
rownames(tmp) = rownames(pData(vumc_g_filt))

pData(vumc_g_filt) = tmp

pData(vumc_g_filt)$DataSet = "VUmc"
pData(trnt_g_filt)$DataSet = "Toronto"
pData(ucsf_g_filt)$DataSet = "UCSF"

all_data = combineArrays(vumc_g_filt, trnt_g_filt, outType = "IlluminaHumanMethylation450k")
all_data = combineArrays(all_data, ucsf_g_filt, outType = "IlluminaHumanMethylation450k")

pData(all_data) = pData(all_data)[,1:3] %>% as.data.frame() %>% left_join(meta) %>% DataFrame()
colnames(all_data) = all_data$Sentrix_Accession

########

bvals = getBeta(all_data)

tmp = bvals[intersect(names(probes), rownames(bvals)), ]

probes_tmp = matrix(probes[match(rownames(tmp), names(probes))], ncol=1)

gcimplow_pred = apply(apply(tmp, 2, function(x) x - probes_tmp) < 0, 2, sum)

