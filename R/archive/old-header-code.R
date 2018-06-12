### old header code

load('results/MSG.QC.filtered.normalized.anno.filt_only.Rdata')
load('results/MSG.master_meta.Rdata')

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