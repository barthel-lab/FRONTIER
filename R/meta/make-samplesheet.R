library(tidyverse)

targ_files = list.files("data/meta", recursive = T, pattern = "Sample_sheet", full.names = T)
targ = lapply(targ_files, read_csv) %>% reduce(full_join)

regex = "^(?:GSM\\d+_)*(\\d{10,12})_(R\\d{2}C\\d{2})_(Grn|Red).idat$"
idat_files = list.files("data/idat", recursive = T, pattern = "idat$", full.names = T)
idat = data.frame(Sentrix_ID = as.numeric(gsub(regex, "\\1", basename(idat_files), perl=T)),
                  Sentrix_Position = gsub(regex, "\\2", basename(idat_files), perl=T),
                  Basename = idat_files,
                  stringsAsFactors = F)
idat = idat %>% filter(!duplicated(paste(Sentrix_ID, Sentrix_Position)))

sheet = targ %>% select(-Basename) %>% left_join(idat) %>%
  mutate(Sentrix_Accession = paste(Sentrix_ID, Sentrix_Position, sep = "_"),
         Sex = ifelse(is.na(Sex), toupper(Gender), Sex),
         Batch = ifelse(is.na(Batch), 1, Batch),
         Sample_Type = ifelse(is.na(Sample_Type), "Initial", Sample_Type)) %>%
  select(Sample_Name, Array, Sentrix_ID, Sentrix_Position, Sentrix_Accession, Dataset, Batch, Basename, Sample_Type,  Age, Sex, Patient, Biopsy, Histology, Grade, IHD_IHC = IDH.IHC)

write.csv(sheet, file="data/meta/Master_Sample_Sheet.csv", row.names = F)
