##################################################
## Project: FRONTIER
## Script purpose: This script creates a master sample sheet (CSV file) that is used by minfi and the methylation pre-processing script
## Date: June 14, 2018
## Author: Floris Barthel
##################################################

library(tidyverse)

setwd(here::here())

## Master sample sheet file name
sheet_file = "data/meta/Master_Sample_Sheet.csv"

## idat file directory
idat_dir = "/projects/verhaak-lab/FRONTIER/data/idat"

## meta file directory
meta_dir = "/projects/verhaak-lab/FRONTIER/data/meta"

## Regular expression for parsing Sentrix ID and Position from filenames
regex = "^(?:GSM\\d+_)*(\\d{10,12})_(R\\d{2}C\\d{2})_(Grn|Red).idat$"

## Find existing sample sheets from various sources
targ_files = list.files(meta_dir, recursive = T, pattern = "Sample_sheet", full.names = T)
targ = lapply(targ_files, read_csv) %>% purrr::reduce(full_join)

## Find .idat files
idat_files = list.files(idat_dir, recursive = T, pattern = "idat$", full.names = T)

## Create a data frame of .idat files and extracted IDs
idat = data.frame(Sentrix_ID = as.numeric(gsub(regex, "\\1", basename(idat_files), perl=T)),
                  Sentrix_Position = gsub(regex, "\\2", basename(idat_files), perl=T),
                  Basename = idat_files,
                  stringsAsFactors = F)
idat = idat %>% filter(!duplicated(paste(Sentrix_ID, Sentrix_Position)))

## Create a master sample sheet by merging idat files and metadata from existing sample sheets
sheet = targ %>% select(-Basename) %>% left_join(idat) %>%
  mutate(Sentrix_Accession = paste(Sentrix_ID, Sentrix_Position, sep = "_"),
         Sex = ifelse(is.na(Sex), toupper(Gender), Sex),
         Batch = ifelse(is.na(Batch), 1, Batch),
         Sample_Type = ifelse(is.na(Sample_Type), "Initial", Sample_Type)) %>%
  select(Sample_Name, Array, Sentrix_ID, Sentrix_Position, Sentrix_Accession, Dataset, Batch, Basename, Sample_Type,  Age, Sex, Patient, Biopsy, Histology, Grade, IHD_IHC = IDH.IHC)

## Write master sample sheet
write.csv(sheet, file=sheet_file, row.names = F)
