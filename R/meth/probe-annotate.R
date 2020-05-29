##################################################
## Project: FRONTIER
## Script purpose: Annotate DNA methylation probes
## Created: May 11, 2020
## Updated: May 13, 2020
## Author: Floris Barthel
##################################################

library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # For official Illumina 450k probe metadata

library(rtracklayer)                                  # For UCSC genome browser tables
library(BSgenome.Hsapiens.UCSC.hg19)                  # Genome FASTA
library(VariantAnnotation)                            # For mapping CpGs to positions in transcripts
library(GenomicFeatures)                              # For loading TxDb objects
library(tidyverse)

## Load GENCODE GTF file
txdb_hg38 = makeTxDbFromGFF('sandbox/gencode.v22.annotation.gtf.gz') ## NOTE: hg38

## Load probe:gene mapping and probe TSS mapping
## From: http://zwdzwd.github.io/InfiniumAnnotation
## Note: also hg38

#manifest_hg38 <- read_tsv('sandbox/hm450.hg38.manifest.tsv')
#manifest_hg19 <- read_tsv('sandbox/hm450.hg19.manifest.tsv')

manifest_gencode_hg38 <- read_tsv('sandbox/hm450.hg38.manifest.gencode.v22.tsv', 
                                  col_types = cols('CpG_chrm' = col_character(),
                                                   'CpG_beg' = col_integer(), 
                                                   'CpG_end' = col_integer(), 
                                                   'distToTSS' = col_character(), 
                                                   'geneNames' = col_character(), 
                                                   'transcriptTypes' = col_character(), 
                                                   'transcriptIDs' = col_character())) %>%
  filter(complete.cases(CpG_chrm, CpG_beg, CpG_end))

manifest_gencode_hg38_gr <- GRanges(seqnames = manifest_gencode_hg38$CpG_chrm, IRanges(manifest_gencode_hg38$CpG_beg + 1, manifest_gencode_hg38$CpG_beg + 1), probe = manifest_gencode_hg38$probeID)

manifest_gencode_hg38_tx <- manifest_gencode_hg38 %>% 
  separate_rows(geneNames, transcriptIDs, transcriptTypes, distToTSS, sep = ";") %>%
  select(probe = probeID, gene = geneNames, transcript_id = transcriptIDs, transcript_type = transcriptTypes, dist_to_tss = distToTSS) %>%
  mutate(dist_to_tss = as.integer(dist_to_tss))

## Map to transcriptome
manifest_gencode_hg38_gr_txmap <- locateVariants(manifest_gencode_hg38_gr, txdb_hg38, AllVariants(), ignore.strand = TRUE)

manifest_gencode_hg38_gr_txmap <- manifest_gencode_hg38_gr_txmap %>% as_tibble() %>% 
  mutate(probe = manifest_gencode_hg38_gr$probe[QUERYID],
         transcript_id = transcripts(txdb_hg38)$tx_name[as.integer(TXID)],
         location = LOCATION) %>%
  select(probe, transcript_id, location) %>%
  mutate(location = factor(location, levels = c("promoter", "fiveUTR", "coding", "spliceSite", "threeUTR", "intron", "intergenic"))) %>% ## Sorted in relative order of importance for methylation
  arrange(location)

## Remove rows that have multiple locations for a single probe:transcript
idx = which(duplicated(sprintf("%s%s", manifest_gencode_hg38_gr_txmap$probe, manifest_gencode_hg38_gr_txmap$transcript_id)))
manifest_gencode_hg38_gr_txmap_filtered <- manifest_gencode_hg38_gr_txmap[-idx,]

## Annotate probe:transcript mappings with location info
manifest_gencode_hg38_tx_anno <- manifest_gencode_hg38_tx %>% left_join(manifest_gencode_hg38_gr_txmap_filtered)

## Get UCSC genome browser gap table for telomeres and centromeres
session <- browserSession()
genome(session) <- "hg19"

tcmap <- getTable(ucscTableQuery(session, track="gap", table="gap")) %>% 
  filter(type %in% c("telomere","centromere")) %>%
  group_by(chrom) %>%
  arrange(chromStart) %>%
  summarize(seqlength = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[unique(chrom)],
            n_tel = sum(type=="telomere"),
            n_cent = sum(type=="centromere"),
            p_start = chromEnd[type=="telomere"][1],
            centromere_start = chromStart[type=="centromere"],
            centromere_end = chromEnd[type=="centromere"],
            q_end = chromStart[type=="telomere"][2]) %>%
  ungroup() %>%
  mutate(chrom = as.character(chrom),
         p_start = ifelse(is.na(p_start), as.integer(0), as.integer(p_start)),
         q_end = ifelse(is.na(q_end), seqlength, q_end))

## Illumina probe metadata
ann450k_hg19 = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as_tibble()

probemap <- ann450k_hg19 %>%
  select(probe = Name, chrom = chr, pos, rel_island = Relation_to_Island) %>%
  left_join(tcmap) %>%
  mutate(arm = case_when(pos < centromere_start ~ "p",
                         pos > centromere_end ~ "q",
                         TRUE ~ NA_character_),
         arm_size = case_when(pos < centromere_start ~ centromere_start - p_start,
                              pos > centromere_end ~ q_end - centromere_end,
                              TRUE ~ NA_integer_),
         dist_to_tel = case_when(pos < centromere_start ~ pos - p_start,
                                 pos > centromere_end ~ q_end - pos,
                                 TRUE ~ NA_integer_),
         dist_to_cent = case_when(pos < centromere_start ~ centromere_start - pos,
                                  pos > centromere_end ~ pos - centromere_end,
                                  TRUE ~ NA_integer_),
         rel_dist_to_tel = dist_to_tel / arm_size,
         rel_dist_to_cent = dist_to_cent / arm_size) %>%
  select(probe, chrom, pos, arm, arm_size, dist_to_tel, dist_to_cent, rel_dist_to_tel, rel_dist_to_cent, rel_island)

nearest_protein_coding_tss <- manifest_gencode_hg38_tx_anno %>%
  filter(transcript_type == "protein_coding" | is.na(transcript_id)) %>%
  arrange(abs(dist_to_tss)) %>%
  group_by(probe) %>%
  #summarize(gene = gene[1], dist_to_tss = dist_to_tss[1]) %>%
  filter(row_number() == 1) %>%
  ungroup()

probemap <- probemap %>% left_join(nearest_protein_coding_tss)

save(probemap, manifest_gencode_hg38_gr_txmap_filtered, file = "sandbox/meth-probe-anno.Rdata")

## END ##
