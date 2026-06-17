#!/usr/bin/env Rscript

library(data.table)



# |   Study   |  N rows   | Bonferroni-adjusted threshold   |
# |-----------|-----------|-----------|
# | QTD000361 | 143317269 | `0.05/143317269` = `3.488763e-10` |
# | QTD000366 | 142812959 | `0.05/142812959` = `3.501083e-10` |
# | QTD000399 | 162075245 | `0.05/162075245` = `3.084987e-10` |

args <- commandArgs(trailingOnly=TRUE)


in_filename <- args[1]
studyid <- unlist(strsplit(in_filename, split='\\.'))[1]

out_filename <- paste0(studyid, '.bonferroni.tsv')

threshold <- as.numeric(args[2])
header_cols <-  c('molecular_trait_id','chromosome','position','ref','alt','variant','ma_samples','maf','pvalue','beta','se','type','ac','an','r2','molecular_trait_object_id','gene_id','median_tpm','rsid')


dat <- fread(in_filename, header=F)
setnames(dat, header_cols)


dat <- dat[pvalue < threshold][]

dat <- dat[type == 'SNP']

fwrite(dat, file=out_filename, quote=F, row.names=F, col.names=T, sep='\t')

