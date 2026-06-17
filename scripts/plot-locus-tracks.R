#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(foreach)
library(ggbio)


mysplit <- function(x) {
    paste(unique(unlist(strsplit(x, split=';'))), collapse=';')
}

if(!file.exists('methQTL/mean-betas.tsv.gz')) {
    # collect beta values for all donors
    betas <- foreach(celltype=c('PBMC','IPSC'), .combine='rbind') %do% {
        dat.tmp <- fread(paste0('MEFFIL/',celltype,'.beta.tsv'))
        dat.tmp <- melt(dat.tmp, measure.vars=grep('^NIH', colnames(dat.tmp), value=T), variable.name='Donor',value.name='beta')
        dat.tmp[, 'celltype' := celltype]
        return(dat.tmp[])
    }
    # Calculate mean beta value per probe per cell type
    mean_betas <- betas[, list('mean_beta'=mean(beta)), by=list(POS,celltype)]

    # Add in position and gene annotation
    EPIC.anno <- fread('DATA/EPIC.anno.GRCh38.tsv', select=c('probeID','chrm','start','GeneNames'))
    setkey(EPIC.anno, probeID)
    setkey(mean_betas, POS)
    mean_betas <- merge(EPIC.anno, mean_betas, by.x='probeID', by.y='POS')

    setkey(mean_betas, chrm, start, celltype)
    mean_betas[celltype=='IPSC', celltype := 'iPSC']
    mean_betas <- dcast(mean_betas, probeID + chrm + start + GeneNames ~ celltype, value.var='mean_beta')   
    rm(betas); gc()
    # Get EPIC annotation gene list from duplicat -> unique list per row
    mean_betas[, 'Genes' := apply(.SD, 1, mysplit), .SDcols='GeneNames']
    mean_betas[Genes == '.', Genes := NA]
    mean_betas[, GeneNames := NULL]
    setcolorder(mean_betas, c('probeID','Genes','chrm','start','PBMC','iPSC'))

    # Exclude probes only passing QC in a single group
    mean_betas <- mean_betas[!is.na(PBMC) & ! is.na(iPSC)]
    fwrite(mean_betas, 'methQTL/mean-betas.tsv.gz', sep='\t')
} else {
    mean_betas <- fread('methQTL/mean-betas.tsv.gz')
}



get_region <- function(DT, CHR, START, STOP) {
    return(DT[chrm==CHR & start %between% c(START, STOP)])
}

plot_region_twotracks <- function(DT) {
    ggplot(DT, aes(x=start, xend=start, y=0, yend=mean_beta)) +
        facet_grid(celltype~.) +
        geom_segment() +
        labs(x='POS', y='% Methylation') +
        geom_point(data=DT[!is.na(eQTL)], aes(x=start, y=mean_beta, fill=eQTL), shape=21, size=3) +
        scale_fill_identity() +
        theme_few() +
        scale_y_continuous(breaks=c(0,1)) +
        scale_x_continuous()
}

mean_betas.long <- melt(mean_betas, measure.vars=c('PBMC','iPSC'), variable.name='celltype', value.name='mean_beta')
mean_betas.long[, celltype := factor(celltype, levels=c('PBMC','iPSC'))]

ipsc_signif <- fread('methQTL/IPSC.cis_qtl_significant.txt')$phenotype_id
pbmc_signif <- fread('methQTL/PBMC.cis_qtl_significant.txt')$phenotype_id
ipsc_unique <- setdiff(ipsc_signif, pbmc_signif)
pbmc_unique <- setdiff(pbmc_signif, ipsc_signif)

mean_betas.long[celltype == 'iPSC' & probeID %in% ipsc_signif, eQTL := 'white']
mean_betas.long[celltype == 'iPSC' & probeID %in% ipsc_unique, eQTL := '#39b3ffff']
mean_betas.long[celltype == 'PBMC' & probeID %in% pbmc_signif, eQTL := 'white']
mean_betas.long[celltype == 'PBMC' & probeID %in% pbmc_unique, eQTL := '#ec137bff']

library(ggbio)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene

tx <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
txdb <- stack(tx)
txdb <- subset(txdb, seqnames %in% c('chrX','chrY',paste0('chr',1:22)))
txdb <- as.data.table(txdb)
setDT(txdb)
txdb[, SYMBOL := unlist(txdb$SYMBOL)]

o <- foreach(symbol=unique(txdb$SYMBOL), .combine='rbind', .errorhandling='remove') %do% {
    dt.tmp <- txdb[SYMBOL==symbol]
    minval <- min(dt.tmp$start)
    maxval <- max(dt.tmp$end)
    chr <- unique(dt.tmp$seqnames)
    strand <- unique(dt.tmp$strand)
    width <- 1 + maxval - minval
    data.table('seqnames'=chr, 'start'=minval, 'end'=maxval, 'width'=width, 'strand'=strand, 'SYMBOL'=symbol)
}

# Convert to gene symbols
tx2 <- makeGRangesFromDataFrame(o, keep.extra=TRUE)
genesymbol <- setNames(tx2, tx2$SYMBOL)




get_plots <- function(min_x, max_x, genename, chrom) {
    wh <- genesymbol[c(genename)]
    wh <- range(wh, ignore.strand = TRUE)
    g.methylation <- plot_region_twotracks(get_region(mean_betas.long, chrom ,min_x, max_x)) + labs(y=paste0(genename, ' locus')) + scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0))
    g.track <- autoplot(txdb, which = wh) + scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0)) + theme_few() + labs(y=paste0(genename, ' locus'), x=paste0(chrom, ' POS'))
    return(list(g.methylation, g.track@ggplot))
}

get_locus_plots <- function(genename, genesymbol) {
    wh <- genesymbol[genename]
    min_x <- min(start(wh))
    max_x <- max(end(wh))
    xtrawidth <- round(max_x - min_x)*0.03
    min_x <- min_x - xtrawidth
    max_x <- max_x + xtrawidth
    chrom <- as.character(wh@seqnames@values)
    wh <- range(wh, ignore.strand = TRUE)
    g.methylation <- plot_region_twotracks(get_region(mean_betas.long, chrom ,min_x, max_x)) + labs(y=paste0(genename, ' locus')) + scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0))
    g.track <- autoplot(Homo.sapiens, which = wh) + scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0)) + theme_few() + labs(y=paste0(genename, ' locus'), x=paste0(chrom, ' POS'))
    return(list(g.methylation, g.track@ggplot))
}


# get locus plots of interest
ADARB2 <- get_locus_plots(genename='ADARB2', genesymbol)
B3GNTL1 <- get_locus_plots(genename='B3GNTL1', genesymbol)
HLADPB2 <- get_locus_plots(genename='HLA-DPB2', genesymbol)
SNTG2 <- get_locus_plots(genename='SNTG2', genesymbol)

# Merge into singular figure
g.all <- plot_grid(ADARB2[[1]], ADARB2[[2]], B3GNTL1[[1]], B3GNTL1[[2]],HLADPB2[[1]], HLADPB2[[2]],SNTG2[[1]], SNTG2[[2]], ncol=1, align='v', axis='lr')
ggsave(g.all, file='PLOTS/METHQTL/Fig6.png', width=50, height=50, units='cm')
ggsave(g.all, file='PLOTS/METHQTL/Fig6.svg', width=50, height=50, units='cm')
ggsave(g.all, file='PLOTS/METHQTL/Fig6.pdf', width=50, height=50, units='cm')

