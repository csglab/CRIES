# setwd('/project/6001417/alkallr/alkallas_projects/gen_intronic_exonic_annot/code/')

library(rtracklayer)
library(data.table)
library(magrittr)

### User arguments ####
args <- commandArgs(trailingOnly = TRUE)
# args <- c('../../shared_data/annotations/Ensembl/Homo_sapiens.GRCh38.101.gtf.gz', 'Ensembl', '<output.directory>')

### Load exon annotation ####
exons <- import(con = args[1], feature = 'exon') %>% sort
exons <- as.data.table(exons)

exons[ , c('gene_id', 'transcript_id') := list(
    paste(gene_id, gene_version, sep = '.'),
    paste(transcript_id, transcript_version, sep = '.')
), ]

### Subset constitutive ensembl_havana exons ####
exons[ , exon_id := paste0(gene_id, '_', start, '-', end), ] # create exon id

# calculate the number of transcripts per exon and the number of transcripts per gene
exons[source == 'ensembl_havana', nTxPerExon := transcript_id %>% unique %>% length, by = 'exon_id'] # N Tx per exon
exons[source == 'ensembl_havana', nTxPerGene := transcript_id %>% unique %>% length, by = 'gene_id'] # N Tx per gene
exons[ , nTx := transcript_id %>% unique %>% length, by = 'gene_id'] # N Tx

# constitutive exons are where nTxPerExon == nTxPerGene
consExons <- makeGRangesFromDataFrame(
    unique(exons[nTxPerExon == nTxPerGene, c(colnames(exons)[1:5], 'gene_id', 'gene_name'), with = F]),
    keep.extra.columns = T
)

consExons$type <- 'exon'
consExons$subtype <- 'consExon'

### Subset constitutive ensembl_havana intronic regions ####
exons <- makeGRangesFromDataFrame(exons, keep.extra.columns = T)

# premature transcript ranges
prematureTx <- split(
    exons[exons$source == 'ensembl_havana'],
    exons[exons$source == 'ensembl_havana']$transcript_id
) %>% range %>% stack(., 'transcript_id')

prematureTx$gene_id <- exons[match(prematureTx$transcript_id, exons$transcript_id), ]$gene_id

# get all exonic regions overlapping pre-RNA
prematureTx.exon.overlaps <- findOverlapPairs(
    query = prematureTx, subject = exons,
    ignore.strand = T
)

# purge all exonic regions from pre-RNA
introns <- setdiff(

    x = split(prematureTx, prematureTx$transcript_id)[prematureTx$transcript_id],

    y = split(
        S4Vectors::second(prematureTx.exon.overlaps), S4Vectors::first(prematureTx.exon.overlaps)$transcript_id
    )[prematureTx$transcript_id],

    ignore.strand = T

) %>% stack(x = ., 'transcript_id')

# get number of transcripts per intronic region per gene using coverage
introns <- as.data.table(introns)
introns[ , seqnames := exons[match(transcript_id, exons$transcript_id), ]$gene_id]

# TODO find out how to restrict coverage to specific regions to avoid unnecessary zeroes
consIntrons <- makeGRangesFromDataFrame(introns, keep.extra.columns = T) %>%
    coverage %>% as(., 'GRanges') %>% as.data.table

consIntrons <- consIntrons[score > 0]

consIntrons <- base::cbind(
    copy(consIntrons)[ , c('seqnames', 'strand') := list(
        exons[match(consIntrons$seqnames, exons$gene_id), ] %>% seqnames %>% as.character,
        exons[match(consIntrons$seqnames, exons$gene_id), ] %>% strand %>% as.character
    )
    ],
    mcols(exons[exons$source == 'ensembl_havana'])[
        match(consIntrons$seqnames, exons[exons$source == 'ensembl_havana']$gene_id),
        c('gene_id', 'gene_name', 'nTxPerGene')
        ] %>% as.data.table
)

# constitutive intronic regions are where coverage == nTxPerGene
consIntrons <- makeGRangesFromDataFrame(
    unique(consIntrons[score == nTxPerGene, c(colnames(consIntrons)[1:5], 'gene_id', 'gene_name'), with = F]),
    keep.extra.columns = T
)

consIntrons$type <- 'intron'
consIntrons$subtype <- 'consIntron'

### Initialize output ####
# which now contains only the constitutive regions of genes that have either
# constitive intronic regions and/or constitutive exons
out <- c(consExons, consIntrons)

### Output ####
out <- sort(out)

out$gene_biotype <- exons[match(out$gene_id, exons$gene_id), ]$gene_biotype
out$gene_id <- paste(out$gene_id, out$gene_name, sep = '_')

seqlevelsStyle(out) <- args[2]

fpath <- head(strsplit(args[1], '/')[[1]], -1) %>% paste(., collapse = '/')
fname <- paste0(tail(strsplit(args[1], '/')[[1]], 1), '_', Sys.Date(), '_constExonsAndIntrons_', args[2], '.gtf')

export(
    object = out[ , -2], # omit gene name
    con = file.path(args[3], fname)
)

w <- warnings()
print(w)

message('Annotations:')
print(out)

message('Done!')

# ### Plot for MMB chapter ####
# library(ggplot2)
#
# chosen.one <- 'ENSG00000250067.12'
#
# to.plot <- exons[exons$transcript_id %in% exons[exons %over% range(exons[exons$gene_id == chosen.one]), ]$transcript_id, ]
#
# # introns
# to.plot.introns <- split(to.plot, to.plot$transcript_id) %>% range %>% stack(., 'transcript_id')
#
# elementMetadata(to.plot.introns) <- elementMetadata(exons)[
#     match(to.plot.introns$transcript_id, exons$transcript_id),
#     c('transcript_id', 'transcript_name', 'gene_id', 'gene_name', 'transcript_biotype', 'source')
# ]
#
# to.plot.introns <- to.plot.introns[order(to.plot.introns$source, order(to.plot.introns), decreasing = T), ]
#
# to.plot.introns <- as.data.table(to.plot.introns)
# to.plot.introns <- to.plot.introns[ , y := 1:length(transcript_id), by = 'gene_id']
#
# # exons
# to.plot.exons <- to.plot[ ,
#     c('transcript_id', 'transcript_name', 'gene_id', 'gene_name', 'transcript_biotype', 'source')
# ]
#
# to.plot.exons$exon.width <- width(to.plot.exons)
#
# to.plot.exons <- as.data.table(to.plot.exons)
# to.plot.exons[ , y := to.plot.introns[match(to.plot.exons$transcript_id, to.plot.introns$transcript_id), ]$y, ]
# to.plot.exons[ , c('yTop', 'yBtm') := list(y + 0.3, y - 0.3), ]
#
# table(to.plot.exons$strand)
# #  +  -  *
# # 39  0  0
#
# # output anot
# to.plot.out <- as.data.table(out[out %over% out[grep(chosen.one, out$gene_id), ], ])
# to.plot.out[ , yTop := -0.3]
# to.plot.out[ , yBtm := +0.3]
# to.plot.out[grepl('consIntron', subtype), transcript_name := paste0(gene_name, '\n', 'cost. introns'), ]
# to.plot.out[grepl('consExon', subtype), transcript_name := paste0(gene_name, '\n', 'cost. exons'), ]
#
# print(to.plot.out)
#
# tx.lev <- c(
#     to.plot.out$transcript_name,
#     sort(to.plot.exons$transcript_name, decreasing = T)
# ) %>% unique
#
# print(tx.lev)
#
# tx.lev <- c(
#     "YJEFN3\ncost. introns", "YJEFN3\ncost. exons", "YJEFN3-203", "YJEFN3-201", "YJEFN3-202",
#     "NDUFA13-203", "AC011448.1-201", "AC011448.1-202"
# )
#
# to.plot.out$transcript_name <- factor(to.plot.out$transcript_name, tx.lev)
#
# to.plot.exons$transcript_name <- factor(to.plot.exons$transcript_name, tx.lev)
#
# to.plot.introns$transcript_name <- factor(to.plot.introns$transcript_name, tx.lev)
#
# # delim
# to.plot.out.delim <- lapply(
#     tx.lev,
#     function(i) {
#         out <- copy(to.plot.out)
#         out[ , transcript_name := i]
#     }
# ) %>% rbindlist
#
# to.plot.out.delim <- melt(
#     to.plot.out.delim,
#     measure.vars = c('start', 'end')
# )
#
# to.plot.out.delim$transcript_name <- factor(to.plot.out.delim$transcript_name, tx.lev)
#
# pdf('../output/annotation_figure_for_bookChapter_new.pdf', width = 7, height = 4.5)
#
# ggp <- ggplot() +
#     # introns
#     geom_segment(mapping = aes(y = y, yend = y, x = start, xend = end, color = source), data = to.plot.introns, inherit.aes = F, size = 1.0) +
#     # exons
#     geom_rect(mapping = aes(xmin = start, xmax = end, ymin = yBtm, ymax = yTop, fill = source), data = to.plot.exons, inherit.aes = F) +
#     # const introns and exons
#     geom_rect(mapping = aes(xmin = start, xmax = end, ymin = yBtm, ymax = yTop), data = to.plot.out[grep(chosen.one, gene_id), ], inherit.aes = F, fill = 'black') +
#     # delim
#     geom_vline(mapping = aes(xintercept = value), data = to.plot.out.delim, linetype = 2, color = 'grey30', size = 0.5) +
#
#     facet_grid(transcript_name ~ ., scales = 'free_y', space = 'free_y') +
#
#     xlab('Chromosome 19 GRCh38 coordinates') +
#
#     theme_minimal() +
#
#     theme(
#         panel.grid = element_blank(), axis.text.y = element_blank(),
#         axis.title.y = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0),
#         axis.title.x = element_text(color = 'black'), axis.text.x = element_text(color = 'black'),
#         axis.ticks.x = element_line(color = 'black'), axis.ticks.length.x = unit(0.5, 'line'), legend.position = "top"
#     ) +
#
#     coord_cartesian(
#         xlim = c(
#             19528300,
#             max(to.plot.introns$end)
#         )
#     ) +
#
#     scale_fill_manual(
#         values = c(havana = 'darkred', 'ensembl_havana' = 'goldenrod'), name = NULL
#     ) +
#
#     scale_color_manual(
#         values = c(havana = 'darkred', 'ensembl_havana' = 'goldenrod'), name = NULL
#     )
#
# print(ggp)
#
# dev.off()
