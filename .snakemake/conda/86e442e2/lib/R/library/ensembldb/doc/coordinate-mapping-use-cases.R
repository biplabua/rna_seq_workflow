## ----biocstyle, echo = FALSE, results = "asis", message = FALSE------------
library(BiocStyle)
BiocStyle::markdown()

## ----load-libs, message = FALSE--------------------------------------------
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Retrieve the genes
gns <- genes(edb, filter = ~ protein_domain_id == "PF00010" & seq_name == "21")

## ----gene-GRanges----------------------------------------------------------
gns

## ----olig-tx-models, message = FALSE---------------------------------------
## Change chromosome naming style to UCSC
seqlevelsStyle(edb) <- "UCSC"


## ----edb-subset, message = FALSE, echo = FALSE-----------------------------
## Subset the EnsDb to speed up vignette processing
edb <- filter(edb, filter = ~ seq_name %in% c("chr21", "chr16"))

## ----retrieve-olig-tx-models, message = FALSE------------------------------
## Retrieve the transcript models for OLIG1 and OLIG2 that encode the
## the protein domain
txs <- getGeneRegionTrackForGviz(
    edb, filter = ~ genename %in% c("OLIG1", "OLIG2") &
             protein_domain_id == "PF00010")

## ----olig-prot-doms, message = FALSE---------------------------------------
pdoms <- proteins(edb, filter = ~ tx_id %in% txs$transcript &
                           protein_domain_source == "pfam",
                  columns = c("protein_domain_id", "prot_dom_start",
                              "prot_dom_end"))
pdoms

## ----olig-proteinToGenome, message = FALSE---------------------------------
pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                     names = pdoms$protein_id)

pdoms_gnm <- proteinToGenome(pdoms_rng, edb)

## ----olig-proteinToGenome-result, message = FALSE--------------------------
pdoms_gnm


## ----olig-proteinToGenome-reorganize, message = FALSE----------------------
## Convert the list to a GRanges with grouping information
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))

pdoms_gnm_grng

## ----olig-plot, message = FALSE, warning = FALSE, fig.align = "center", fig.width = 8, fig.height = 2, fig.cap = "Transcripts of genes OLIG1 and OLIG2 encoding a helix-loop-helix protein domain. Shown are all transcripts of the genes OLIG1 and OLIG2 that encode a protein with a helix-loop-helix protein domain (PF00010). Genomic positions encoding protein domains defined in Pfam are shown in light blue.", fig.pos = "h!"----
library(Gviz)

## Define the individual tracks:
## - Ideagram
ideo_track <- IdeogramTrack(genome = "hg38", chromosome = "chr21")
## - Genome axis
gaxis_track <- GenomeAxisTrack()
## - Transcripts
gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                              name = "", geneSymbol = TRUE, size = 0.5)
## - Protein domains
pdom_track <- AnnotationTrack(pdoms_gnm_grng, group = pdoms_gnm_grng$grp,
                              id = pdoms_gnm_grng$id, groupAnnotation = "id",
                              just.group = "right", shape = "box",
                              name = "Protein domains", size = 0.5)

## Generate the plot
plotTracks(list(ideo_track, gaxis_track, gene_track, pdom_track))


## ----sim2-fetch, message = FALSE-------------------------------------------
## Fetch all SIM2 transcripts encoding PF00010
txs <- getGeneRegionTrackForGviz(edb, filter = ~ genename == "SIM2" &
                                          protein_domain_id == "PF00010")
## Fetch all Pfam protein domains within these transcripts
pdoms <- proteins(edb, filter = ~ tx_id %in% txs$transcript &
                           protein_domain_source == "pfam",
                  columns = c("protein_domain_id", "prot_dom_start",
                              "prot_dom_end"))


## ----sim2-plot, message = FALSE, warning = FALSE, fig.align = "center", fig.width = 8, fig.height = 2, fig.cap = "Transcripts of the gene SIM2 encoding a helix-loop-helix domain. Shown are all transcripts of SIM2 encoding a protein with a helix-loop-helix protein domain (PF00010). Genomic positions encoding protein domains defined in Pfam are shown in light blue.", echo = FALSE, fig.pos = "h!"----
pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                     names = pdoms$protein_id)
pdoms_gnm <- proteinToGenome(pdoms_rng, edb)

## Convert the list to a GRanges with grouping information
pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))

gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                              name = "", geneSymbol = TRUE, size = 0.5)
## - Protein domains
pdom_track <- AnnotationTrack(pdoms_gnm_grng, group = pdoms_gnm_grng$grp,
                              id = pdoms_gnm_grng$id, groupAnnotation = "id",
                              just.group = "right", shape = "box",
                              name = "Protein domains", size = 0.5)

## Generate the plot
plotTracks(list(ideo_track, gaxis_track, gene_track, pdom_track))


## ----ex2-map, message = FALSE, warning = FALSE-----------------------------
gnm_pos <- GRanges("chr16", IRanges(89920138, width = 1))
prt_pos <- genomeToProtein(gnm_pos, edb)
prt_pos

## ----ex2-select, message = FALSE-------------------------------------------
select(edb, keys = ~ protein_id == names(prt_pos[[1]]), columns = "SYMBOL")

## ----ex2-plot, message = FALSE, warning = FALSE, fig.cap = "Transcripts overlapping, or close to, the genomic position of interest. Shown are all transcripts The genomic position of the variant is highlighted in red.", fig.width = 8, fig.height = 4, fig.pos = "h!"----
## Get transcripts overlapping the genomic position.
txs <- getGeneRegionTrackForGviz(edb, filter = GRangesFilter(gnm_pos))

## Get all transcripts within the region from the start of the most 5'
## and end of the most 3' exon.
all_txs <- getGeneRegionTrackForGviz(
    edb, filter = GRangesFilter(range(txs), type = "within"))

## Plot the data
## - Ideogram
ideo_track <- IdeogramTrack(genome = "hg38", chromosome = "chr16")
## - Genome axis
gaxis_track <- GenomeAxisTrack()
## - Transcripts
gene_track <- GeneRegionTrack(all_txs, showId = TRUE, just.group = "right",
                              name = "", geneSymbol = TRUE, size = 0.5)
## - highlight the region.
hl_track <- HighlightTrack(list(gaxis_track, gene_track), range = gnm_pos)

## Generate the plot
plotTracks(list(ideo_track, hl_track))


## ----ex2-res, message = FALSE----------------------------------------------
## Get the amino acid sequences for the 3 transcripts
prt_seq <- proteins(edb, return.type = "AAStringSet",
                    filter = ~ protein_id == names(prt_pos[[1]]))
## Extract the amino acid at position 294
subseq(prt_seq, start = 294, end = 294)

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

