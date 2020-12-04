## ----biocstyle, echo = FALSE, results = "asis", message = FALSE------------
library(BiocStyle)
BiocStyle::markdown() 

## ----load-libs, message = FALSE--------------------------------------------
library(ensembldb)
library(EnsDb.Hsapiens.v86)

edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == "X") 

## ----genomeToTranscript-define---------------------------------------------
gnm <- GRanges("X:107716399-107716401")

## ----genomeToTranscript-ex1-plot, message = FALSE--------------------------
library(Gviz)
## Since we're using Ensembl chromosome names we have to set:
options(ucscChromosomeNames = FALSE)

## Define a genome axis track
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))
gtx <- GeneRegionTrack(gnm_gns, name = "tx", geneSymbol = TRUE,
                       showId = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm)
## plot the region
plotTracks(list(ht))
 

## ----genomeToTranscript-ex1-map, message = FALSE---------------------------
## Map genomic coordinates to within-transcript coordinates
gnm_tx <- genomeToTranscript(gnm, edbx) 

## ----genomeToTranscript-ex1-object-----------------------------------------
gnm_tx 

## ----genomeToTranscript-ex2, message = FALSE-------------------------------
gnm_1 <- gnm
strand(gnm_1) <- "-"
gnm_2 <- gnm
strand(gnm_2) <- "+"
gnm <- c(gnm_1, gnm_2)

genomeToTranscript(gnm, edbx)

## ----genomeToProtein-ex1, message = FALSE----------------------------------
gnm <- GRanges("X", IRanges(start = c(630898, 644636, 644633, 634829),
                            width = c(5, 1, 1, 3)))
gnm_prt <- genomeToProtein(gnm, edbx)
 

## ----genomeToProtein-ex1-res1----------------------------------------------
gnm_prt[[1]]

## ----genomeToProtein-ex1-res2----------------------------------------------
gnm_prt[[2]]

## ----genomeToProtein-ex1-res3----------------------------------------------
gnm_prt[[3]] 

## ----genomeToProtein-ex1-res3-2, message = FALSE---------------------------
prt <- proteins(edbx, filter = ProteinIdFilter(names(gnm_prt[[3]])))

nchar(prt$protein_sequence) 

## ----genomeToProtein-ex1-res4----------------------------------------------
gnm_prt[[4]] 

## ----proteinToTranscript-ex1, message = FALSE------------------------------
GAGE10 <- proteins(edbx, filter = ~ genename == "GAGE10")
GAGE10

## Define the IRanges object.
GAGE10_prt <- IRanges(start = 5, end = 9, names = GAGE10$protein_id)

## ----proteinToTranscript-ex1-map, message = FALSE--------------------------
GAGE10_tx <- proteinToTranscript(GAGE10_prt, edbx) 

## ----proteinToTranscript-ex1-res-------------------------------------------
GAGE10_tx

## ----proteinToTranscript-ex2, message = FALSE------------------------------
ids <- c("O15266", "Q9HBJ8", "donotexistant")
prt <- IRanges(start = c(13, 43, 100), end = c(21, 80, 100))
names(prt) <- ids

prt_tx <- proteinToTranscript(prt, edbx, idType = "uniprot_id") 

## ----proteinToTranscript-ex2-res1------------------------------------------
prt_tx[[1]] 

## ----proteinToTranscript-ex2-res2------------------------------------------
prt_tx[[2]] 

## ----proteinToTranscript-ex2-res3------------------------------------------
prt_tx[[3]] 

## ----proteinToGenome-gage10-define, message = FALSE------------------------
## Define the IRanges object.
GAGE10_prt <- IRanges(start = 5, end = 9, names = "ENSP00000385415")
 

## ----proteinToGenome-gage10-map, message = FALSE---------------------------
GAGE10_gnm <- proteinToGenome(GAGE10_prt, edbx) 

## ----proteinToGenome-gage10-res--------------------------------------------
GAGE10_gnm 

## ----proteinToGenome-uniprot-ids, message = FALSE--------------------------
## Define the IRanges providing Uniprot IDs.
uni_rng <- IRanges(start = c(2, 12, 8), end = c(2, 15, 17),
                   names = c("D6RDZ7", "O15266", "H7C2F2"))

## We have to specify that the IDs are Uniprot IDs
uni_gnm <- proteinToGenome(uni_rng, edbx, idType = "uniprot_id") 

## ----proteinToGenome-uniprot-cds_ok----------------------------------------
uni_gnm[[3]]

## ----proteinToGenome-uniprot-counts----------------------------------------
## To how many Ensembl proteins was each Uniprot ID mapped?
lengths(uni_gnm) 

## ----proteinToGenome-uniprot-multi-----------------------------------------
uni_gnm[["O15266"]] 

## ----proteinToGenome-SYP-fetch-domains, message = FALSE--------------------
SYP <- proteins(edbx, filter = ~ genename == "SYP",
                columns = c("protein_id", "tx_id",
                            listColumns(edbx, "protein_domain")),
                return.type = "AAStringSet")

SYP

## ----proteinToGenome-SYP-single-protein, message = FALSE-------------------
## How many proteins are annotated to SYP?
unique(mcols(SYP)$protein_id)

## Reduce the result to a single protein
SYP <- SYP[names(SYP) == "ENSP00000263233"]

## List the available protein domains and additional annotations
mcols(SYP) 

## ----proteinToGenome-SYP-map, message = FALSE------------------------------
SYP_rng <- IRanges(start = mcols(SYP)$prot_dom_start,
                   end = mcols(SYP)$prot_dom_end)
mcols(SYP_rng) <- mcols(SYP)

## Map the domains to the genome. We set "id" to the name
## of the metadata columns containing the protein IDs
SYP_gnm <- proteinToGenome(SYP_rng, edbx, id = "protein_id") 

## ----proteinToGenome-SYP-second--------------------------------------------
SYP_gnm[[2]] 

## ----proteinToGenome-SYP-plot, message = FALSE-----------------------------
library(Gviz)

## Define a genome axis track
gat <- GenomeAxisTrack()

## Get the transcript ID:
txid <- SYP_gnm[[1]]$tx_id[1]

## Get a GRanges for the transcript
trt <- getGeneRegionTrackForGviz(edbx, filter = TxIdFilter(txid))

## Define a GRanges for the mapped protein domains and add
## metadata columns with the grouping of the ranges and the
## IDs of the corresponding protein domains, so they can be
## identified in the plot
dmns <- unlist(GRangesList(SYP_gnm))
dmns$grp <- rep(1:length(SYP_rng), lengths(SYP_gnm))
dmns$id <- rep(mcols(SYP_rng)$protein_domain_id, lengths(SYP_gnm))

## Since we're using Ensembl chromosome names we have to set
options(ucscChromosomeNames = FALSE)

## Plotting the transcript and the mapped protein domains.
plotTracks(list(gat,
                GeneRegionTrack(trt, name = "tx"),
                AnnotationTrack(dmns, group = dmns$grp,
                                id = dmns$id,
                                groupAnnotation = "id",
                                just.group = "above",
                                shape = "box",
                                name = "Protein domains")),
           transcriptAnnotation = "transcript")
 

## ----transcriptToGenome-map, message = FALSE-------------------------------
rng_tx <- IRanges(start = c(501, 1), width = c(5, 5),
                  names = c("ENST00000486554", "ENST00000381578"))

rng_gnm <- transcriptToGenome(rng_tx, edbx) 

## ----transcriptToGenome-res-1----------------------------------------------
rng_gnm

## ----pkp2-cdsToTranscript--------------------------------------------------
## Define the position within the CDS of the transcript
pkp2_cds <- IRanges(start = c(1643, 1881), width = c(1, 1),
                    name = rep("ENST00000070846", 2))

## Convert cds-relative to transcript-relative coordinates
pkp2 <- cdsToTranscript(pkp2_cds, EnsDb.Hsapiens.v86)

pkp2 

## ----pkp2-transcriptToGenome-----------------------------------------------
pkp2_gnm <- transcriptToGenome(pkp2, EnsDb.Hsapiens.v86)

pkp2_gnm 

## ----pkp2-variant-pos-validate---------------------------------------------
library(BSgenome.Hsapiens.NCBI.GRCh38)

getSeq(BSgenome.Hsapiens.NCBI.GRCh38, pkp2_gnm) 

## ----transcriptToPrptein-map, message = FALSE------------------------------
rng_tx <- IRanges(start = c(501, 1, 200), width = c(5, 5, 4),
                  names = c("ENST00000486554", "ENST00000381578",
                            "ENST00000431238"))
rng_prt <- transcriptToProtein(rng_tx, edbx) 

## ----transcriptToProtein-res-----------------------------------------------
rng_prt 

## ----sessionInfo-----------------------------------------------------------
sessionInfo() 

