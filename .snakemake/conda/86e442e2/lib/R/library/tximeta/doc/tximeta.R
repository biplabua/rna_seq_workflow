## ------------------------------------------------------------------------
dir <- system.file("extdata/salmon_dm", package="tximportData")
files <- file.path(dir, "SRR1197474", "quant.sf") 
file.exists(files)
coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
coldata

## ----eval=FALSE----------------------------------------------------------
#  library(tximeta)
#  se <- tximeta(coldata)

## ------------------------------------------------------------------------
indexDir <- file.path(dir, "Dm.BDGP6.22.98_salmon-0.14.1")
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
              "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")
gtfPath <- file.path(dir,"Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
suppressPackageStartupMessages(library(tximeta))
makeLinkedTxome(indexDir=indexDir,
                source="Ensembl",
                organism="Drosophila melanogaster",
                release="98",
                genome="BDGP6.22",
                fasta=fastaFTP,
                gtf=gtfPath,
                write=FALSE)

## ------------------------------------------------------------------------
library(tximeta)
se <- tximeta(coldata)

## ----echo=FALSE----------------------------------------------------------
dir2 <- system.file("extdata", package="tximeta")
tab <- read.csv(file.path(dir2, "hashtable.csv"),
                stringsAsFactors=FALSE)
release.range <- function(tab, source, organism) {
  tab.idx <- tab$organism == organism & tab$source == source
  rels <- tab$release[tab.idx]
  if (organism == "Mus musculus" & source == "GENCODE") {
    paste0("M", range(as.numeric(sub("M","",rels))))
  } else if (source == "RefSeq") {
    paste0("p", range(as.numeric(sub(".*p","",rels))))
  } else {
    range(as.numeric(rels))
  }
}
dat <- data.frame(
  source=rep(c("GENCODE","Ensembl","RefSeq"),c(2,3,2)),
  organism=c("Homo sapiens","Mus musculus",
             "Drosophila melanogaster")[c(1:2,1:3,1:2)]
)
rng <- t(sapply(seq_len(nrow(dat)), function(i)
  release.range(tab, dat[i,1], dat[i,2])))
dat$releases <- paste0(rng[,1], "-", rng[,2])
knitr::kable(dat)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)

## ------------------------------------------------------------------------
assayNames(se)

## ------------------------------------------------------------------------
rowRanges(se)

## ------------------------------------------------------------------------
seqinfo(se)

## ------------------------------------------------------------------------
se.exons <- addExons(se)
rowRanges(se.exons)[[1]]

## ------------------------------------------------------------------------
gse <- summarizeToGene(se)
rowRanges(gse)

## ------------------------------------------------------------------------
library(org.Dm.eg.db)
gse <- addIds(gse, "REFSEQ", gene=TRUE)
mcols(gse)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))
# here there is a single sample so we use ~1.
# expect a warning that there is only a single sample...
suppressWarnings({dds <- DESeqDataSet(gse, ~1)})
# ... see DESeq2 vignette

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(edgeR))
cts <- assays(gse)[["counts"]]
normMat <- assays(gse)[["length"]]
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))
# ... see edgeR User's Guide for further steps

## ----eval=FALSE----------------------------------------------------------
#  y <- se # rename the object to 'y'
#  library(fishpond)
#  # if inferential replicates existed in the data,
#  # analysis would begin with:
#  #
#  # y <- scaleInfReps(y)
#  # ... see Swish vignette in the fishpond package

## ------------------------------------------------------------------------
gse <- summarizeToGene(se, countsFromAbundance="lengthScaledTPM")
library(limma)
y <- DGEList(assays(gse)[["counts"]])
# see limma User's Guide for further steps

## ------------------------------------------------------------------------
metadata(gse)$countsFromAbundance 

## ------------------------------------------------------------------------
names(metadata(se))
str(metadata(se)[["quantInfo"]])
str(metadata(se)[["txomeInfo"]])
metadata(se)[["tximetaInfo"]]
metadata(se)[["txdbInfo"]]

## ------------------------------------------------------------------------
file <- file.path(dir, "SRR1197474.plus", "quant.sf")
file.exists(file)
coldata <- data.frame(files=file, names="SRR1197474", sample="1",
                      stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
se <- tximeta(coldata)

## ------------------------------------------------------------------------
indexDir <- file.path(dir, "Dm.BDGP6.22.98.plus_salmon-0.14.1")

## ------------------------------------------------------------------------
fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz",
              "ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz",
              "extra_transcript.fa.gz")
#gtfFTP <- "ftp://path/to/custom/Drosophila_melanogaster.BDGP6.22.98.plus.gtf.gz"

## ------------------------------------------------------------------------
gtfPath <- file.path(dir,"Drosophila_melanogaster.BDGP6.22.98.plus.gtf.gz")

## ------------------------------------------------------------------------
tmp <- tempdir()
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
makeLinkedTxome(indexDir=indexDir,
                source="Ensembl", organism="Drosophila melanogaster",
                release="98", genome="BDGP6.22",
                fasta=fastaFTP, gtf=gtfPath,
                jsonFile=jsonFile)

## ------------------------------------------------------------------------
se <- tximeta(coldata)

## ------------------------------------------------------------------------
rowRanges(se)
seqinfo(se)

## ------------------------------------------------------------------------
library(BiocFileCache)
if (interactive()) {
  bfcloc <- getTximetaBFC()
} else {
  bfcloc <- tempdir()
}
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)

## ------------------------------------------------------------------------
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
loadLinkedTxome(jsonFile)

## ------------------------------------------------------------------------
se <- tximeta(coldata)

## ------------------------------------------------------------------------
if (interactive()) {
  bfcloc <- getTximetaBFC()
} else {
  bfcloc <- tempdir()
}
bfc <- BiocFileCache(bfcloc)
bfcinfo(bfc)
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
bfcinfo(bfc)

## ------------------------------------------------------------------------
library(devtools)
session_info()

