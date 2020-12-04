## ----biocstyle, echo = FALSE, results = "asis", message = FALSE------------
library(BiocStyle)
BiocStyle::markdown()

## ----load-libs, warning=FALSE, message=FALSE-------------------------------
library(EnsDb.Hsapiens.v86)

## Making a "short cut"
edb <- EnsDb.Hsapiens.v86
## print some informations for this package
edb

## For what organism was the database generated?
organism(edb)


## ----no-network, echo = FALSE, results = "hide"----------------------------
## Disable code chunks that require network connection - conditionally
## disable this on Windows only. This is to avoid TIMEOUT errors on the
## Bioconductor Windows build maching (issue #47).
use_network <- TRUE


## ----filters---------------------------------------------------------------
supportedFilters(edb)


## ----transcripts-----------------------------------------------------------
Tx <- transcripts(edb, filter = GeneNameFilter("BCL2L11"))

Tx

## as this is a GRanges object we can access e.g. the start coordinates with
head(start(Tx))

## or extract the biotype with
head(Tx$tx_biotype)


## ----transcripts-filter-expression-----------------------------------------
## Use a filter expression to perform the filtering.
transcripts(edb, filter = ~ gene_name == "ZBTB16")


## ----transcripts-filter, message = FALSE-----------------------------------
library(magrittr)

edb %>% filter(~ symbol == "BCL2" & tx_biotype != "protein_coding") %>%
    transcripts


## ----filter-Y--------------------------------------------------------------
edb_y <- addFilter(edb, SeqNameFilter("Y"))

## All subsequent filters on that EnsDb will only work on features encoded on
## chromosome Y
genes(edb_y)

## Get all lincRNAs on chromosome Y
genes(edb_y, filter = ~ gene_biotype == "lincRNA")


## ----list-columns----------------------------------------------------------
## list all database tables along with their columns
listTables(edb)

## list columns from a specific table
listColumns(edb, "tx")


## ----transcripts-example2--------------------------------------------------
Tx <- transcripts(edb,
		  columns = c(listColumns(edb , "tx"), "gene_name"),
		  filter = TxBiotypeFilter("nonsense_mediated_decay"),
		  return.type = "DataFrame")
nrow(Tx)
Tx


## ----cdsBy-----------------------------------------------------------------
yCds <- cdsBy(edb, filter = SeqNameFilter("Y"))
yCds


## ----genes-GRangesFilter---------------------------------------------------
## Define the filter
grf <- GRangesFilter(GRanges("11", ranges = IRanges(114129278, 114129328),
			     strand = "+"), type = "any")

## Query genes:
gn <- genes(edb, filter = grf)
gn

## Next we retrieve all transcripts for that gene so that we can plot them.
txs <- transcripts(edb, filter = GeneNameFilter(gn$gene_name))


## ----granges-zbtb16, message = FALSE, echo = FALSE-------------------------
plot(3, 3, pch = NA, xlim = c(start(gn), end(gn)), ylim = c(0, length(txs)),
     yaxt = "n", ylab = "")
## Highlight the GRangesFilter region
rect(xleft = start(grf), xright = end(grf), ybottom = 0, ytop = length(txs),
     col = "red", border = "red")
for(i in 1:length(txs)) {
    current <- txs[i]
    rect(xleft = start(current), xright = end(current), ybottom = i-0.975,
         ytop = i-0.125, border = "grey")
    text(start(current), y = i-0.5, pos = 4, cex = 0.75, labels = current$tx_id)
}


## ----transcripts-GRangesFilter---------------------------------------------
transcripts(edb, filter = grf)


## ----biotypes--------------------------------------------------------------
## Get all gene biotypes from the database. The GeneBiotypeFilter
## allows to filter on these values.
listGenebiotypes(edb)

## Get all transcript biotypes from the database.
listTxbiotypes(edb)


## ----genes-BCL2------------------------------------------------------------
## We're going to fetch all genes which names start with BCL.
BCLs <- genes(edb,
	      columns = c("gene_name", "entrezid", "gene_biotype"),
	      filter = GeneNameFilter("BCL", condition = "startsWith"),
	      return.type = "DataFrame")
nrow(BCLs)
BCLs


## ----example-AnnotationFilterList------------------------------------------
## determine the average length of snRNA, snoRNA and rRNA genes encoded on
## chromosomes X and Y.
mean(lengthOf(edb, of = "tx",
	      filter = AnnotationFilterList(
		  GeneBiotypeFilter(c("snRNA", "snoRNA", "rRNA")),
		  SeqNameFilter(c("X", "Y")))))

## determine the average length of protein coding genes encoded on the same
## chromosomes.
mean(lengthOf(edb, of = "tx",
	      filter = ~ gene_biotype == "protein_coding" &
		  seq_name %in% c("X", "Y")))


## ----example-first-two-exons-----------------------------------------------
## Extract all exons 1 and (if present) 2 for all genes encoded on the
## Y chromosome
exons(edb, columns = c("tx_id", "exon_idx"),
      filter = list(SeqNameFilter("Y"),
                    ExonRankFilter(3, condition = "<")))


## ----transcriptsBy-X-Y-----------------------------------------------------
TxByGns <- transcriptsBy(edb, by = "gene", filter = SeqNameFilter(c("X", "Y")))
TxByGns


## ----exonsBy-RNAseq, message = FALSE, eval = FALSE-------------------------
#  ## will just get exons for all genes on chromosomes 1 to 22, X and Y.
#  ## Note: want to get rid of the "LRG" genes!!!
#  EnsGenes <- exonsBy(edb, by = "gene", filter = AnnotationFilterList(
#  					  SeqNameFilter(c(1:22, "X", "Y")),
#  					  GeneIdFilter("ENSG", "startsWith")))
#  

## ----toSAF-RNAseq, message = FALSE, eval = FALSE---------------------------
#  ## Transforming the GRangesList into a data.frame in SAF format
#  EnsGenes.SAF <- toSAF(EnsGenes)
#  

## ----disjointExons, message = FALSE, eval = FALSE--------------------------
#  ## Create a GRanges of non-overlapping exon parts.
#  DJE <- disjointExons(edb, filter = AnnotationFilterList(
#  			      SeqNameFilter(c(1:22, "X", "Y")),
#  			      GeneIdFilter("ENSG%", "startsWith")))
#  

## ----transcript-sequence-AnnotationHub, message = FALSE, eval = FALSE------
#  library(EnsDb.Hsapiens.v86)
#  edb <- EnsDb.Hsapiens.v86
#  
#  ## Get the TwoBit with the genomic sequence matching the Ensembl version
#  ## using the AnnotationHub package.
#  dna <- ensembldb:::getGenomeTwoBitFile(edb)
#  
#  ## Get start/end coordinates of all genes.
#  genes <- genes(edb)
#  ## Subset to all genes that are encoded on chromosomes for which
#  ## we do have DNA sequence available.
#  genes <- genes[seqnames(genes) %in% seqnames(seqinfo(dna))]
#  
#  ## Get the gene sequences, i.e. the sequence including the sequence of
#  ## all of the gene's exons and introns.
#  geneSeqs <- getSeq(dna, genes)
#  

## ----transcript-sequence-extractTranscriptSeqs, message = FALSE, eval = FALSE----
#  ## get all exons of all transcripts encoded on chromosome Y
#  yTx <- exonsBy(edb, filter = SeqNameFilter("Y"))
#  
#  ## Retrieve the sequences for these transcripts from the TwoBitile.
#  library(GenomicFeatures)
#  yTxSeqs <- extractTranscriptSeqs(dna, yTx)
#  yTxSeqs
#  
#  ## Extract the sequences of all transcripts encoded on chromosome Y.
#  yTx <- extractTranscriptSeqs(dna, edb, filter = SeqNameFilter("Y"))
#  
#  ## Along these lines, we could use the method also to retrieve the coding sequence
#  ## of all transcripts on the Y chromosome.
#  cdsY <- cdsBy(edb, filter = SeqNameFilter("Y"))
#  extractTranscriptSeqs(dna, cdsY)
#  

## ----extractTranscriptSeqs-BSGenome, warning = FALSE, message = FALSE------
library(BSgenome.Hsapiens.NCBI.GRCh38)
bsg <- BSgenome.Hsapiens.NCBI.GRCh38

## Get the genome version
unique(genome(bsg))
unique(genome(edb))

## Extract the full transcript sequences.
yTxSeqs <- extractTranscriptSeqs(
  bsg, exonsBy(edb, "tx", filter = SeqNameFilter("Y")))

yTxSeqs

## Extract just the CDS
Test <- cdsBy(edb, "tx", filter = SeqNameFilter("Y"))
yTxCds <- extractTranscriptSeqs(
  bsg, cdsBy(edb, "tx", filter = SeqNameFilter("Y")))
yTxCds


## ----seqlevelsStyle, message = FALSE---------------------------------------
## Change the seqlevels style form Ensembl (default) to UCSC:
seqlevelsStyle(edb) <- "UCSC"

## Now we can use UCSC style seqnames in SeqNameFilters or GRangesFilter:
genesY <- genes(edb, filter = ~ seq_name == "chrY")
## The seqlevels of the returned GRanges are also in UCSC style
seqlevels(genesY)


## ----seqlevelsStyle-2, message = FALSE-------------------------------------
seqlevelsStyle(edb) <- "UCSC"

## Getting the default option:
getOption("ensembldb.seqnameNotFound")

## Listing all seqlevels in the database.
seqlevels(edb)[1:30]

## Setting the option to NA, thus, for each seqname for which no mapping is available,
## NA is returned.
options(ensembldb.seqnameNotFound=NA)
seqlevels(edb)[1:30]

## Resetting the option.
options(ensembldb.seqnameNotFound = "ORIGINAL")


## ----seqlevelsStyle-restore------------------------------------------------
seqlevelsStyle(edb) <- "Ensembl"


## ----gviz-plot, message=FALSE----------------------------------------------
## Loading the Gviz library
library(Gviz)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## Retrieving a Gviz compatible GRanges object with all genes
## encoded on chromosome Y.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "Y",
				start = 20400000, end = 21400000)
## Define a genome axis track
gat <- GenomeAxisTrack()

## We have to change the ucscChromosomeNames option to FALSE to enable Gviz usage
## with non-UCSC chromosome names.
options(ucscChromosomeNames = FALSE)

plotTracks(list(gat, GeneRegionTrack(gr)))

options(ucscChromosomeNames = TRUE)


## ----message=FALSE---------------------------------------------------------
seqlevelsStyle(edb) <- "UCSC"
## Retrieving the GRanges objects with seqnames corresponding to UCSC chromosome names.
gr <- getGeneRegionTrackForGviz(edb, chromosome = "chrY",
                                start = 20400000, end = 21400000)
seqnames(gr)
## Define a genome axis track
gat <- GenomeAxisTrack()
plotTracks(list(gat, GeneRegionTrack(gr)))


## ----gviz-separate-tracks, message=FALSE, warning=FALSE--------------------
protCod <- getGeneRegionTrackForGviz(edb, chromosome = "chrY",
				     start = 20400000, end = 21400000,
				     filter = GeneBiotypeFilter("protein_coding"))
lincs <- getGeneRegionTrackForGviz(edb, chromosome = "chrY",
				   start = 20400000, end = 21400000,
				   filter = GeneBiotypeFilter("lincRNA"))

plotTracks(list(gat, GeneRegionTrack(protCod, name = "protein coding"),
		GeneRegionTrack(lincs, name = "lincRNAs")), transcriptAnnotation = "symbol")

## At last we change the seqlevels style again to Ensembl
seqlevelsStyle <- "Ensembl"


## ----pplot-plot, message = FALSE, eval = FALSE-----------------------------
#  library(ggbio)
#  
#  ## Create a plot for all transcripts of the gene SKA2
#  autoplot(edb, ~ gene_name == "SKA2")
#  

## ----pplot-plot-2, message = FALSE, eval = FALSE---------------------------
#  ## Get the chromosomal region in which the gene is encoded
#  ska2 <- genes(edb, filter = ~ gene_name == "SKA2")
#  strand(ska2) <- "*"
#  autoplot(edb, GRangesFilter(ska2), names.expr = "gene_name")
#  

## ----AnnotationDbi, message = FALSE----------------------------------------
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

## List all available columns in the database.
columns(edb)

## Note that these do *not* correspond to the actual column names
## of the database that can be passed to methods like exons, genes,
## transcripts etc. These column names can be listed with the listColumns
## method.
listColumns(edb)

## List all of the supported key types.
keytypes(edb)

## Get all gene ids from the database.
gids <- keys(edb, keytype = "GENEID")
length(gids)

## Get all gene names for genes encoded on chromosome Y.
gnames <- keys(edb, keytype = "GENENAME", filter = SeqNameFilter("Y"))
head(gnames)


## ----select, message = FALSE, warning=FALSE--------------------------------
## Use the /standard/ way to fetch data.
select(edb, keys = c("BCL2", "BCL2L11"), keytype = "GENENAME",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))

## Use the filtering system of ensembldb
select(edb, keys = ~ gene_name %in% c("BCL2", "BCL2L11") &
		tx_biotype == "protein_coding",
       columns = c("GENEID", "GENENAME", "TXID", "TXBIOTYPE"))


## ----mapIds, message = FALSE-----------------------------------------------
## Use the default method, which just returns the first value for multi mappings.
mapIds(edb, keys = c("BCL2", "BCL2L11"), column = "TXID", keytype = "GENENAME")

## Alternatively, specify multiVals="list" to return all mappings.
mapIds(edb, keys = c("BCL2", "BCL2L11"), column = "TXID", keytype = "GENENAME",
       multiVals = "list")

## And, just like before, we can use filters to map only to protein coding transcripts.
mapIds(edb, keys = list(GeneNameFilter(c("BCL2", "BCL2L11")),
                        TxBiotypeFilter("protein_coding")), column = "TXID",
       multiVals = "list")


## ----AnnotationHub-query, message = FALSE, eval = use_network--------------
library(AnnotationHub)
## Load the annotation resource.
ah <- AnnotationHub()

## Query for all available EnsDb databases
query(ah, "EnsDb")


## ----AnnotationHub-query-2, message = FALSE, eval = use_network------------
ahDb <- query(ah, pattern = c("Xiphophorus Maculatus", "EnsDb", 87))
## What have we got
ahDb


## ----AnnotationHub-fetch, message = FALSE, eval = FALSE--------------------
#  ahEdb <- ahDb[[1]]
#  
#  ## retriebe all genes
#  gns <- genes(ahEdb)
#  

## ----edb-from-ensembl, message = FALSE, eval = FALSE-----------------------
#  library(ensembldb)
#  
#  ## get all human gene/transcript/exon annotations from Ensembl (75)
#  ## the resulting tables will be stored by default to the current working
#  ## directory
#  fetchTablesFromEnsembl(75, species = "human")
#  
#  ## These tables can then be processed to generate a SQLite database
#  ## containing the annotations (again, the function assumes the required
#  ## txt files to be present in the current working directory)
#  DBFile <- makeEnsemblSQLiteFromTables()
#  
#  ## and finally we can generate the package
#  makeEnsembldbPackage(ensdb = DBFile, version = "0.99.12",
#                       maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>",
#                       author = "J Rainer")
#  

## ----gtf-gff-edb, message = FALSE, eval = FALSE----------------------------
#  ## Load the AnnotationHub data.
#  library(AnnotationHub)
#  ah <- AnnotationHub()
#  
#  ## Query all available files for Ensembl release 77 for
#  ## Mus musculus.
#  query(ah, c("Mus musculus", "release-77"))
#  
#  ## Get the resource for the gtf file with the gene/transcript definitions.
#  Gtf <- ah["AH28822"]
#  ## Create a EnsDb database file from this.
#  DbFile <- ensDbFromAH(Gtf)
#  ## We can either generate a database package, or directly load the data
#  edb <- EnsDb(DbFile)
#  
#  
#  ## Identify and get the TwoBit object with the genomic DNA sequence matching
#  ## the EnsDb annotation.
#  Dna <- getGenomeTwoBitFile(edb)
#  ## We next retrieve the sequence of all exons on chromosome Y.
#  exons <- exons(edb, filter = SeqNameFilter("Y"))
#  exonSeq <- getSeq(Dna, exons)
#  

## ----EnsDb-from-Y-GRanges, message = FALSE, eval = use_network-------------
## Generate a sqlite database from a GRanges object specifying
## genes encoded on chromosome Y
load(system.file("YGRanges.RData", package = "ensembldb"))
Y

## Create the EnsDb database file
DB <- ensDbFromGRanges(Y, path = tempdir(), version = 75,
		       organism = "Homo_sapiens")

## Load the database
edb <- EnsDb(DB)
edb


## ----EnsDb-from-GTF, message = FALSE, eval = FALSE-------------------------
#  library(ensembldb)
#  
#  ## the GTF file can be downloaded from
#  ## ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/
#  gtffile <- "Homo_sapiens.GRCh37.75.gtf.gz"
#  ## generate the SQLite database file
#  DB <- ensDbFromGtf(gtf = gtffile)
#  
#  ## load the DB file directly
#  EDB <- EnsDb(DB)
#  
#  ## alternatively, build the annotation package
#  ## and finally we can generate the package
#  makeEnsembldbPackage(ensdb = DB, version = "0.99.12",
#                       maintainer = "Johannes Rainer <johannes.rainer@eurac.edu>",
#                       author = "J Rainer")
#  

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

