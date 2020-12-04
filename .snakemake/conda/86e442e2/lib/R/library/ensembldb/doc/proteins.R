## ----biocstyle, echo = FALSE, results = "asis", message = FALSE------------
library(BiocStyle)
library(ensembldb)
BiocStyle::markdown() 

## ----doeval, echo = FALSE, results = "hide"--------------------------------
## Globally switch off execution of code chunks
evalMe <- TRUE
haveProt <- FALSE
## evalMe <- .Platform$OS.type == "unix"
 

## ----loadlib, message = FALSE, eval = evalMe-------------------------------
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
## Evaluate whether we have protein annotation available
hasProteinData(edb)
 

## ----restrict9, message = FALSE, echo = FALSE------------------------------
## silently subsetting to chromosome 11
edb <- filter(edb, filter = ~ seq_name == "11") 

## ----listCols, message = FALSE, eval = evalMe------------------------------
listTables(edb)
 

## ----haveprot, echo = FALSE, results = "hide", eval = evalMe---------------
## Use this to conditionally disable eval on following chunks
haveProt <- hasProteinData(edb) & evalMe
 

## ----a_transcripts, eval = haveProt----------------------------------------
## Get also protein information for ZBTB16 transcripts
txs <- transcripts(edb, filter = GeneNameFilter("ZBTB16"),
		   columns = c("protein_id", "uniprot_id", "tx_biotype"))
txs
 

## ----a_transcripts_coding_noncoding, eval = haveProt-----------------------
## Subset to transcripts with tx_biotype other than protein_coding.
txs[txs$tx_biotype != "protein_coding", c("uniprot_id", "tx_biotype",
					  "protein_id")]
 

## ----a_transcripts_coding, eval = haveProt---------------------------------
## List the protein IDs and uniprot IDs for the coding transcripts
mcols(txs[txs$tx_biotype == "protein_coding",
	  c("tx_id", "protein_id", "uniprot_id")])
 

## ----a_transcripts_coding_up, eval = haveProt------------------------------
## List all uniprot mapping types in the database.
listUniprotMappingTypes(edb)

## Get all protein_coding transcripts of ZBTB16 along with their protein_id
## and Uniprot IDs, restricting to protein_id to uniprot_id mappings based
## on "DIRECT" mapping methods.
txs <- transcripts(edb, filter = list(GeneNameFilter("ZBTB16"),
				      UniprotMappingTypeFilter("DIRECT")),
		   columns = c("protein_id", "uniprot_id", "uniprot_db"))
mcols(txs)
 

## ----a_genes_protdomid_filter, eval = haveProt-----------------------------
## Get all genes encoded on chromosome 11 which protein contains 
## a certain protein domain.
gns <- genes(edb, filter = ~ prot_dom_id == "PS50097" & seq_name == "11")
length(gns)

sort(gns$gene_name)
 

## ----a_2_annotationdbi, message = FALSE, eval = haveProt-------------------
## Show all columns that are provided by the database
columns(edb)

## Show all key types/filters that are supported
keytypes(edb)
 

## ----a_2_select, message = FALSE, eval = haveProt--------------------------
select(edb, keys = "ZBTB16", keytype = "GENENAME",
       columns = "UNIPROTID")
 

## ----a_2_select_nmd, message = FALSE, eval = haveProt----------------------
## Call select, this time providing a GeneNameFilter.
select(edb, keys = GeneNameFilter("ZBTB16"),
       columns = c("TXBIOTYPE", "UNIPROTID", "PROTEINID"))
 

## ----b_proteins, message = FALSE, eval = haveProt--------------------------
## Get all proteins and return them as an AAStringSet
prts <- proteins(edb, filter = GeneNameFilter("ZBTB16"),
		 return.type = "AAStringSet")
prts
 

## ----b_proteins_mcols, message = FALSE, eval = haveProt--------------------
mcols(prts)
 

## ----b_proteins_prot_doms, message = FALSE, eval = haveProt----------------
## Get also protein domain annotations in addition to the protein annotations.
pd <- proteins(edb, filter = GeneNameFilter("ZBTB16"),
	       columns = c("tx_id", listColumns(edb, "protein_domain")),
	       return.type = "AAStringSet")
pd
 

## ----b_proteins_prot_doms_2, message = FALSE, eval = haveProt--------------
## The number of protein domains per protein:
table(names(pd))

## The mcols
mcols(pd)
 

## ----sessionInfo-----------------------------------------------------------
sessionInfo() 

