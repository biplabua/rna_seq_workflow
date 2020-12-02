library(DESeq2)
library(tximport)
library(readr)
samples <- read.csv(snakemake@input[["meta_data"]], row.names = 1)
files <- file.path(snakemake@input[["file_path"]], "quant.sf")
names(files) <- samples$id
txi <- tximport(files = files, type = "salmon", txOut = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ treatment)
ds <- DESeq(ddsTxi, test="Wald")
df <- results(ds, contrast = c("treatment", "control", "treated"))
write.csv(df, file = snakemake@output[["out_file"]])