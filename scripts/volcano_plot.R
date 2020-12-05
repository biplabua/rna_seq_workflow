library(ggplot2)
library(ggrepel)
de <- read.csv(snakemake@input[["log2fc"]], row.names = 1)
de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"
de <- de[order(de$padj), ] 
de$genelabels <- ""
de$genelabels[1:10] <- rownames(de)[1:10]
png(snakemake@output[["volcano_plot"]])
ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=genelabels)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()
dev.off() 
