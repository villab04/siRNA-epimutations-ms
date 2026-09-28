library(DESeq2)
library(clusterProfiler)
library(AnnotationDbi)
library(GO.db)
library(enrichplot)

setwd("/home/bioc1877-local/Documents/MA_pombe/mRNApombe/raw_files/bam")
#Load counts from count features from mRNA seq
counts <- read.csv("counts.csv")
#Filter only data columns and eliminate Genid, Chr, Start, End, Strand and Length
data.counts <- counts[,c(7,8,9,10,11,12,13,14,15)]
#Other columns are: Genid, Chr, Start, End, Strand and Length
other.col <- counts[,c(1,2,3,4,5,6)]
sample <- colnames(data.counts)
condition <- as.factor(c("ancester", "control", "control", "control", "control", "ethanol",
               "ethanol", "ethanol", "ethanol"))
transfer <- as.factor(c("0", "10", "5", "10", "5", "10", "5", "10", "5"))
#col data is the metadata
coldata <- data.frame(sample = sample,
                        condition = condition,
                        transfer = transfer)
head(data.counts)
rownames(data.counts) <- paste(counts$Geneid)

dds_mRNA = DESeqDataSetFromMatrix(countData= data.counts, colData=coldata, design = ~condition)
dds_mRNA$condition <- factor(dds_mRNA$condition, levels = c("control","treated"))
dds_mRNA <- estimateSizeFactors(dds_mRNA)

normalized_counts <- counts(dds_mRNA, normalized=TRUE)
norm_counts <- as.data.frame(normalized_counts)
norm_counts <- rownames_to_column(norm_counts, var = "genID")

#save as csv the normalized counts for gene expression
write.csv(norm_counts, file = "/home/bioc1877-local/Documents/MA_pombe/mRNApombe/normalized_counts.csv",
            row.names = FALSE, quote = FALSE)
#READ AGAIN ALL THE FILES AND CREATE THE FILES TO ANALYSE WITH DESEQ
counts <- read.csv("counts.csv")
data.counts <- counts[,c(8,9,10,11,12,13,14,15)]
other.col <- counts[,c(1,2,3,4,5,6)]
sample <- colnames(data.counts)
condition <- as.factor(c("control", "control", "control", "control", "ethanol",
                         "ethanol", "ethanol", "ethanol"))
transfer <- as.factor(c("10", "5", "10", "5", "10", "5", "10", "5"))
coldata <- data.frame(sample = sample,
                      condition = condition,
                      transfer = transfer)

dds_mRNA = DESeqDataSetFromMatrix(countData= data.counts, colData=coldata, design = ~condition)
dds_mRNA <- estimateSizeFactors(dds_mRNA)
dse <- DESeq(dds_mRNA)
res = results(dse)
plotMA(res, ylim=c(-2,1))
data = cbind(other.col, data.frame(res))
data_clean <- na.omit(data)
gene_list <- data_clean$log2FoldChange
names(gene_list) <- data_clean$Geneid
gene_list <- sort(gene_list, decreasing = TRUE)

write.table(names(gene_list), file = "sig.genes_geneid", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################################################
#CREATE THE DATABASE TO CONTINUE WITH ENRICHMENT ANALYSIS
install.packages("/home/bioc1877-local/Documents/MA_pombe/sRNApombe/org.Spombe.eg.db")
library(org.Spombe.eg.db)

# Prepare TERM2GENE mapping: GO term <-> gene ID
go_df <- AnnotationDbi::select(org.Spombe.eg.db, 
                               keys=keys(org.Spombe.eg.db, keytype="GID"),
                               columns=c("GO", "GID"),
                               keytype="GID")
go_df <- na.omit(go_df)
go_df <- unique(go_df[, c("GO", "GID")])

# Prepare TERM2NAME (GO term <-> GO term description)
go_terms <- AnnotationDbi::select(GO.db, keys=unique(go_df$GO), columns=c("TERM"), keytype="GOID")

TERM2GENE <- go_df[, c("GO", "GID")]
colnames(TERM2GENE) <- c("term", "gene")

TERM2NAME <- go_terms[, c("GOID", "TERM")]
colnames(TERM2NAME) <- c("term", "name")
########################################################################################

gse <- GSEA(
  geneList=gene_list,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  minGSSize = 10,
)

gse.plot <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

ggsave("~/Documents/MA_pombe/sRNApombe/cluster_profiler_mRNA.svg", plot = gse.plot, 
       width = 17, height = 20, units = "cm", dpi = 300)

write.csv(gse, file = "gse.csv", row.names = FALSE)
