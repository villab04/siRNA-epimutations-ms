library(tidyverse)
library(tidyr)
library(dplyr)
library(DESeq2)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ggVennDiagram)

setwd("/home/bioc1877-local/Documents/MA_pombe/sRNApombe")
list.files()

bed.file <- read.table("ANC.clean.fa.bed", sep="\t")
head(bed.file)
filen <- "ANC.clean.fa.bed"


#This code filters mapped small RNA reads to keep only those likely 22–26 nucleotides 
#long, then writes a cleaner BED file with their genomic location, count, score, strand,
#and name. It assumes the read name contains the count and length information.

#TRying it in one sample before making a loop.
Input <- cbind(bed.file, read.table(text=bed.file[,4], sep="_")) #separate the lenght
Input1 <- cbind(Input,  read.table(text=Input[,7], sep="-")) # separate depth
filtered <- Input1[Input1[, 9] > 21 & Input1[, 9] < 27, ] #filtrate everything that is in between 21 and 27
filtered_input <- filtered[, c(1:3,11,5,6,4)] # write out the final columns
colnames(filtered_input) <- c("Chrm", "Start", "End", "Count", "Score", "Strand", "Name")
outfile <- sub("\\.clean\\.fa\\.bed$", ".filtered.bed", filen)
write.table(filtered_input, outfile, sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE) # print the final "filtered" file



FilterFile <- function(filein) {
  bed.file <- read.table(filein, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  Input <- cbind(bed.file, read.table(text = bed.file[, 4], sep = "_"))
  Input1 <- cbind(Input, read.table(text = Input[, 7], sep = "-"))
  filtered <- Input1[Input1[, 9] > 21 & Input1[, 9] < 27, ]
  filtered_input <- filtered[, c(1:3, 11, 5, 6, 4)]
  #colnames(filtered_input) <- c("Chrm", "Start", "End", "Count", "Score", "Strand", "Name")
  outfile <- sub("\\.clean\\.fa\\.bed$", ".filtered.bed", filein)
  write.table(filtered_input, outfile, sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  message("Saved: ", outfile)
}


files <- list.files(pattern = "clean\\.fa\\.bed$")
for (f in files) {
  FilterFile(f)
  }


################################################################################
#Create the matrix with the count files for window data
files <- list.files(pattern="*\\.windows\\.bed")
files <- files[!grepl("L2T9", files)] #remove L2T9 line because its and outlier
#files <- files[!grepl("_R2", files)] #remove duplicates only when necesary
#create a list of the data frames
dfs <- lapply(files, function(f) {
  df <- read.table(f, header=FALSE, sep="\t",
                   col.names=c("chrom","start","end","count"))
  sample <- sub("\\.counts\\.bed","",f)
  #creates a list
  df$sample <- sample
  return(df)
})

# reshape to matrix
df_all <- bind_rows(dfs)
df_all$count[df_all$count == "."] <- 0
df_all$count <- as.numeric(df_all$count)

matrix <- reshape(
  df_all,
  idvar = c("chrom", "start", "end"),
  timevar = "sample",
  direction = "wide"
)

#Prepare the data for DEseq
mat <- as.matrix(matrix[ , -(1:3)])

rownames(mat) <- paste(matrix$chrom, matrix$start, matrix$end, sep=":")

coldata <- data.frame(
  row.names = colnames(mat),
  condition = c("ANC", 
                "control", "control", "control", "control", "control", "control",
                "control", "control", "control", "control", "control", "control", "control", 
                "control", "control", "control", "control", "control", "control", "control", 
                "control", "control", "control", "control", 
                "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%",
                "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%"))


#convert to deseq object
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData   = coldata,
                              design    = ~1)

###################################
#Normalize the data
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
norm_counts <- as.data.frame(normalized_counts)
norm_counts$coord <- rownames(norm_counts)

rownames(norm_counts) <- NULL
  norm_counts <- norm_counts %>%
  separate(coord, into = c("chrom", "start", "end"), sep = ":", convert = TRUE)

  norm_counts <- norm_counts %>%
    dplyr::select(chrom, start, end, everything())
  
write.csv(norm_counts, "normalized_counts_windows.csv", row.names = FALSE)
###############################################################################
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)

# Percent variance explained
percentVar <- round(100 * attr(pcaData, "percentVar"))

#plot the results
ggplot(pcaDATA, aes(PC1, PC2, color=condition, label =name)) +
  geom_point(size=3) +
  geom_text(vjust = -0.7, size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


############################################################################
############################################################################
#Create the matrix with the count files for CDS count data
cds.files <- list.files(pattern="*.cds.bed_sum$")

#read the files
cds.dfs <- lapply(cds.files, function(f) {
  cds.df <- read.table(f, header=FALSE, sep="\t",
                   col.names=c("chrom","start","end", "name", "count"))
  sample <- sub("*.cds.bed_sum$","",f) #remove the pattern of the filename
  #creates a list
  cds.df$sample <- sample
  return(cds.df)
})

# reshape to matrix
#stacks all those data frames into one big table
cds.df_all <- bind_rows(cds.dfs)

cds.df_all <- cds.df_all %>%
  mutate(count = ifelse(count == ".", 0, as.numeric(count)))%>%
  group_by(chrom, start, end, name, sample) %>%
  summarise(count = sum(count), .groups = "drop")

dt <- as.data.table(cds.df_all)

#dcast need to be used because the data set is larger
cds.matrix <- dcast(
  dt,
  chrom + start + end + name ~ sample,
  value.var = "count",
  fun.aggregate = sum
)


#Prepare the data for DEseq
cds.mat <- as.matrix(cds.matrix[ , -(1:4)])

rownames(cds.mat) <- paste(cds.matrix$chrom, cds.matrix$start,
                           cds.matrix$end, cds.matrix$name, sep="_")

cds.coldata <- data.frame(
  row.names = colnames(cds.mat),
  condition = c("ANC", 
                "control", "control", "control", "control", "control", "control", "control",
                "control", "control", "control", "control", "control", "control", "control", 
                "control", "control", "control", "control", "control", "control", "control", 
                "control", "control", "control", "control",
                "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%",
                "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%", "4%"),
  
  generation = c("ANC", "1", "10", "2", "2", "3", "4", "5", "6", "7", "7", "8", "8", "9",
                 "1", "10", "2", "3", "3", "4", "5", "6", "7", "7", "8", "9", 
                 "1", "10", "2", "3", "4", "4", "5", "6", "6", "7", "8", "9", "1", "10", "2",
                 "2", "3", "4", "4", "5", "6", "6", "7", "8", "9"))

#convert to deseq object
cds.dds <- DESeqDataSetFromMatrix(countData = cds.mat,
                              colData   = cds.coldata,
                              design    = ~1)
##############################################################################
#Normalize the data
cds.dds <- estimateSizeFactors(cds.dds)
sizeFactors(cds.dds)
cds.normalized_counts <- counts(cds.dds, normalized=TRUE)
cds.norm_counts <- as.data.frame(cds.normalized_counts)
cds.norm_counts$coord <- rownames(cds.norm_counts)

rownames(cds.norm_counts) <- NULL
cds.norm_counts <- cds.norm_counts %>%
  # Use regex to safely extract 4 parts
  extract(coord, into = c("chrom", "start", "end", "name"), 
          regex = "^([^_]+)_([^_]+)_([^_]+)_(.+)$", convert = TRUE) %>%
  select(chrom, start, end, name, everything())

write.csv(cds.norm_counts, "cds.normalized_counts_windows.csv", row.names = FALSE)

dse.cds <- DESeq(cds.dds)
res.cds = results(dse.cds)
plotMA(res.cds, ylim=c(-2,1))
data_clean.cds <- na.omit(res.cds)
gene_list <- data_clean.cds$log2FoldChange
names(gene_list) <- data_clean.cds$Geneid
gene_list <- sort(gene_list, decreasing = TRUE)

###################################################################################

#transform the data (variance stabilized values)
cds.vsd <- varianceStabilizingTransformation(cds.dds, blind=TRUE, fitType = "local")   
cds.pcaDATA <- plotPCA(cds.vsd, intgroup= c("condition", "generation"), returnData=TRUE)
cds.percentVar <- round(100 * attr(cds.pcaDATA, "percentVar"))

#plot the results
ggplot(cds.pcaDATA, aes(PC1, PC2, color=generation, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",cds.percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",cds.percentVar[2],"% variance")) + 
  coord_fixed()

ggplot(cds.pcaDATA, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


###########################################################################

#Create the matrix with the count files for filtered data
all.files <- list.files(pattern = "L.*\\.filtered\\.bed$")
filtered.files <- all.files[!grepl("_R2", all.files)]

filtered.dfs <- lapply(filtered.files, function(f) {
  filtered.df <- read.table(f, header=FALSE, sep="\t",
                            col.names=c("chrom","start","end", "count", 
                                        "score", "strand", "code"))
  # Replace '.' with 0 in the 'count' column
  filtered.df$count[filtered.df$count == "."] <- 0
  filtered.df$count <- as.integer(filtered.df$count)
  sample <- sub("\\.cov\\.bed","",f)
  #creates a list
  filtered.df$sample <- sample
  return(filtered.df)
})

# reshape to matrix
filtered.df_all <- bind_rows(filtered.dfs)

filtered.df_all_summary <- filtered.df_all %>%
  group_by(chrom, start, end, sample) %>%
  summarise(count = sum(count), .groups = "drop")

filtered.matrix <- pivot_wider(filtered.df_all_summary,
                          id_cols=c(chrom,start,end),
                          names_from=sample,
                          values_from=count,
                          values_fill = 0)
#--------------------------------------------------------------
#create a database for use in the genomic_distributions script-
sirna_df <- as.data.frame(filtered.matrix) 
df_control <- sirna_df[, c("chrom", "start", "end", 
                   colnames(sirna_df)[grepl("L1|L2", colnames(sirna_df))])]

df_ethanol <- sirna_df[, c("chrom", "start", "end", 
                           colnames(sirna_df)[grepl("L6|L7", colnames(sirna_df))])]

df_control$total_count <- rowSums(df_control[, 4:ncol(df_control)])
df_control <- df_control[, c("chrom", "start", "end", "total_count")]

df_ethanol$total_count <- rowSums(df_ethanol[, 4:ncol(df_ethanol)])
df_ethanol <- df_ethanol[, c("chrom", "start", "end", "total_count")]

sirna_df$total_count <- rowSums(sirna_df[, 4:ncol(sirna_df)])
df_sirnall <- sirna_df[, c("chrom", "start", "end", "total_count")]

sirna_windows_matrix <- as.data.frame(matrix) 
sirna_windows_matrix$total_count <- rowSums(sirna_windows_matrix[, 4:ncol(sirna_windows_matrix)])
df_sirn_windowsall <- sirna_windows_matrix[, c("chrom", "start", "end", "total_count")]
df_sirn_windowsall <- df_sirn_windowsall %>% filter(total_count > 0)
#-------------------------------------------------------------

#Prepare the data for DEseq
filtered.mat <- as.matrix(filtered.matrix[ , -(1:4)])

rownames(filtered.mat) <- paste(filtered.matrix$chrom, filtered.matrix$start,
                                filtered.matrix$end, sep="_")

filtered.coldata <- data.frame(
  row.names = colnames(filtered.mat),
  condition = c("control", "ethanol", "ethanol", "ethanol", 
                 "ethanol", "ethanol", "control", "control",
                "control", "control", "control", "control",
                "ethanol", "control", "control",
                "ethanol", "control", "control", "control",
                "ethanol", "ethanol", "ethanol", "control",
                "control", "ethanol", "control", "ethanol",
                "ethanol", "control", "ethanol", "ethanol",
                "ethanol", "ethanol", "ethanol", "ethanol",
                "control", "ethanol", "ethanol", "ethanol"),
  generation = c("1", "9", "8", "1", "2", "4", "10", "2",
                 "8", "8", "4", "5", "2", "5", "4",
                 "10", "7", "2", "6", "3", "5", "9", "6",
                 "1", "3", "7", "6", "6", "3", "4", "7",
                 "9", "10", "5", "7", "9", "8", "1", "3"))

#convert to deseq object
filtered.dds <- DESeqDataSetFromMatrix(countData = filtered.mat,
                                  colData   = filtered.coldata,
                                  design    = ~condition)

dds <- DESeq(filtered.dds)
res <- results(dds)
plotMA(res, ylim=c(-10,10))
res <- res[!is.na(res$padj), ]
sig_res <- res[res$padj < 0.05, ]

genes_df <- data.frame(
  region = sub("^(.*)_(\\d+)_\\d+$", "\\1", rownames(df.sig_res)),  # prefix (e.g. mating_type_region)
  start  = as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", rownames(df.sig_res))), # start number
  end    = as.numeric(sub(".*_(\\d+)$", "\\1", rownames(df.sig_res))),      # end number
  stringsAsFactors = FALSE
)

write.table(genes_df, file="constant_sRNA_deseq_results.bed", sep ="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Do the enrichment analysis of constat sRNA difference bewteen treatments
#Import GFF3 file to math the coordinates with gene names
gff <- import("~/Documents/MA_pombe/reference/Schizosaccharomyces_pombe_all_chromosomes.gff3")

#select only the "type=gene"
genes_gr <- gff[gff$type == "gene"]

res <- as.data.frame(res)

res$chr <- sub("^(.*)_(\\d+)_\\d+$", "\\1", rownames(res))
res$start <- as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", rownames(res))) # start number
res$end <- as.numeric(sub(".*_(\\d+)$", "\\1", rownames(res)))      # end number

#convert the data frame in Granges
df_gr <- GRanges(
  seqnames = res$chr,
  ranges = IRanges(start = res$start, end = res$end))
#Find the hits of the coordinates with genes

hits <- findOverlaps(df_gr, genes_gr)

result <- data.frame(
res[queryHits(hits), ],
  gene_id = mcols(genes_gr)$ID[subjectHits(hits)])
#prepare the genelist for cluster profiler

universe_genes <- unique(result$gene_id)
sig_res <- unique(result$gene_id[result$padj < 0.05])

#run GSEA with cluster profiler, The preparation of CLuster profiler 
#in in the clusterProfiler_mRNA.R script

ora <- enricher(
  gene = sig_res,
  universe = universe_genes,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10)

#0 enriched terms found

###Create Venn diagram bewteen mRNA and siRNA constant differences

mRNA_genes <- read.table("~/Documents/MA_pombe/mRNApombe/raw_files/bam/sig.genes_geneid")
sRNA_genes <- read.table("~/Documents/MA_pombe/sRNApombe/constant_sRNA_geneid")

mRNA_genes <- unique(mRNA_genes[,1])
sRNA_genes <- unique(sRNA_genes[,1])

x <- list(
  "Significant mRNA" = mRNA_genes,
  "Constant sRNA" = sRNA_genes)


Venn.diag <- ggVennDiagram(x[1:2], label_alpha = 0, label = "count", alpha = 0.5) +
  theme(legend.position = "none")

N <- 12456  # replace with your actual gene universe
mat1 <- matrix(c(27, 6074,
                115, N - 6216),
              nrow = 2, byrow = TRUE)

fisher.test(mat1, alternative = "greater")

ggsave("~/Documents/MA_pombe/sRNApombe/Venn_diag.svg", plot = Venn.diag, 
       width = 10, height = 10, units = "cm", dpi = 300)



constant_changes <- data.frame(
  type = c("protein coding", "tRNA gene", "lncRNA", "centromeric repeats", "rRNA", "snoRNA"),
  all.genesvalue = c(5, 8, 71, 29, 23, 2),
  intersect.genesvalue = c(4, 0, 14, 0, 9, 0)
)

constant_changes <- constant_changes %>%
  arrange(desc(all.genesvalue)) %>%
  mutate(type = factor(type, levels = type))

constant_changes_long <- constant_changes |>
  pivot_longer(
    cols = c(all.genesvalue, intersect.genesvalue),
    names_to = "group",
    values_to = "value"
  )|>
  mutate(
    group = recode(
      group,
      all.genesvalue = "Genes with different siRNA counts",
      intersect.genesvalue = "Overlap of siRNA DEG with DEG"
    )
  )

siRNAgenes.plot <- ggplot(constant_changes_long, aes(x = type, y = value, fill = group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c(
      "Genes with different siRNA counts" = "#6A5ACD",
      "Overlap of siRNA DEG with DEG" = "#20B2AA"
    ), name = NULL) +
  labs( y = "Number of genes")+
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.size = unit(0.5, "cm"))

ggsave("~/Documents/MA_pombe/sRNApombe/siRNAgenes.plot.svg", plot = siRNAgenes.plot, width = 25, height = 30, units = "cm", dpi = 300)

pie.plotintersect <- ggplot(constant_changes, aes(x = "", y = intersect.genesvalue, fill = type)) +
  geom_col(color = NA) +
  geom_text(
    aes(label = intersect.genesvalue),
    position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill = NULL)+
  theme_minimal()

ggsave("~/Documents/MA_pombe/sRNApombe/pie.plotintersect.svg", plot = pie.plotintersect, width = 20, height = 20, units = "cm", dpi = 300)
