library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
library(dplyr)
library(AnnotationForge)
library(BiocManager)
library(GO.db)
library(enrichplot)
library(biomaRt)

setwd("/home/bioc1877-local/Documents/MA_pombe/sRNApombe")
ANC <- readLines("ANC_geneid")
#CONSTANT <- readLines("constant_sRNA_geneid")
#ALL.EPI <- readLines("all_epi.bed_geneid")#this were created with the old code, probably wrong
#CONTROL <- readLines("all_control.bed_geneid")#this were created with the old code, probably wrong
#ETHANOL <- readLines("all_ethanol.bed_geneid")#this were created with the old code, probably wrong

#Load the geneID from the new code(august 2026)
load("epimutations_GeneID.RData")
ALL.EPI <- c(control_geneID, ethanol_geneID)

fungi_mart <- useMart("fungi_mart", host = "https://fungi.ensembl.org")
datasets <- listDatasets(fungi_mart)
View(datasets)
ensembl <- useDataset("spombe_eg_gene", mart = fungi_mart)

gene_df <- getBM(attributes = c("ensembl_gene_id", 
                                "external_gene_name", 
                                "description", 
                                "go_id", 
                                "namespace_1003"),
                 mart = ensembl)

gene_df <- gene_df[gene_df$go_id != "", ]
head(gene_df)

gene_info <- gene_df %>%
  dplyr::select(
    GID = ensembl_gene_id,
    SYMBOL = external_gene_name,
    GENENAME = description
  ) %>%
  dplyr::distinct()


go <- gene_df %>%
  dplyr::select(
    GID = ensembl_gene_id,
    GO = go_id,
    EVIDENCE = namespace_1003
  ) %>%
  dplyr::mutate(EVIDENCE = ifelse(is.na(EVIDENCE) | EVIDENCE == "", "IEA", EVIDENCE)) %>%
  dplyr::distinct()


makeOrgPackage(
  gene_info = gene_info,
  go = go,
  version = "1.0",
  maintainer = "Mariana Villalba <mariana.villalbadelapena@bioch.ox.ac.uk>",
  author = "Mariana Villalba",
  outputDir = ".",
  tax_id = "4896",  # NCBI taxonomy ID for S. pombe
  genus = "Schizosaccharomyces",
  species = "pombe"
)

devtools::install("./org.Spombe.eg.db", repos = NULL, type = "source")
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

gene_set_sizes <- table(TERM2GENE$term)
summary(gene_set_sizes)
all.genome <- TERM2GENE$gene
# Now run enrichment with the gene list
ANC_enrich<- enricher(
  gene = ANC,
  universe = all.genome,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 1
)

ANC.plot <- barplot(ANC_enrich) + theme(
  axis.text.y = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  legend.position = "none")

ggsave(filename = "ORA_anc_enrich.svg", plot = ANC.plot,
       width = 17, height = 23, units = "cm")

allepi_enrich <- enricher(
  gene = ALL.EPI,
  universe = all.genome,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 1
)

barplot(allepi_enrich)

p <- barplot(allepi_enrich, color = "p.adjust")
enrich.plot <- p +
  scale_fill_distiller(
    palette = "Spectral",
    direction = 1)+
  theme(
      axis.title = element_text(size = 20),     
      axis.text = element_text(size = 20),        
      legend.title = element_text(size = 18),     
      legend.text = element_text(size = 15),      
      text = element_text(size = 20)             
    )

ggsave(filename = "enrich_plot.svg", plot = enrich.plot,
       width = 15, height = 15, units = "cm")

all_epi_enrich_df <- as.data.frame(allepi_enrich)
write.csv(all_epi_enrich_df, file = "all_epi_enrich_results.csv", row.names = FALSE)

control_enrich <- enricher(
  gene = control_geneID,
  universe = all.genome,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 1
)

control.ORA.plot <- barplot(control_enrich,  x = "Count")+
    theme(
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20))

ggsave(filename = "ORA_control_enrich.svg", plot = control.ORA.plot,
       width = 17, height = 23, units = "cm")

ethanol_enrich <- enricher(
  gene = ethanol_geneID,
  universe = all.genome,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 1
)

ethanol.ora.plot <- barplot(ethanol_enrich)+
  theme(
  axis.text.y = element_text(size = 18),
  axis.text.x = element_text(size = 18),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20))

ggsave(filename = "ORA_ethanol_enrich.svg", plot = ethanol.ora.plot,
       width = 17, height = 23, units = "cm")

control_geneID <- unique(unlist(strsplit(control_enrich$geneID, "/")))
ethanol_geneID <- unique(unlist(strsplit(ethanol_enrich$geneID, "/")))
