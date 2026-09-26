#Genomic Distributions and regioneR
library(regioneR)
library("GenomicDistributions")
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)
library(ggplot2)

setwd("~/Documents/MA_pombe/sRNApombe")
fastaSource = "/home/bioc1877-local/Documents/MA_pombe/reference/Schizosaccharomyces_pombe_all_chromosomes.fa"
gtfSource = "/home/bioc1877-local/Documents/MA_pombe/reference/Schizosaccharomyces_pombe_all_chromosomes.gff3"

#query = rtracklayer::import("all_epi.bed")
#query = rtracklayer::import("all_epi_intersect.bed")
#query_control = rtracklayer::import("all_control.bed")
#query_ethanol = rtracklayer::import("all_ethanol.bed")
transposons = rtracklayer::import("~/Documents/MA_pombe/reference/transposons.bed")
#load the new coordinates
load("epimutations_coord.RData")
control.epicoord <- control.epicoord %>% filter(chrom != "mating_type_region", chrom != "mitochondrial")
query_control <- makeGRangesFromDataFrame(
  control.epicoord,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE)

ethanol.epicoord <- ethanol.epicoord %>% filter(chrom != "mating_type_region", chrom != "mitochondrial")
query_ethanol <- makeGRangesFromDataFrame(
  control.epicoord,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE)

all.epicoord <- all.epicoord %>% filter(chrom != "mating_type_region", chrom != "mitochondrial")
query_all <- makeGRangesFromDataFrame(
  all.epicoord,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = FALSE)

#These queries comes from the filter.bed script and are the siRNA count only filtered by size
query_srna_all = GRanges(df_sirn_windowsall)
query_srna_all <- query_srna_all[seqnames(query_srna_all) %in% c("I", "II", "III")]
# Get Chromosome sites
SPChromSizes = getChromSizesFromFasta(source = fastaSource)
bins  = getGenomeBins(SPChromSizes)
#Get the TSS
SPTSSs = getTssFromGTF(source=gtfSource, filterProteinCoding = FALSE)
#Get the features for gene models

features = c("gene", "exon", "three_prime_UTR", "five_prime_UTR", "mRNA", "promoter", 
             "repeat_region", "regional_centromere", "rRNA", "CDS")
SPGeneModels = getGeneModelsFromGTF(source=gtfSource, features=features, convertEnsemblUCSC=TRUE,
                                    filterProteinCoding = FALSE)

for(feature in names(SPGeneModels)) {
  # Fix seqlevels
  seqlevels(SPGeneModels[[feature]]) <- sub("^chr", "", seqlevels(SPGeneModels[[feature]]))
  # Fix seqnames
  seqnames(SPGeneModels[[feature]]) <- sub("^chr", "", as.character(seqnames(SPGeneModels[[feature]])))
}


partitionList = genomePartitionList(SPGeneModels$gene, 
                                    SPGeneModels$exon,
                                    SPGeneModels$three_prime_utr, 
                                    SPGeneModels$five_prime_utr)

partitionList[["CDS"]] <- SPGeneModels$CDS
partitionList[["repeat_region"]] <- SPGeneModels$repeat_region
partitionList[["rRNA"]] <- SPGeneModels$rRNA
partitionList[["regional_centromere"]] <- SPGeneModels$regional_centromere

#Chromosome distribution plots
#control and ethanol separated
queryList = GRangesList(control=query_control, ethnaol=query_ethanol)
x2 = calcChromBins(queryList, bins)
plotChromBins(x2)
#sRNAs raw 
gp_raw = calcPartitions(query_srna_all, partitionList, bpProportion=TRUE)
plotPartitions(gp_raw)
#############################
#control and ethanol separated
TSSdist2 = calcFeatureDist(queryList, SPTSSs)
plotFeatureDist(TSSdist2, featureName="TSS", , size = 1000, nbins = 20, 
                tile = TRUE, labelOrder = "center")

TSSdisttrans = calcFeatureDist(queryList, transposons)
plotFeatureDist(TSSdisttrans, featureName="Transposons", , size = 4000, nbins = 10, 
                tile = TRUE, labelOrder = "center")
#### caclulate frequency
gp = calcPartitions(query_all, partitionList, bpProportion=TRUE)
plotPartitions(gp)
gp_clean <- gp %>%
  filter(partition != "intergenic")
gp_clean$names <- c("Promoter Core", "Promoter Prox", "Exon", "Intron", "CDS",
                    "Repeat Region", "rRNA", "Centromere")
#These values were calculated based on the permutated results (obs/exp)
gp_clean$fold_change <- c(log2(942/777.68), log2(1858/168.776), log2(1130/809.80), 
                          log2(721/893.20), log2(685/1011.15), log2(23/5.23), log2(64/4.026),log2(84/8.53))
gp_clean$pvalue <- c(9.999e-05, 0.000999, 0.000999, 0.000999, 9.999e-05, 9.999e-05, 9.999e-05, 9.999e-05)

#gp_clean <- gp_clean  %>% 
 #arrange(desc(frequency)) %>%           
  #mutate(names = factor(names, levels = names)) 


bubble_plot_epimutations <- ggplot(gp_clean, aes(y = fold_change, x = names, 
                     size = frequency, color= pvalue)) + geom_point() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  labs(x = "Genomic Partition", y = "log2(Fold change)")+
  scale_size(name = "Frequency", range = c(3, 15))+
  scale_color_distiller(
    palette = "Spectral",
    direction = -1)+
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size= 12))
  
ggsave("bubble_plot.svg", plot = bubble_plot_epimutations, width = 18, height = 12, units = "cm", dpi = 300)
       
        #########################
         #raw sirnas bubble plot
       ##########################
#load ancestor raw counts
ancestor.file <- list.files(pattern = "A.*\\.filtered\\.bed$")
filtered.df <- read.table(ancestor.file, header=FALSE, sep="\t",
                          col.names=c("chrom","start","end", "count", 
                                      "score", "strand", "code"))
filtered.df$count <- as.integer(filtered.df$count)
filtered.df <- filtered.df %>% filter(count != 0)


filtered.gr <- GRanges(
  seqnames = filtered.df$chrom,
  ranges = IRanges(
    start = filtered.df$start,
    end   = filtered.df$end),
  count  = filtered.df$count)
filtered.gr <- filtered.gr[seqnames(filtered.gr) %in% c("I", "II", "III")]

partition_gr <- unlist(GRangesList(partitionList), use.names = FALSE)
partition_gr$partition <- rep(names(partitionList), lengths(partitionList))

hits <- findOverlaps(filtered.gr, partition_gr, ignore.strand = TRUE)

df <- data.frame(
  partition = partition_gr$partition[subjectHits(hits)],
  total_count = filtered.gr$count[queryHits(hits)]
)

gp_weighted <- aggregate(
  total_count ~ partition,
  data = df,
  sum)

names(gp_weighted)[2] <- "bpOverlap"   # or "frequency" if that is your preferred name
gp_weighted$frequency <- gp_weighted$bpOverlap / sum(gp_weighted$bpOverlap)

plotPartitions(gp_weighted)

#These values were calculated based on the permutated results
gp_weighted$fold_change <- c(log2(1742359/25516969), log2(38823178/20381238), log2(2439329/22460979),
                             log2(26826713/19563252),log2(36689266/42582684), log2(2363572/214416),  
                            log2(56181/125421.8),  log2(24554460/98241.76))
                        
gp_weighted$pvalue <- c(9.999e-05, 0.001, 9.999e-05, 0.1, 0.1, 0.022, 0.03, 0.000999)

gp_weighted <- gp_weighted %>%
  mutate(
    partition = recode(
      partition,
      "exon" = "Exon",
      "intron" = "Intron",
      "promoterCore" = "Prom Core",
      "promoterProx" = "Prom Prox",
      "regional_centromere" = "Centromere",
      "repeat_region" = "Rep Region"))

bubble_plot_raw <- ggplot(gp_weighted, aes(y = fold_change, x = partition, 
                            size = frequency, color= pvalue)) + geom_point() +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1)+
  labs(x = "Genomic Partition", y = "log2(Fold change)")+
  scale_size(name = "Frequency", range = c(3, 15))+
  scale_color_distiller(
    palette = "Spectral",
    direction = -1)+
  theme_minimal() +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size= 20))

ggsave("bubble_plot_raw.svg", plot = bubble_plot_raw, width = 40, height = 12, units = "cm", dpi = 300)

#control and ethanol separated
gp2 = calcPartitions(queryList, partitionList, bpProportion=TRUE)
plotPartitions(gp2)
###########################################################################################
promoterCore_data <- (partitionList[[1]])
promoterProx_data <- (partitionList[[2]])
exon_data <- (partitionList[[3]])
intron_data <- (partitionList[[4]])
CDS_data <- (partitionList[[5]])
repeat_region_data <- (partitionList[[6]])
rRNA_data <- (partitionList[[7]])
centromere_data <- (partitionList[[8]])

##############################################################
#RunregioneR
gffSource <- rtracklayer::import("/home/bioc1877-local/Documents/MA_pombe/reference/Schizosaccharomyces_pombe_all_chromosomes.gff3", format = "gff3")
transdf <- rtracklayer::import("~/Documents/MA_pombe/reference/transposons.bed")
common_seqs <- intersect(seqlevels(query), seqlevels(gffSource))
gffSource <- keepSeqlevels(gffSource, common_seqs, pruning.mode="coarse")
seqlevels(gffSource)

SPChromGR <- GRanges(
  seqnames = names(SPChromSizes),
  ranges = IRanges(start = 1L, end = as.integer(SPChromSizes))
)

SPChromGR
seqlengths(SPChromGR) <- width(SPChromGR)

pt <- permTest(A= query, ntimes = 10000, randomize.function = randomizeRegions, B=gffSource, 
               evaluate.function = numOverlaps, genome = SPChromGR, alternative = "auto")
plot(pt)

window = 10*mean(width(query))
step = mean(width(query))/2
lz <- localZScore(A=query, pt=pt, B=gffSource, window = window, step = step)
plot(lz)


#Do the permutation with the promoter core
pt_promoterCore <- permTest(A=query_all , ntimes = 10000, randomize.function = randomizeRegions, 
                     B=promoterCore_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                     alternative = "auto", verbose = TRUE)

plot(pt_promoterCore)
summary(pt_promoterCore)
mean(pt_promoterCore[[1]]$permuted)
pt_promoterCore[[1]]$observed
#Do the permutation with the promoter prox
pt_promoterProx <- permTest(A= query_all, ntimes = 1000, randomize.function = randomizeRegions, 
                            B=promoterProx_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                            alternative = "auto", verbose = TRUE)

plot(pt_promoterProx)
summary(pt_promoterProx)
mean(pt_promoterProx[[1]]$permuted)
pt_promoterProx[[1]]$observed
#Do the permutation with the exon
pt_exon <- permTest(A= query_all, ntimes = 1000, randomize.function = randomizeRegions, 
                    B= exon_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                    alternative = "auto", verbose = TRUE)

plot(pt_exon)
summary(pt_exon)
mean(pt_exon[[1]]$permuted)
pt_exon[[1]]$observed
#Do the permutation with the intron
pt_intron <- permTest(A= query_all, ntimes = 1000, randomize.function = randomizeRegions, 
                    B=intron_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                    alternative = "auto", verbose = TRUE)
plot(pt_intron)
summary(pt_intron)
mean(pt_intron[[1]]$permuted)
pt_intron[[1]]$observed
#Do the permutation with the CDS
pt_CDS <- permTest(A= query_all, ntimes = 10000, randomize.function = randomizeRegions, 
                   B=CDS_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                   alternative = "auto", verbose = TRUE)
plot(pt_CDS)
summary(pt_CDS)
mean(pt_CDS[[1]]$permuted)
pt_CDS[[1]]$observed
#Do the permutation with the repeat region
pt_repeat_region <- permTest(A= query_all, ntimes = 10000, randomize.function = randomizeRegions, 
                      B=repeat_region_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                      alternative = "auto", verbose = TRUE)
plot(pt_repeat_region)
summary(pt_repeat_region)
mean(pt_repeat_region[[1]]$permuted)
pt_repeat_region[[1]]$observed
#Do the permutation with the rRNA
pt_rRNA <- permTest(A= query_all, ntimes = 10000, randomize.function = randomizeRegions, 
                    B=rRNA_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                    alternative = "auto", verbose = TRUE)
plot(pt_rRNA)
summary(pt_rRNA)
mean(pt_rRNA[[1]]$permuted)
pt_rRNA[[1]]$observed
#Do the permutation with the centromere
pt_centromere <- permTest(A= query_all, ntimes = 10000, randomize.function = randomizeRegions, 
                      B=centromere_data, evaluate.function = numOverlaps, genome = SPChromGR, 
                      alternative = "auto", verbose = TRUE)

plot(pt_centromere)
summary(pt_centromere)
mean(pt_centromere[[1]]$permuted)
pt_centromere[[1]]$observed
####################################################################################################
# Count raw siRNA data
# Build indexed region sets


  weighted_overlap_stat <- function(A, B, ...) {
    hits <- findOverlaps(A, B, ignore.strand = TRUE)
    if (length(hits) == 0) return(0)
    sum(mcols(A)$total_count[queryHits(hits)])
  }
  
  randomize_keep_counts <- function(A, genome, ...) {
    w <- mcols(A)$total_count
    A_rand <- randomizeRegions(A, genome = genome, ...)
    mcols(A_rand)$total_count <- w
    A_rand
  }

# Run permTest for each region set
  promoterCore_pt <- permTest(A = filtered.gr, B = promoterCore_data, ntimes = 10000, 
                             randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                             genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")

  
   plot(promoterCore_pt)
  summary(promoterCore_pt)
  mean(promoterCore_pt[[1]]$permuted)
  promoterCore_pt[[1]]$observed
  
  promoterProx_pt <- permTest(A = filtered.gr, B = promoterProx_data, ntimes = 10000, 
                              randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                              genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(promoterProx_pt)
  summary(promoterProx_pt)
  mean(promoterProx_pt[[1]]$permuted)
  promoterProx_pt[[1]]$observed
  
  exon_data_pt <- permTest(A = filtered.gr, B = exon_data, ntimes = 10000, 
                              randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                              genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(exon_data_pt)
  summary(exon_data_pt)
  mean(exon_data_pt[[1]]$permuted)
  exon_data_pt[[1]]$observed
  
 intron_data_pt <- permTest(A = filtered.gr, B = intron_data, ntimes = 10000, 
                           randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                           genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(intron_data_pt)
  summary(intron_data_pt)
  mean(intron_data_pt[[1]]$permuted)
  intron_data_pt[[1]]$observed
  
  CDS_data_pt <- permTest(A = filtered.gr, B = CDS_data, ntimes = 10000, 
                             randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                             genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(CDS_data_pt)
  summary(CDS_data_pt)
  mean(CDS_data_pt[[1]]$permuted)
  CDS_data_pt[[1]]$observed
  
  
  repeat_region_pt <- permTest(A = filtered.gr, B = repeat_region_data, ntimes = 10000, 
                          randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                          genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(repeat_region_pt)
  summary(repeat_region_pt)
  mean(repeat_region_pt[[1]]$permuted)
  repeat_region_pt[[1]]$observed
  ####
  rRNA_pt <- permTest(A = filtered.gr, B = rRNA_data, ntimes = 1000, 
                               randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                               genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(rRNA_pt)
  summary(rRNA_pt)
  mean(rRNA_pt[[1]]$permuted)
  rRNA_pt[[1]]$observed
  
  centromere_pt <- permTest(A = filtered.gr, B = centromere_data, ntimes = 10000, 
                               randomize.function = randomize_keep_counts, evaluate.function = weighted_overlap_stat, 
                               genome = SPChromGR, force.parallel = FALSE, verbose = TRUE, alternative = "auto")
  
  plot(centromere_pt)
  summary(centromere_pt)
  mean(centromere_pt[[1]]$permuted)
  centromere_pt[[1]]$observed
  

  