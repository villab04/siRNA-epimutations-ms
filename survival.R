library(survival)
library(survminer)
library(fuzzyjoin)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(cowplot)
setwd("/home/bioc1877-local/Documents/MA_pombe/sRNApombe")
load("BinaryTablesAll.RData")
#Define the function
EpiLength <- function(dataIn){
  
  status <- dataIn[1]
  epimutlength <- c()
  lengthStore <- 1
  censor <- c()
  countable <- status != 0
  
  for(i in 2:length(dataIn)){
    if(!countable){
      if(dataIn[i] != 0){
        status <- dataIn[i]
        lengthStore <- 1
        countable <- TRUE
      }
    } else {
      if(dataIn[i] == status){
        lengthStore <- lengthStore + 1
      } else {
        epimutlength <- c(epimutlength, lengthStore)
        censor <- c(censor, 0)
        status <- dataIn[i]
        lengthStore <- 1}}
  }
  
  if(countable){
    epimutlength <- c(epimutlength, lengthStore)
    censor <- c(censor, 1)
  }
  
  return(rbind(epimutlength, censor))
}

head(BinaryTablesAll)

EpiMutLengths_L1<-c()
for(i in 1:nrow(BinaryTablesAll[[1]]))
  {EpiMutLengths_L1 <-cbind(EpiMutLengths_L1,EpiLength(BinaryTablesAll[[1]][i,]))}
EpiMutLengths_L1<-as.data.frame(t(EpiMutLengths_L1))

sum(EpiMutLengths_L1[["epimutlength"]] != 0)


EpiMutLengths_L2<-c()
for(i in 1:nrow(BinaryTablesAll[[2]]))
{EpiMutLengths_L2 <-cbind(EpiMutLengths_L2,EpiLength(BinaryTablesAll[[2]][i,]))}
EpiMutLengths_L2<-as.data.frame(t(EpiMutLengths_L2))

sum(EpiMutLengths_L2[["epimutlength"]] != 0)

EpiMutLengths_L6<-c()
for(i in 1:nrow(BinaryTablesAll[[3]]))
{EpiMutLengths_L6 <-cbind(EpiMutLengths_L6,EpiLength(BinaryTablesAll[[3]][i,]))}
EpiMutLengths_L6<-as.data.frame(t(EpiMutLengths_L6))

sum(EpiMutLengths_L6[["epimutlength"]] != 0)

EpiMutLengths_L7<-c()
for(i in 1:nrow(BinaryTablesAll[[4]]))
{EpiMutLengths_L7 <-cbind(EpiMutLengths_L7,EpiLength(BinaryTablesAll[[4]][i,]))}
EpiMutLengths_L7<-as.data.frame(t(EpiMutLengths_L7))

sum(EpiMutLengths_L7[["epimutlength"]] != 0)

SurvData<-list(EpiMutLengths_L1, EpiMutLengths_L2, EpiMutLengths_L6, EpiMutLengths_L7)
for(i in 1:4){colnames(SurvData[[i]])<-c("time","status")
SurvData[[i]][,2]<-SurvData[[i]][,2]*(-1)+1
SurvData[[i]][,1]<-SurvData[[i]][,1]+1
}

for(i in 1:2){SurvData[[i]]<-data.frame(SurvData[[i]],rep("Cont",length=nrow(SurvData[[i]])))}
for(i in 3:4){SurvData[[i]]<-data.frame(SurvData[[i]],rep("Ethanol",length=nrow(SurvData[[i]])))}
for(i in 1:4){colnames(SurvData[[i]])[3]<-"condition"}

TotalSurvData<-c()
for(i in 1:4){TotalSurvData<-rbind(TotalSurvData,SurvData[[i]])}

EpiMutsurvival<-survfit(Surv(time,status)~condition, data=TotalSurvData)

p <- ggsurvplot(
  EpiMutsurvival,
  conf.int = TRUE,
  risk.table = TRUE,
  pval = TRUE)

 gg <- p$plot
 
 midlife_plot <- gg + theme_minimal() +
   labs(x = "Time", y = "Survival\n Probability")+
   scale_color_manual(values=c("#298c8c","#ffbb6f"))+
   scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
   theme(legend.position = "none",
         text = element_text(size = 15))
 
midlife_table <- p$table + labs(title =NULL, y="Condition") +
   scale_y_discrete(
     labels = c("Ethanol 4%", ))+ 
    theme_gray()+
    theme(text = element_text(size = 15))
 
 p$plot <- p$plot + theme_classic()
 p$table <- p$table + theme_classic()
################################################################################
################################################################################
#combine all the list into a matrix
load("BinaryTablesAll.RData")
mat_epi <- do.call(cbind, BinaryTablesAll)
colnames(mat_epi) <- gsub("\\.windows\\.bed_norm", "", colnames(mat_epi))
head(mat_epi)
colSums(mat_epi != 0)

#remove all the rows where every column is 0
filtered_mat <- mat_epi[apply(mat_epi, 1, function(row) any(row != 0)), ]
all.epi <- do.call(rbind, strsplit(rownames(filtered_mat), ":"))
write.table(all.epi, file="all_epi.bed", sep ="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


#This is what I need to get the list of sites with epimutations
#Filter only the control lines L1 and L2
control <- filtered_mat[, grepl("L1|L2", colnames(filtered_mat))]
control <- control[apply(control, 1, function(row) any(row != 0)), ]
control_list <- do.call(rbind, strsplit(rownames(control), ":"))
write.table(control_list, file="all_control.bed", sep ="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

sum(control[["epimutlength"]] != 0)

ethanol <- filtered_mat[, grepl("L6|L7", colnames(filtered_mat))]
ethanol <- ethanol[apply(ethanol, 1, function(row) any(row != 0)), ]
ethanol_list <- do.call(rbind, strsplit(rownames(ethanol), ":"))
write.table(ethanol_list, file="all_ethanol.bed", sep ="\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Now make it a loop for the control lines
transfers <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
control <- list(c("L1T1", "L2T1"),  c("L1T2", "L2T2"),
                c("L1T3", "L2T3"),  c("L1T4", "L2T4"),
                c("L1T5", "L2T5"),  c("L1T6", "L2T6"),
                c("L1T7", "L2T7"),  c("L1T8", "L2T8"),
                c("L1T9", "L2T9"),  c("L1T10", "L2T10"))

#Create a file per transfer in the control
#the code finds genomic regions that have a nonzero value 
#in at least one of the two relevant control samples, converts their 
#row names into BED-style columns, and writes them to a separate file.
for (i in seq_along(control)) {
  cols <- intersect(control[[i]], colnames(filtered_mat))      
  if (length(cols) == 0) next                                  
  tx <- filtered_mat[, cols, drop = FALSE]
  keep <- rowSums(!is.na(tx) & tx != 0) > 0                    
  tx <- tx[keep, , drop = FALSE]
  if (nrow(tx) == 0) next
  
  out <- do.call(rbind, strsplit(rownames(tx), ":", fixed = TRUE))
  write.table(out, file = paste0("control_T", transfers[i], ".bed"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#Now a loop for the ethanol lines
ethanol <- list(c("L6T1", "L7T1"),  c("L6T2", "L7T2"),
                c("L6T3", "L7T3"),  c("L6T4", "L7T4"),
                c("L6T5", "L7T5"),  c("L6T6", "L7T6"),
                c("L6T7", "L7T7"),  c("L6T8", "L7T8"),
                c("L6T9", "L7T9"),  c("L6T10", "L7T10"))

for (i in seq_along(ethanol)) {
  cols <- intersect(ethanol[[i]], colnames(filtered_mat))      
  if (length(cols) == 0) next                                  
  txe <- filtered_mat[, cols, drop = FALSE]
  keep <- rowSums(!is.na(txe) & txe != 0) > 0                    
  txe <- txe[keep, , drop = FALSE]
  if (nrow(txe) == 0) next
  out <- do.call(rbind, strsplit(rownames(txe), ":", fixed = TRUE))
  write.table(out, file = paste0("ethanol_T", transfers[i], ".bed"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

##Take the mutations that survive for at least one transfer. 
pairwise <- list(c("L1T1", "L1T2"), c("L1T2", "L1T3"), c("L1T3", "L1T4"), c("L1T4", "L1T5"), c("L1T5", "L1T6"),
                 c("L1T6", "L1T7"), c("L1T7", "L1T8"), c("L1T8", "L1T9"), c("L1T9", "L1T10"),
                 c("L2T1", "L2T2"), c("L2T2", "L2T3"), c("L2T3", "L2T4"), c("L2T4", "L2T5"), c("L2T5", "L2T6"),
                 c("L2T6", "L2T7"), c("L2T7", "L2T8"), 
                 c("L6T1", "L6T2"), c("L6T2", "L6T3"), c("L6T3", "L6T4"), c("L6T4", "L6T5"), c("L6T5", "L6T6"),
                 c("L6T6", "L6T7"), c("L6T7", "L6T8"), c("L6T8", "L6T9"), c("L6T9", "L6T10"),
                 c("L7T1", "L7T2"), c("L7T2", "L7T3"), c("L7T3", "L7T4"), c("L7T4", "L7T5"), c("L7T5", "L7T6"),
                 c("L7T6", "L7T7"), c("L7T7", "L7T8"), c("L7T8", "L7T9"), c("L7T9", "L7T10"))

pairwise2 <- list(c("L1T1", "L1T2", "L1T3"), c("L1T2", "L1T3", "L1T4"), c("L1T3", "L1T4", "L1T5"), 
                 c("L1T4", "L1T5", "L1T6"), c("L1T5", "L1T6", "L1T7"), c("L1T6", "L1T7", "L1T8"), 
                 c("L1T7", "L1T8", "L1T9"), c("L1T8", "L1T9", "L1T10"), 
                 c("L2T1", "L2T2", "L2T3"), c("L2T2", "L2T3", "L2T4"), c("L2T3", "L2T4", "L2T5"), 
                 c("L2T4", "L2T5", "L2T6"), c("L2T5", "L2T6", "L2T7"), c("L2T6", "L2T7", "L2T8"),
                 c("L6T1", "L6T2", "L6T3"), c("L6T2", "L6T3", "L6T4"), c("L6T3", "L6T4", "L6T5"), 
                 c("L6T4", "L6T5", "L6T6"), c("L6T5", "L6T6", "L6T7"), c("L6T6", "L6T7", "L6T8"), 
                 c("L6T7", "L6T8", "L6T9"), c("L6T8", "L6T9", "L6T10"),
                 c("L7T1", "L7T2", "L7T3"), c("L7T2", "L7T3", "L7T4"), c("L7T3", "L7T4", "L7T5"),
                 c("L7T4", "L7T5", "L7T6"), c("L7T5", "L7T6", "L7T7"),
                 c("L7T6", "L7T7", "L7T8"), c("L7T7", "L7T8", "L7T9"), c("L7T8", "L7T9", "L7T10"))

##Take the mutations that survive for at least one transfer. 
for (i in seq_along(pairwise)) {
  pair <- pairwise[[i]]
  Lx <- as.data.frame(mat_epi[, pair, drop = FALSE])
  Lx <- Lx[apply(Lx, 1, function(row) any(row != 0)), ]
  Lx <- dplyr::filter(Lx, .data[[pair[1]]] == .data[[pair[2]]])
  print(paste("Comparing", pair[1], "vs", pair[2]))
  print(colSums(Lx != 0))
  if (nrow(Lx) > 0) {
    Lx_list <- do.call(rbind, strsplit(rownames(Lx), ":"))
    out_file <- paste0("stable_", pair[2], ".bed")
    write.table(Lx_list, file = out_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

#Take the transitions at least one transfer
for (i in seq_along(pairwise)) {
  pair <- pairwise[[i]]
  Lx <- as.data.frame(mat_epi[, pair, drop = FALSE])
  Lx <- dplyr::filter(Lx, .data[[pair[1]]] != .data[[pair[2]]])
  print(paste("Comparing", pair[1], "vs", pair[2]))
  if (nrow(Lx) > 0) {
    Lx_list <- do.call(rbind, strsplit(rownames(Lx), ":"))
    out_file <- paste0("transition_", pair[2], ".bed")
    write.table(Lx_list, file = out_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

#Take the transitions at that last at least one generation
#I dont understand what I intended with this, but I think is wrong
for (i in seq_along(pairwise2)) {
  pair <- pairwise2[[i]]
  Lx <- as.data.frame(mat_epi[, pair, drop = FALSE])
  Lx <- dplyr::filter(Lx, .data[[pair[1]]] != .data[[pair[2]]] & .data[[pair[3]]])
  print(paste("Comparing", pair[1], "vs", pair[2], "and", pair[3]))
  if (nrow(Lx) > 0) {
    Lx_list <- do.call(rbind, strsplit(rownames(Lx), ":"))
    out_file <- paste0("stabletransition_", pair[2], ".bed")
    write.table(Lx_list, file = out_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}


#Check for gains compared with previous transfer
for (i in seq_along(pairwise)) {
  pair <- pairwise[[i]]
  Lx <- as.data.frame(mat_epi[, pair, drop = FALSE])
  Lx <- dplyr::filter(Lx, .data[[pair[1]]] == 0 | .data[[pair[1]]] == -1)
  Lx <- dplyr::filter(Lx, .data[[pair[2]]] == 1)
  print(paste("Comparing", pair[1], "vs", pair[2]))
  print(colSums(Lx == 1))
  if (nrow(Lx) > 0) {
    Lx_list <- do.call(rbind, strsplit(rownames(Lx), ":"))
    out_file <- paste0("gains_", pair[2], ".bed")
    write.table(Lx_list, file = out_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

sum(mat_epi[, "L1T1"] == 1)
sum(mat_epi[, "L2T1"] == 1)
sum(mat_epi[, "L6T1"] == 1)
sum(mat_epi[, "L7T1"] == 1)
#Check for losses compared with previous transfer
for (i in seq_along(pairwise)) {
  pair <- pairwise[[i]]
  Lx <- as.data.frame(mat_epi[, pair, drop = FALSE])
  Lx <- dplyr::filter(Lx, .data[[pair[1]]] %in% c(0,1) & .data[[pair[2]]] == -1 |
                          .data[[pair[1]]] == 1 & .data[[pair[2]]] == 0)
  print(paste("Comparing", pair[1], "vs", pair[2]))
  print(dim(Lx))
  if (nrow(Lx) > 0) {
    Lx_list <- do.call(rbind, strsplit(rownames(Lx), ":"))
    out_file <- paste0("loss_", pair[2], ".bed")
    write.table(Lx_list, file = out_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

sum(mat_epi[, "L1T1"] == -1)
sum(mat_epi[, "L2T1"] == -1)
sum(mat_epi[, "L6T1"] == -1)
sum(mat_epi[, "L7T1"] == -1)

#create a database with the count of epimutations per 10kb windows
df_number <- as.data.frame(rowSums(mat_epi != 0))
colnames(df_number) <- c("count")
df_complete <- as.data.frame(do.call(rbind, strsplit(rownames(df_number), ":")))
colnames(df_complete) <- c("chrom", "start", "end")
df_complete$count <- df_number$count
df_complete <- df_complete %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                                  "mating_type_region", "mitochondrial"))
df_complete$start <- as.numeric(df_complete$start)
df_complete$end <- as.numeric(df_complete$end)

norm_counts <- read.csv("normalized_counts_windows.csv")
sum_counts <- norm_counts %>% rowwise() %>%  mutate(s = sum(c_across(4:53)))
sum_counts <- sum_counts %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                                  "mating_type_region", "mitochondrial"))
#create counts from Ancestor
anc_counts <- norm_counts %>% dplyr::select(c(chrom, start, end, count.ANC.windows.bed))

sum_count <- data.frame(
  chrom = sum_counts$chrom,
  start = as.numeric(sum_counts$start),
  end   = sum_counts$end,
  sum  = as.numeric(sum_counts$s)
)
#############################################################################
site <- as.data.frame(L1_order['I:1314000:1315000',])
colnames(site) <- c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
rownames(site) <- "value"
site_t <- as.data.frame(t(site))
site_t$transfer <- transfer

site_t$transfer <- factor(site_t$transfer, levels = paste0("T", 1:10))

methods_plot2 <- ggplot(site_t, aes(x = transfer, y = value)) +
  geom_point() + geom_line(aes(group = 1)) +
  labs(x="Transfer")+ theme_minimal() + 
  theme(text = element_text(size= 25),
        axis.title.y = element_blank())

ggsave(filename = "methods_plot2.svg", plot = methods_plot2,
       width = 15, height = 7, units = "cm")

##############################################################################
#Create databases of the count of epigenetic changes of control and ethanol separated
control1 <- mat_epi[, grepl("L1|L2", colnames(mat_epi))]
control_number <- as.data.frame(rowSums(control1 != 0))
colnames(control_number) <- c("count")
control_complete <- as.data.frame(do.call(rbind, strsplit(rownames(control_number), ":")))
colnames(control_complete) <- c("chrom", "start", "end")
control_complete$count <- control_number$count
control_complete <- control_complete %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                                    "mating_type_region", "mitochondrial"))
control_complete$start <- as.integer(control_complete$start)
control_complete$end <- as.integer(control_complete$end)

ethanol1 <- mat_epi[, grepl("L5|L6", colnames(mat_epi))]
ethanol_number <- as.data.frame(rowSums(ethanol1 != 0))
colnames(ethanol_number) <- c("count")
ethanol_complete <- as.data.frame(do.call(rbind, strsplit(rownames(ethanol_number), ":")))
colnames(ethanol_complete) <- c("chrom", "start", "end")
ethanol_complete$count <- ethanol_number$count
ethanol_complete <- ethanol_complete %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                                              "mating_type_region", "mitochondrial"))
ethanol_complete$start <- as.integer(ethanol_complete$start)
ethanol_complete$end <- as.integer(ethanol_complete$end)

#Know the amount of epigenetic changes fall inside heterochromatin regions
heterochromatin <- read.table("~/Documents/MA_pombe/reference/centromeresandtelomeres.bed")
colnames(heterochromatin) <- c("chrom", "start", "end")
        
gr_windowscontrol <- GRanges(
  seqnames = control_complete$chrom,
  ranges   = IRanges(start = control_complete$start,
                     end   = control_complete$end),
  count = control_complete$count)

gr_windowsethanol <- GRanges(
  seqnames = ethanol_complete$chrom,
  ranges   = IRanges(start = ethanol_complete$start,
                     end   = ethanol_complete$end),
  count = ethanol_complete$count)

gr_hetero  <- rtracklayer::import("~/Documents/MA_pombe/reference/centromeresandtelomeres.bed")

    #count the number windows of overlap with heterochromatin regions and the number of non overlap
hits_control <- findOverlaps(gr_windowscontrol, gr_hetero, ignore.strand = TRUE)
overlap_windows_control <- gr_windowscontrol[unique(queryHits(hits_control))]
sum(mcols(overlap_windows_control)$count != 0)
nonoverlap_windows_control <- gr_windowscontrol[ !overlapsAny(gr_windowscontrol, gr_hetero) ]
sum(mcols(nonoverlap_windows_control)$count != 0)
sum(mcols(gr_windowscontrol)$count != 0)

hits_ethanol <- findOverlaps(gr_windowsethanol, gr_hetero, ignore.strand = TRUE)
overlap_windows_ethanol <- gr_windowsethanol[unique(queryHits(hits_ethanol))]
sum(mcols(overlap_windows_ethanol)$count != 0)
nonoverlap_windows_ethanol <- gr_windowsethanol[ !overlapsAny(gr_windowsethanol, gr_hetero) ]
sum(mcols(nonoverlap_windows_ethanol)$count != 0)
sum(mcols(gr_windowsethanol)$count != 0)


# extract windows that overlap at least one heterochromatin region
windows_overlapping <- gr_windows[unique(queryHits(hits))]

# write back to bed
rtracklayer::export(windows_overlapping, "windows_overlapping_heterochromatin.bed")

density_graph <- ggplot() +
  geom_density(data = ethanol_number, aes(x = log10(count)), color = FALSE, fill = "#298c8c", alpha = 1) +
  geom_density(data = control_number, aes(x = log10(count)), color = FALSE, fill = "#ffbb6f", alpha = 0.8)+
  labs(x = "log10 Epimutation count", y = "Density") +
  theme_minimal()+
  theme(text = element_text(size = 15))

ggsave(filename = "density_plot.svg", plot = density_graph, units = "cm", width = 10, height = 15)

###########################################################################
centromere <- read.table("~/Documents/MA_pombe/reference/centromere.bed")
colnames(centromere) <- c("chrom", "start", "end")


transposons <- read.table("~/Documents/MA_pombe/reference/transposons.bed")
transposons <- as.data.frame(transposons)
colnames(transposons) <- c("chrom", "start", "end")

df_chrom <- data.frame(chrom, start, end)

siRNA_cons <- read.table("~/Documents/MA_pombe/sRNApombe/constant_sRNA_deseq_results.bed")
colnames(siRNA_cons) <- c("chrom", "start", "end", "value")
rownames(siRNA_cons) <- NULL
siRNA_cons <- siRNA_cons %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                                    "mating_type_region", "mitochondrial"))

chipseq <- read.table("~/Documents/MA_pombe/reference/peaks_peaks.narrowPeak")
chipseq <- chipseq %>% dplyr::select(V1, V2, V3)
colnames(chipseq) <- c("chrom", "start", "end")
chipseq <- chipseq %>% filter(!chrom %in% c("chr_II_telomeric_gap", 
                                               "mating_type_region", "mitochondrial"))

load("epimutation_number_results.RData")

wide_df <- result_df %>%
  dplyr::select(window, line, number_of_epimutation) %>%
  pivot_wider(
    names_from = line,
    values_from = number_of_epimutation,
    values_fill = 0)

wide_df$total <- rowSums(wide_df[, c("L1", "L2", "L6", "L7")],na.rm = TRUE) 
wide_df <- tidyr::separate(wide_df, col = window,into = c("chrom", "start", "end"), sep = ":")
counts.df <- merge(wide_df, anc_counts, by = c("chrom", "start", "end"), all.y = TRUE)
counts.df <- counts.df %>% mutate(total = replace_na(total, 0))

  anI   <- anc_counts   %>% filter(chrom == "I")
  dcI   <- counts.df %>% filter(chrom == "I")
  siI   <- siRNA_cons   %>% filter(chrom == "I")
  chipI <- chipseq      %>% filter(chrom == "I")
  cenI  <- centromere   %>% filter(chrom == "I")
  teI   <- transposons  %>% filter(chrom == "I")
  
  anII   <- anc_counts   %>% filter(chrom == "II")
  dcII   <- counts.df %>% filter(chrom == "II")
  siII   <- siRNA_cons   %>% filter(chrom == "II")
  chipII <- chipseq      %>% filter(chrom == "II")
  cenII  <- centromere   %>% filter(chrom == "II")
  teII   <- transposons  %>% filter(chrom == "II")
  
  anIII   <- anc_counts   %>% filter(chrom == "III")
  dcIII   <- counts.df %>% filter(chrom == "III")
  siIII   <- siRNA_cons   %>% filter(chrom == "III")
  chipIII <- chipseq      %>% filter(chrom == "III")
  cenIII  <- centromere   %>% filter(chrom == "III")
  teIII   <- transposons  %>% filter(chrom == "III")
  

  p1 <- ggplot(anI) +
    geom_line(aes(x = start, y = log10(count.ANC.windows.bed+1)), color = "#999999")+
    theme_minimal()+ labs(y = "log10\n(count+1)")+
    theme(
      axis.title.x = element_blank(),  
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 9)
    )
  p2 <- ggplot(dcI, aes(x = as.numeric(start), y = log10(total+1))) +
    geom_line(color = "tomato3", linewidth = 0.5) +
    theme_minimal()+ labs(y = "log10\n(count+1)") +
    theme(
      axis.title.x = element_blank(),  
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 9))
  
   p3 <- ggplot()+ geom_rect(data = cenI,
              aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
              fill = "#61D04F",  color = "#61D04F", linewidth = 0.3) +
    geom_rect(data = teI,
              aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
              fill = "blue",  color = "blue", linewidth = 0.3) +
    geom_rect(data = siI,
              aes(xmin = start, xmax = end, ymin = 0, ymax = 0.3),
              fill = "black",  color = "black", linewidth = 0.4) +
    geom_rect(data = chipI,
              aes(xmin = start, xmax = end, ymin = 0.35, ymax = 0.65),
              fill = "orange", color = "orange", linewidth = 0.3) + theme_minimal() +
     theme(
       axis.title.y  = element_blank(),  
       axis.text.y = element_blank(),
       axis.ticks.y = element_blank() 
     )
  
 chrom1 <-  p1 / p2 / p3

   pII <- ggplot(anII) +
     geom_line(aes(x = start, y = log10(count.ANC.windows.bed)), color = "#999999")+
     theme_minimal()+ labs(y = "log10\n(count+1)") +
     theme(
       axis.title.x = element_blank(),  
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y = element_text(size = 9))
   
   p2II <- ggplot(dcII, aes(x = as.numeric(start), y = log10(total+1))) +
     geom_line(color = "tomato3", linewidth = 0.5) +
     theme_minimal()+ labs(y = "log10\n(count+1)") +
     theme(
       axis.title.x = element_blank(),  
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y = element_text(size = 9))
   
   p3II <- ggplot()+ geom_rect(data = cenII,
                             aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
                             fill = "#61D04F", color = "#61D04F", linewidth = 0.3 ) +
     geom_rect(data = teII,
               aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
               fill = "blue", color = "blue", linewidth = 0.3) +
     geom_rect(data = siII,
               aes(xmin = start, xmax = end, ymin = 0, ymax = 0.3),
               fill = "black", color = "black", linewidth = 0.4) +
     geom_rect(data = chipII,
               aes(xmin = start, xmax = end, ymin = 0.35, ymax = 0.65),
               fill = "orange", color = "orange", linewidth = 0.3) + theme_minimal() +
     theme(
       axis.title.y  = element_blank(),  
       axis.text.y = element_blank(),
       axis.ticks.y = element_blank() 
     )

   chrom2 <-  pII / p2II / p3II

   pIII <- ggplot(anIII) +
     geom_line(aes(x = start, y = log10(count.ANC.windows.bed)), color = "#999999")+
     theme_minimal()+ labs(y = "log10\n(count+1)") +
     theme(
       axis.title.x = element_blank(),  
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y = element_text(size = 9),
     )
   p2III <- ggplot(dcIII, aes(x = as.numeric(start), y = log10(total+1))) +
     geom_line(color = "tomato3", linewidth = 0.5) +
     theme_minimal()+ labs(y = "log10\n(count+1)") +
     theme(
       axis.title.x = element_blank(),  
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.y = element_text(size = 9))
   
   p3III <- ggplot()+ geom_rect(data = cenIII,
                               aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
                               fill = "#61D04F", color = "#61D04F", linewidth = 0.3) +
     geom_rect(data = teIII,
               aes(xmin = start, xmax = end, ymin = 0.7, ymax = 1),
               fill = "blue", color = "blue", linewidth = 0.3) +
     geom_rect(data = siIII,
               aes(xmin = start, xmax = end, ymin = 0, ymax = 0.3),
               fill = "black", color = "black", linewidth = 0.4) +
     geom_rect(data = chipIII,
               aes(xmin = start, xmax = end, ymin = 0.35, ymax = 0.65),
               fill = "orange", color = "orange", linewidth = 0.3) + theme_minimal() + labs(x = "Genomic position")+
     theme(
       axis.title.y  = element_blank(),  
       axis.text.y = element_blank(),
       axis.ticks.y = element_blank()
     )
   
chrom3 <-  pIII / p2III / p3III

legend_df <- data.frame(
  label = c("Ancestor", "Epi changes",
            "Transposons", "Centromeres", "H3K9me3", "Consistent siRNA differences"),
  x = seq(0.02, 0.50, length.out = 6)  # evenly spaced, closer together
)
legend_plot <- ggplot(legend_df) +
  geom_rect(aes(xmin = x, xmax = x + 0.01, ymin = 0.017, ymax = 0.05),
            fill = c("#999999","tomato3","blue","#61D04F","orange","black")) +
  geom_text(aes(x = x + 0.012, y = 0.035, label = label),
            hjust = 0, size = 3.5) +   # <-- smaller text
  coord_cartesian(xlim = c(0, 0.72), expand = FALSE) +  # <-- extra room on right
  theme_void()

my_plot <- ggdraw() +
  draw_plot(chrom3, x = 0, y = 0.07, width = 1, height = 0.31) +
  draw_plot(chrom2, x = 0, y = 0.37, width = 1, height = 0.31) +
  draw_plot(chrom1, x = 0, y = 0.67, width = 1, height = 0.31)+
  draw_plot_label(
  label = c("Chromosome III", "Chromosome II", "Chromosome I"),
    x = c(0, 0, 0), y = c(0.39, 0.69, 0.98), size = 10)+
  draw_plot(legend_plot, x = 0.1, y = 0, width = 0.9, height = 0.04)
  
ggsave(filename = "example-plot.svg", plot = my_plot,
       width = 30, height = 25, units = "cm")

Figure2 <- ggdraw() +
  draw_plot(my_plot, x=0, y=0.04, width = 1, height = 0.70)+
  draw_plot_label(
    label = c("A", "B"),
    x= c(0,0), y=c(0.97, 0.77)
  )

ggsave(filename = "Figure2.svg", plot = my_plot,
       width = 25, height = 30, units = "cm")

#transposon extraction
trans_counts <- fuzzy_inner_join(
  norm_counts, transposons,
  by = c("chrom" = "chrom",
         "start" = "start",
         "end" = "end"),
  match_fun = list(`==`, `>=`, `<=`)
)

trans_epicounts <- fuzzy_inner_join(
  df_complete, transposons,
  by = c("chrom" = "chrom",
         "start" = "start",
         "end" = "end"),
  match_fun = list(`==`, `>=`, `<=`)
)

#"Add the transposon name"
#This is also for trans_counts!!!
trans_epicounts <- trans_epicounts %>%
  mutate(region_label = case_when(
    start.x >= 5191060 & end.x <= 5195145 & chrom.x == "I" ~ "Tf2-7",
    start.x >= 1564101 & end.x <= 1568217 & chrom.x == "I" ~ "Tf2-2",
    start.x >= 5195711 & end.x <= 5199712 & chrom.x == "I" ~ "Tf2-8",
    start.x >= 3361547 & end.x <= 3365548 & chrom.x == "I" ~ "Tf2-4",
    start.x >= 4022492 & end.x <= 4026493 & chrom.x == "I" ~ "Tf2-6",
    start.x >= 2927204 & end.x <= 2931205 & chrom.x == "I" ~ "Tf2-3",
    start.x >= 1465847  & end.x <= 1469848 & chrom.x == "I" ~ "Tf2-1",
    start.x >= 3996423 & end.x <= 4000424 & chrom.x == "I"  ~ "Tf2-5",
    start.x >= 4414657 & end.x <= 4418658 & chrom.x == "II" ~ "Tf2-11",
    start.x >= 1965390  & end.x <= 1969390 & chrom.x == "II" ~ "Tf2-10",
    start.x >= 1812748 & end.x <= 1816749 & chrom.x == "II" ~ "Tf2-9",
    start.x >= 778133 & end.x <= 782134 & chrom.x == "III" ~ "Tf2-12",
    start.x >= 2320320 & end.x <= 2324320 & chrom.x == "III" ~ "Tf2-13",
    start.x >= 2926800 & end.x <= 2931725 & chrom.x == "I" ~ "SPTF.3",
    start.x >= 777729  & end.x <= 782654 & chrom.x == "III" ~ "CU329672",
    start.x >= 1812344 & end.x <= 1817267 & chrom.x == "II" ~ "SPTF.6",
    start.x >= 1964875 & end.x <= 1969789 & chrom.x == "II" ~ "CU329671",
    start.x >= 3996019 & end.x <= 4000691 & chrom.x == "I" ~ "CU329670",
    start.x >= 3361143 & end.x <= 3366062 & chrom.x == "I" ~ "SPTF.5",
    start.x >= 1964875 & end.x <= 1969789 & chrom.x == "II" ~ "CU329671",
    start.x >= 3996019 & end.x <= 4000691 & chrom.x == "I" ~ "CU329670",
    start.x >= 3361143 & end.x <= 3366062 & chrom.x == "I" ~ "SPTF.5",
    start.x >= 2939937 & end.x <= 2942286 & chrom.x == "I" ~ "SPTF.4",
    start.x >= 1465327 & end.x <= 1470252 & chrom.x == "I" ~ "SPTF.1",
    start.x >= 1563817 & end.x <= 1567862 & chrom.x == "I" ~ "SPTF.2",
    start.x >= 254621 & end.x <= 256351 & chrom.x == "III" ~ "SPTF.8",
    start.x >= 2319916 & end.x <= 2324840 & chrom.x == "III" ~ "SPTF.9",
    start.x >= 4414197 & end.x <= 4419057 & chrom.x == "II" ~ "SPTF.7",
    start.x >= 4021972 & end.x <= 4026897 & chrom.x == "I" ~ "CU329670",
    start.x >= 5194995 & end.x <= 5200227 & chrom.x == "I" ~ "CU329670",
    TRUE ~ "Other"
  ))
#change the name of the columns to remove counts windows ed
trans_counts <- trans_counts %>%
  rename_with(~ str_remove_all(., "count\\.|\\.windows\\.bed"), starts_with("count."))


trans_counts_long <- trans_counts %>%
  # Remove unwanted columns
  dplyr::select(
    -starts_with("chrom.y"),
    -starts_with("start.y"),
    -starts_with("end.y"),
    -matches("*_R2$")
  ) %>%
  # Pivot to long format
  pivot_longer(
    cols = -c(chrom.x, start.x, end.x, region_label),
    names_to = "sample",
    values_to = "count"
  )

#create the line column
trans_counts_long <- trans_counts_long %>%
  mutate(line = str_extract(sample, "L\\d+"))
#create the transfer column
trans_counts_long <- trans_counts_long %>%
  mutate(transfer = str_extract(sample, "(?<=T)\\d+"))

region_sample_summary <- trans_counts_long %>%
  group_by(region_label, sample, line, transfer) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))
region_sample_summary$transfer <- as.numeric(region_sample_summary$transfer)

region_sample_summary <- na.omit(region_sample_summary)

 ggplot(region_sample_summary, aes(x = transfer, y = log10(mean_value), color= region_label)) +
  geom_point() +
  facet_wrap(~ line ) +
  labs(x = "transfer", y = "log sum count") +
  theme_bw()
 
 df1 <- trans_counts_long %>% select(chrom.x, start.x, end.x)
 distinct(df1)
 df1 <- df1 %>%
   unite("genomic_position", chrom.x, start.x, end.x, sep = ":")
 
 df <- df %>%
   mutate(highlight = df$region %in% df1$genomic_position)
###################################################################
#MAKE THE SINGLE TRANSPOSONS GRAPH WITH THE COUNT OF EPIMUTATIONS

unique(trans_epicounts$region_label)
count_Tf21 <- trans_epicounts %>% filter(region_label == "Tf2-1") #zero epimutations  
count_Tf22 <- trans_epicounts %>% filter(region_label == "Tf2-2") # zero epimutations
count_Tf23 <- trans_epicounts %>% filter(region_label == "Tf2-3") # two(two windows)
count_Tf24 <- trans_epicounts %>% filter(region_label == "Tf2-4") # two (two windows)
count_Tf25 <- trans_epicounts %>% filter(region_label == "Tf2-5") # four (two windows)
count_Tf26 <- trans_epicounts %>% filter(region_label == "Tf2-6") # two (two windows)
count_Tf27 <- trans_epicounts %>% filter(region_label == "Tf2-7") # zero
count_Tf28 <- trans_epicounts %>% filter(region_label == "Tf2-8") # two (two windows)
count_Tf29 <- trans_epicounts %>% filter(region_label == "Tf2-9") # two(two windows)
count_Tf210 <- trans_epicounts %>% filter(region_label == "Tf2-10") #zero
count_Tf211 <- trans_epicounts %>% filter(region_label == "Tf2-11")
count_Tf212 <- trans_epicounts %>% filter(region_label == "Tf2-12") #13
count_Tf213 <- trans_epicounts %>% filter(region_label == "Tf2-13") #3

count_SPTF.1 <- trans_epicounts %>% filter(region_label == "SPTF.1") 
count_SPTF.2 <- trans_epicounts %>% filter(region_label == "SPTF.2") 
count_SPTF.3 <- trans_epicounts %>% filter(region_label == "SPTF.3") # 25 
count_SPTF.4 <- trans_epicounts %>% filter(region_label == "SPTF.4")
count_SPTF.5 <- trans_epicounts %>% filter(region_label == "SPTF.5") # 22
count_SPTF.6 <- trans_epicounts %>% filter(region_label == "SPTF.6") # 2
count_SPTF.7 <- trans_epicounts %>% filter(region_label == "SPTF.7") 
count_SPTF.8 <- trans_epicounts %>% filter(region_label == "SPTF.8") 
count_SPTF.9 <- trans_epicounts %>% filter(region_label == "SPTF.9") 

count_CU329670 <- trans_epicounts %>% filter(region_label == "CU329670") #1 epimutations
count_CU329671 <- trans_epicounts %>% filter(region_label == "CU329671") #zero epimutations
count_CU329672 <- trans_epicounts %>% filter(region_label == "CU329672") #22 epimutations

transposons.names <- c("Tf2-1", "Tf2-2", "Tf2-3", "Tf2-4", "Tf2-5", "Tf2-6", 
                       "Tf2-7", "Tf2-8", "Tf2-9", "Tf2-10", "Tf2-11", "Tf2-12", 
                       "Tf2-13", "SPTF.1", "SPTF.2", "SPTF.3", "SPTF.4", "SPTF.5", 
                       "SPTF.6", "SPTF.7", "SPTF.8", "SPTF.9", "CU329670", "CU329671", 
                       "CU329672")  
epimut <- c(0, 0, 1, 1, 2, 1, 0, 1, 1, 0, 0, 13, 3, 0, 0, 25, 0, 22, 2, 0, 0, 0, 1, 0, 22)

trans_count <- data.frame(
  transposons = transposons.names,
  epimut = epimut
)

trans_count <- trans_count %>%
  mutate(
    family = case_when(
      grepl("^Tf2", transposons) ~ "Tf2",
      grepl("^SPTF", transposons) ~ "SPTF",
      TRUE ~ "Other"
    )
  )

ggplot(trans_count, aes(x = reorder(transposons, -epimut), 
                        y = epimut, fill = family)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("Epimut count") +
  xlab("Transposon") +
  scale_fill_manual(values = c("Tf2" = "#E41A1C", 
                               "SPTF" = "#377EB8", 
                               "Other" = "#4DAF4A"))


coords <- list(
  SPTF.3    = "I:2927000:2928000",
  CU329672  = "III:778000:779000",
  SPTF.5    = "I:3365000:3366000",
  TF12      = "III:9000:10000"
)

extract_long <- function(vec, transposon_name) {
  
  df <- data.frame(
    name = names(vec),
    value = as.numeric(vec)
  ) %>%
    separate(name, into = c("L", "T"), sep = "T") %>%
    mutate(
      Line     = as.numeric(gsub("L", "", L)),
      Transfer = as.numeric(T),
      transposon = transposon_name
    ) %>%
    select(transposon, Line, Transfer, Value = value)
  
  return(df)
}

final_db <- do.call(
  rbind,
  lapply(names(coords), function(tp) {
    vec <- mat_epi[ coords[[tp]], ]
    extract_long(vec, tp)
  })
)
 
final_db$Transfer <- as.factor(final_db$Transfer)
final_db$Line <- as.factor(final_db$Line)
final_db$Value <- as.factor(final_db$Value)

final_SPTF.3 <- final_db %>% filter(final_db$transposon == "SPTF.3")
final_CU329672 <- final_db %>% filter(final_db$transposon == "CU329672")
final_SPTF.5 <- final_db %>% filter(final_db$transposon == "SPTF.5")
final_TF.12 <- final_db %>% filter(final_db$transposon == "TF12")

ggplot(final_SPTF.3, aes(x = Transfer, y = Value)) +
  geom_point() +            
  geom_line(aes(group = Line), linewidth = 0.7) + 
  facet_wrap(~ Line, ncol = 1) +                  
  theme_bw() +
  labs(
    title = "SPTF.3",
    x = "Transfer",
    y = "Value"
  )

ggplot(final_CU329672, aes(x = Transfer, y = Value)) +
  geom_point() +            
  geom_line(aes(group = Line), linewidth = 0.7) + 
  facet_wrap(~ Line, ncol = 1) +                  
  theme_bw() +
  labs(
    title = "CU329672",
    x = "Transfer",
    y = "Value"
  )

ggplot(final_SPTF.5, aes(x = Transfer, y = Value)) +
  geom_point() +            
  geom_line(aes(group = Line), linewidth = 0.7) + 
  facet_wrap(~ Line, ncol = 1) +                  
  theme_bw() +
  labs(
    title = "SPTF.5",
    x = "Transfer",
    y = "Value"
  )

ggplot(final_TF.12, aes(x = Transfer, y = Value)) +
  geom_point() +            
  geom_line(aes(group = Line), linewidth = 0.7) + 
  facet_wrap(~ Line, ncol = 1) +                  
  theme_bw() +
  labs(
    title = "TF.12",
    x = "Transfer",
    y = "Value"
  )
####################################################################
#look into the number of epimutations
data.number <- read.csv("number_epimutations.csv")
shared.number <- read.csv("number_sharedepimutations.csv")
transition_data <- read.csv("transitions_data.csv")

ggplot(transition_data, aes(condition, Transitions))+
  geom_boxplot()+
  theme_bw()

ggplot(transition_data, aes(condition, stable_transitions))+
  geom_boxplot()+
  theme_bw()

cent_transitions.plot <- ggplot(transition_data, aes(condition, centromere, fill = condition))+
  geom_boxplot()+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  labs(y = "Number of transitions\nin the centromere")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")

model1 <- lm(data = transition_data, centromere ~ condition)
summary(model1)

ggplot(transition_data, aes(condition, transposons))+
  geom_boxplot()+
  theme_bw()

model2 <- lm(data = transition_data, transposons ~ condition)
summary(model2)

epichanges.plot <- ggplot(data.number, aes(condition, total, fill = condition))+
  geom_boxplot()+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  labs(y = "Number of\nepigentic changes")+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")

epimutations.graph <- ggplot(data.number, aes(condition, stable, fill = condition))+
  geom_boxplot()+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  labs(y = "Number of\n epimutations")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 15),
        legend.position = "none",
        )

ggplot(data.number, aes(condition, proportion))+
  geom_boxplot()+
  theme_bw()

transition.graph <- ggplot(data.number, aes(condition, transition, fill = condition))+
  geom_boxplot()+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  labs(y = "Number of\ntransitions")+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")
  

model.transition <- lm(transition~ condition, data.number) 
summary(model.total)

ggplot(data.number, aes(condition, gains))+
  geom_boxplot()+
  theme_bw()

model.gains <- lm(gains ~ condition, data.number)
summary(model.gains)

ggplot(data.number, aes(condition, lost))+
  geom_boxplot()+
  theme_bw()

model.lost <- lm(lost ~ condition, data.number)
summary(model.lost)

ggplot(shared.number, aes(condition, shared_epimutation))+
  geom_boxplot()+
  theme_bw()

model.total <- lm(proportion ~ condition, data.number)

#Put all the graphs together
Figure2 <- ggdraw() +
  draw_plot(epichanges.plot, x = 0, y = 0.65, width = 0.45, height = 0.35) +
  draw_plot(epimutations.graph, x = 0.5, y = 0.65, width = 0.45, height = 0.35) +
  draw_plot(transition.graph, x = 0, y = 0.31, width = 0.45, height = 0.35)+
  draw_plot(cent_transitions.plot, x = 0.5, y = 0.31, width = 0.45, height = 0.35)+
  draw_plot(midlife_plot, x = 0.01, y = 0, width = 0.93, height = 0.3)+
  #draw_plot(midlife_table, x = 0.5, y = 0, width = 0.45, height = 0.3)+
  draw_plot_label(
    label = c("A", "B"),
    x = c(0, 0), y = c(0.98, 0.35), size = 20)

ggsave(filename = "Figure2.svg", plot = Figure2,
       width = 25, height = 25, units = "cm")
