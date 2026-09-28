library(tidyverse)
library(ggpubr)
setwd("/home/bioc1877-local/Documents/MA_pombe/sRNApombe")

# Function to calculate epimutations based on residuals from a LOESS fit
#Create the function GET EPIMUTATIONS
################################################################################
GetEpimutations <- function(fileIn, control) {
  # Read test (descendant) and control (ancestor) count files
  fileTest <- get(fileIn)
  fileCont <- get(control)
  # Create IDs for merging
  fileTest$ID <- paste(fileTest$chrom, fileTest$start, fileTest$end, sep=":")
  fileCont$ID <- paste(fileCont$chrom, fileCont$start, fileCont$end, sep=":")
  # Merge by ID (only shared windows kept)
  merged <- merge(fileTest[,c("ID","normalized_count")], 
                  fileCont[,c("ID","normalized_count")], 
                  by="ID", suffixes=c(".test",".cont"))
  # Fit a LOESS model: log2(test counts) ~ log2(control counts)
  Mod <- loess(log2(merged$normalized_count.test + 1) ~ log2(merged$normalized_count.cont + 1))$residuals
  # Convert residuals to z-scores
  #ZMod <- (Mod - mean(Mod)) / sd(Mod)
  ZMod <- (Mod - median(Mod)) / mad(Mod)
  # Add results back
  merged$Residual <- Mod
  merged$Zscore <- ZMod
  return(merged)
}
################################################################################

#Group files into lineage sets
# Get list of files for each lineage
#.* → any characters, any number of times
#$end of the file
L1 <- list.files(pattern="L1.*\\.windows.bed$")
L2 <- list.files(pattern="L2.*\\.windows.bed$")
L6 <- list.files(pattern="L6.*\\.windows.bed$")
L7  <- list.files(pattern="L7.*\\.windows.bed$")

#Load the csv of the normalized counts
norm_counts <- read.csv("normalized_counts_windows.csv")
norm_counts <- norm_counts %>% filter(count.ANC.windows.bed < 5791)
# Identify all sample columns (those starting with "count.")
sample_cols <- grep("^count\\.", names(norm_counts), value = TRUE)
# Create a list of sample-specific data.frames
sample_dfs <- lapply(sample_cols, function(col) {
  data.frame(
    chrom = norm_counts$chrom,
    start = norm_counts$start,
    end   = norm_counts$end,
    normalized_count = norm_counts[[col]]  # consistent column name
  )
})

# Name each element in the list according to the original sample column
names(sample_dfs) <- sample_cols
for(nm in names(sample_dfs)) {
  # Create a clean variable name by removing "count." and adding "_norm"
  df_name <- paste0(sub("^count\\.", "", nm), "_norm")
  # Assign the data.frame to the environment
  assign(df_name, sample_dfs[[nm]])
}

######### extract the replicates to find the threshold value
# All *_norm objects in your environment
all_norms <- ls(pattern = "_norm$")
# Keep only their names (as character vectors)
L1T2 <- grep("^L1T2.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L1T7 <- grep("^L1T7.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L1T8 <- grep("^L1T8.*\\.windows\\.bed_norm$", all_norms, value = TRUE)

L2T3 <- grep("^L2T3.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L2T7 <- grep("^L2T7.*\\.windows\\.bed_norm$", all_norms, value = TRUE)

L6T4 <- grep("^L6T4.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L6T6 <- grep("^L6T6.*\\.windows\\.bed_norm$", all_norms, value = TRUE)

L7T2 <- grep("^L7T2.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L7T4 <- grep("^L7T4.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L7T6 <- grep("^L7T6.*\\.windows\\.bed_norm$", all_norms, value = TRUE)

# Combine into list of all lineages
LineagesAll <- list(L1, L2, L6, L7) 


# Path to the control (ancestor) count file
Control <- grep("^ANC.*\\.windows\\.bed_norm$", all_norms, value = TRUE)

# Empty list to store pairwise results
PairsAllOut <- list()


BinaryTables <- list()   # initialize once

cutoffs <- c(1, 2, 3, 4, 5)

BinaryTables <- list()

for (j in 1:length(LineagesAll)) {
  EpimutLineages <- lapply(LineagesAll[[j]], function(x) GetEpimutations(x, Control))
  OverallTab <- do.call(cbind, lapply(EpimutLineages, function(x) x[, ncol(x)]))
  rownames(OverallTab) <- EpimutLineages[[1]][, 5]
  colnames(OverallTab) <- LineagesAll[[j]]
  
  BinaryTab <- NULL
  for (thr in cutoffs) {
    RangeTab <- matrix(0, nrow = nrow(OverallTab), ncol = ncol(OverallTab))
    for (i in 1:ncol(OverallTab)) {
      a <- which(OverallTab[, i] > thr)
      b <- which(OverallTab[, i] < -thr)
      RangeTab[a, i] <- 1
      RangeTab[b, i] <- -1
    }
    colnames(RangeTab) <- paste0(colnames(OverallTab), "_gt_", thr)
    BinaryTab <- cbind(BinaryTab, RangeTab)
  }
  BinaryTables[[j]] <- BinaryTab
}
################################################################################
#Calculate the Jaccard index
BinaryL7T4 <- as.data.frame(BinaryTab)
r1 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_1 ==  BinaryL7T4$L7T4.windows.bed_norm_gt_1 & BinaryL7T4$L7T4_R2.windows.bed_norm_gt_1 != 0)
u1 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_1 != 0 | BinaryL7T4$L7T4.windows.bed_norm_gt_1 != 0)
j1 <- r1/u1

r2 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_2 ==  BinaryL7T4$L7T4.windows.bed_norm_gt_2 & BinaryL7T4$L7T4_R2.windows.bed_norm_gt_2 != 0)
u2 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_2 != 0 | BinaryL7T4$L7T4.windows.bed_norm_gt_2 != 0)
j2 <- r2/u2


r3 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_3 ==  BinaryL7T4$L7T4.windows.bed_norm_gt_3 & BinaryL7T4$L7T4_R2.windows.bed_norm_gt_3 != 0)
u3 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_3 != 0 | BinaryL7T4$L7T4.windows.bed_norm_gt_3 != 0)
j3 <- r3/u3

r4 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_4 ==  BinaryL7T4$L7T4.windows.bed_norm_gt_4 & BinaryL7T4$L7T4_R2.windows.bed_norm_gt_4 != 0)
u4 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_4 != 0 | BinaryL7T4$L7T4.windows.bed_norm_gt_4 != 0)
j4 <- r4/u4

r5 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_5 ==  BinaryL7T4$L7T4.windows.bed_norm_gt_5 & BinaryL7T4$L7T4_R2.windows.bed_norm_gt_5 != 0)
u5 <- sum(BinaryL7T4$L7T4_R2.windows.bed_norm_gt_5 != 0 | BinaryL7T4$L7T4.windows.bed_norm_gt_5 != 0)
j5 <- r5/u5


c(j1, j2, j3, j4, j5)
c(r1, r2, r3, r4, r5)



cut_db <-read.csv("cut_off.csv")

db1 <- cut_db %>%
  group_by(Treshold) %>%
  summarise(mean_jaccard = mean(Jaccard),
            sd_jaccard = sd(Jaccard))

jaccard_plot <- ggplot(data = db1, aes(x = Treshold, y = mean_jaccard))+
  geom_errorbar(aes(ymin=mean_jaccard-sd_jaccard, ymax=mean_jaccard+sd_jaccard), width=.3)+
  geom_point(size = 3)+
  labs(x="Z score treshold", y="Jaccard index")+
  theme_minimal()+
  theme(text = element_text(size = 15))

#create the supplemetary figure 1  
sup_fig1 <- ggdraw() +
  draw_plot(jaccard_plot, x = 0, y = 0, width = 0.45, height = 0.95) +
  draw_plot(method.plot, x = 0.5, y = 0, width = 0.45, height = 0.95) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0, 0.5), y = c(0.96, 0.96), size = 20)

ggsave(filename = "sup_fig1.svg", plot = sup_fig1,
       width = 20, height = 15, units = "cm")
 
#################################################################################
###################################################################################
sort_by_timepoint <- function(x) {
  timepoint <- as.numeric(sub(".*T([0-9]+).*", "\\1", x))
  x[order(timepoint)]}

L1 <- grep("^L1.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L1 <- L1[!grepl("_R2", L1)]
L1 <- sort_by_timepoint(L1)

L2 <- grep("^L2.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L2 <- L2[!grepl("_R2", L2)]
L2 <- sort_by_timepoint(L2)

L6 <- grep("^L6.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L6 <- L6[!grepl("_R2", L6)]
L6 <- sort_by_timepoint(L6)

L7 <- grep("^L7.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L7 <- L7[!grepl("_R2", L7)]
L7 <- sort_by_timepoint(L7)

# Combine into list of all lineages

LineagesAll <- list(L1, L2, L6, L7) 
# Path to the control (ancestor) count file
Control <- grep("^ANC.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
# Empty list to store pairwise results
PairsAllOut <- list()
BinaryTables <- list()   # initialize once

#Loop over each lineage and calculate epimutations
for (j in 1:length(LineagesAll)) {
  EpimutLineages <- list()
  # For each sample in the lineage, compute epimutations
  for (i in 1:length(LineagesAll[[j]])) {
    EpimutLineages[[i]] <- GetEpimutations(LineagesAll[[j]][i], Control)
  }
  # Collect residuals for all replicates in the lineage
  OverallTab <- c()
  for (i in 1:length(LineagesAll[[j]])) {
    OverallTab <- cbind(OverallTab, EpimutLineages[[i]][,ncol(EpimutLineages[[i]])])
  }
  # Set row and column names for the data
  row.names(OverallTab) <- EpimutLineages[[1]][,5]
  colnames(OverallTab) <- LineagesAll[[j]]
  
  # Convert to binary presence/absence of epimutation:
  # 1: positive residual between 1 and 2
  # -1: negative residual between -1 and -2
  # 0: otherwise
  BinaryTab <- OverallTab
  for (i in 1:ncol(OverallTab)) {
    a <- which(OverallTab[, i] > 3.5)
    b <- which(OverallTab[, i] < -3.5)
    BinaryTab[a,i] <- 1
    BinaryTab[b,i] <- -1
    BinaryTab[-c(a,b),i] <- 0
  }
  
  # Assign generation numbers for each sample
  Gen <- as.numeric(sub(".*T([0-9]+)\\..*", "\\1", colnames(BinaryTab)))
  
  # Compute pairwise Spearman correlations between samples
  OutPairs <- c()
  for (i in 1:ncol(BinaryTab)) {
    for (k in 1:ncol(BinaryTab)) {
      if (!(i == k)) {
        COR <- cor.test(BinaryTab[,i], BinaryTab[,k], method="spearman")$estimate
        OutPairs <- rbind(OutPairs, c(Gen[i], Gen[k], COR))
      }
    }
  }
  
  # Store output for this lineage
  PairsAllOut[[j]] <- OutPairs
}
#################################################################################
#Plot correlation vs generations of separation
# Combine all pairwise correlation results
allpoints <- c()
for (i in 1:2) {
  allpoints <- rbind(allpoints, cbind(PairsAllOut[[i]][,2] - PairsAllOut[[i]][,1], 
                                      PairsAllOut[[i]][,3]))
}

# Keep only positive generation differences
allpoints <- allpoints[which(allpoints[,1] > 0),]

df <- data.frame(
  generations = allpoints[,1],
  correlation = allpoints[,2]
)

lmfit <- lm(correlation ~ generations, data = df)
pred <- predict(lmfit, newdata = df, se.fit = TRUE)

plotdf <- data.frame(
  generations = df$generations,
  fit   = pred$fit,
  upper = pred$fit + pred$se.fit,
  lower = pred$fit - pred$se.fit
)

# GGplot
divergence_control <- ggplot(df, aes(x = generations, y = correlation)) +
  geom_point() +
  geom_line(data = plotdf, aes(y = fit), color = "#298c8c", size = 1) +
  geom_ribbon(data = plotdf,
              aes(x = generations, ymin = lower, ymax = upper),
              fill = "#298c8c", alpha = 0.2, inherit.aes = FALSE) +
  labs(
    x = "Transfer of separation",
    y = "Correlation"
  ) +
  theme_minimal()+
  theme(text = element_text(size = 30))

ggsave(filename = "divergence_ethanol.svg", plot = divergence_control,
       width = 15, height = 10, units = "cm")



# Save plot
dev.copy(pdf, "EthanolEpimutationscorrelations.pdf")
dev.off()

# Save R workspace
allpoints_allchr <- allpoints
save.image("allchrdata.Rdata")
################################################################################
#STOP AT THE BINARY TABLE FOR SURVIVAL ANALYSIS
BinaryTablesAll <- list()
ZScoreTablesAll <- list()
ResidualTablesAll <- list()

for (j in 1:length(LineagesAll)) {
  EpimutLineages <- list()
  for (i in 1:length(LineagesAll[[j]])) {
    EpimutLineages[[i]] <- GetEpimutations(LineagesAll[[j]][i], Control)
  }
  
  # Combine the last column of each result into one table
  OverallTab <- c()
  for (i in 1:length(LineagesAll[[j]])) {
    OverallTab <- cbind(OverallTab, EpimutLineages[[i]][, ncol(EpimutLineages[[i]])])
  }
  row.names(OverallTab) <- EpimutLineages[[1]][, 1]
  colnames(OverallTab) <- LineagesAll[[j]]
  
  ZScoreTablesAll[[j]] <- OverallTab
  
  ResidualTab <- c()
  for (i in 1:length(LineagesAll[[j]])) {
    ResidualTab <- cbind(ResidualTab, EpimutLineages[[i]][, "Residual"])
  }
  row.names(ResidualTab) <- EpimutLineages[[1]][, 1]
  colnames(ResidualTab) <- LineagesAll[[j]]
  
  ResidualTablesAll[[j]] <- ResidualTab
  
  BinaryTab <- OverallTab
  for (i in 1:ncol(OverallTab)) {
    a <- which(OverallTab[, i] > 3.5)
    b <- which(OverallTab[, i] < -3.5)
    BinaryTab[a, i] <- 1
    BinaryTab[b, i] <- -1
    BinaryTab[-c(a, b), i] <- 0
  }
  BinaryTablesAll[[j]] <- BinaryTab
}

save(BinaryTablesAll, file = "BinaryTablesAll.RData")
##############################################################################
load("/home/bioc1877-local/Documents/MA_pombe/sRNApombe/BinaryTablesAll.RData")
mat_epi <- do.call(cbind, BinaryTablesAll)
zscore.df <- do.call(cbind, ZScoreTablesAll)
zscore.df <- as.data.frame(zscore.df)

residual.df <- do.call(cbind, ResidualTablesAll)
residual.df <- as.data.frame(residual.df)

colnames(mat_epi) <- gsub("\\.windows\\.bed_norm", "", colnames(mat_epi))
filtered_mat <- mat_epi[apply(mat_epi, 1, function(row) any(row != 0)), ]
filtered_mat.df <- as.data.frame(filtered_mat)

epi_windows.df <- filtered_mat.df %>%
  pivot_longer(
    cols = all_of(samp.names),
    names_to = "sample",
    values_to = "epi_value")

epi_windows.df %>% filter(epi_value == 1)%>%
  summarise(n())

hist(residual.df$L1T7.windows.bed_norm, breaks = 40)
qqnorm(residual.df$L6T9.windows.bed_norm)
qqline(residual.df$L6T9.windows.bed_norm)
##############################################################################
##############################################################################
#Load the csv of the normalized cds counts
cds.norm_counts <- read_csv("cds.normalized_counts_windows.csv")

# Identify all sample columns (those starting with "count.")
sample_cols <- setdiff(names(cds.norm_counts), c("chrom", "start", "end", "name"))
# Create a list of sample-specific data.frames
sample_dfs <- lapply(sample_cols, function(col) {
  data.frame(
    chrom = cds.norm_counts$chrom,
    start = cds.norm_counts$start,
    end   = cds.norm_counts$end,
    normalized_count = cds.norm_counts[[col]]  # consistent column name
  )
})

# Name each element in the list according to the original sample column
names(sample_dfs) <- sample_cols
samples <- names(sample_dfs)

for(nm in names(sample_dfs)) {
  # Create a clean variable name by removing "count." and adding "_norm"
  df_name <- nm
  # Assign the data.frame to the environment
  assign(df_name, sample_dfs[[nm]], envir = .GlobalEnv)
}

######### extract the replicates
# Keep only their names (as character vectors)
sort_by_timepoint <- function(x) {
  timepoint <- as.numeric(sub(".*T([0-9]+).*", "\\1", x))
  x[order(timepoint)]}

L1 <- grep("^L1.*", samples, value = TRUE)
L1 <- L1[!grepl("_R2", L1)]
L1 <- sort_by_timepoint(L1)

L2 <- grep("^L2.*", samples, value = TRUE)
L2 <- L2[!grepl("_R2", L2)]
L2 <- sort_by_timepoint(L2)

L6 <- grep("^L6.*", samples, value = TRUE)
L6 <- L6[!grepl("_R2", L6)]
L6 <- sort_by_timepoint(L6)

L7 <- grep("^L7.*", samples, value = TRUE)
L7 <- L7[!grepl("_R2", L7)]
L7 <- sort_by_timepoint(L7)

# Combine into list of all lineages
LineagesAll <- list(L6, L7) 


# Path to the control (ancestor) count file
Control <- grep("^ANC", samples, value = TRUE)

# Empty list to store pairwise results
PairsAllOut <- list()
BinaryTables <- list()   # initialize once

#Loop over each lineage and calculate epimutations
for (j in 1:length(LineagesAll)) {
  EpimutLineages <- list()
  # For each sample in the lineage, compute epimutations
  for (i in 1:length(LineagesAll[[j]])) {
    EpimutLineages[[i]] <- GetEpimutations(LineagesAll[[j]][i], Control)
  }
  # Collect residuals for all replicates in the lineage
  OverallTab <- c()
  for (i in 1:length(LineagesAll[[j]])) {
    OverallTab <- cbind(OverallTab, EpimutLineages[[i]][,ncol(EpimutLineages[[i]])])
  }
  # Set row and column names for the data
  row.names(OverallTab) <- EpimutLineages[[1]][,5]
  colnames(OverallTab) <- LineagesAll[[j]]
  
  # Convert to binary presence/absence of epimutation:
  # 1: positive residual between 1 and 2
  # -1: negative residual between -1 and -2
  # 0: otherwise
  BinaryTab <- OverallTab
  for (i in 1:ncol(OverallTab)) {
    a <- which(OverallTab[, i] > 3)
    b <- which(OverallTab[, i] < -3)
    BinaryTab[a,i] <- 1
    BinaryTab[b,i] <- -1
    BinaryTab[-c(a,b),i] <- 0
  }
  
  # Assign generation numbers for each sample
  Gen <- as.numeric(sub(".*T([0-9]+)*", "\\1", colnames(BinaryTab)))
  
  # Compute pairwise Spearman correlations between samples
  OutPairs <- c()
  for (i in 1:ncol(BinaryTab)) {
    for (k in 1:ncol(BinaryTab)) {
      if (!(i == k)) {
        COR <- cor.test(BinaryTab[,i], BinaryTab[,k], method="spearman")$estimate
        OutPairs <- rbind(OutPairs, c(Gen[i], Gen[k], COR))
      }
    }
  }
  
  # Store output for this lineage
  PairsAllOut[[j]] <- OutPairs
}
#################################################################################
#Plot correlation vs generations of separation
# Combine all pairwise correlation results
allpoints <- c()
for (i in 1:2) {
  allpoints <- rbind(allpoints, cbind(PairsAllOut[[i]][,2] - PairsAllOut[[i]][,1], 
                                      PairsAllOut[[i]][,3]))
}

# Keep only positive generation differences
allpoints <- allpoints[which(allpoints[,1] > 0),]

# Fit linear model: correlation ~ generations of separation
lmfit <- lm(allpoints[,2] ~ allpoints[,1])
summary(lmfit)
# Plot results
plot(allpoints[,1], allpoints[,2], pch=18, col="black",
     ylab="correlation_epimutations",
     xlab="generations of separation")
abline(lmfit)

##############################################################################
# Do the correlation analysis
cds.norm_counts <- read_csv("cds.normalized_counts_windows.csv")
sRNA.norm_counts <- as.data.frame(cds.norm_counts)
colnames(sRNA.norm_counts)[4] <- "genID"
sRNA.norm_counts <- sRNA.norm_counts %>% dplyr::select(colnames(norm.counts))

sRNA.norm_counts <- sRNA.norm_counts[order(row.names(sRNA.norm_counts)), ]
sRNA_summary <- sRNA.norm_counts %>%
  group_by(genID) %>%         
  summarise(across(where(is.numeric), mean, na.rm = TRUE))


sRNA_long_summary <- sRNA_summary %>%
  pivot_longer(
    cols = -c(genID),  # keep first three columns
    names_to = "sample",           # new column for sample names
    values_to = "counts"            # new column for the values
  )


norm.counts <- read.csv("/home/bioc1877-local/Documents/MA_pombe/mRNApombe/normalized_counts.csv")
colnames(norm.counts) <- gsub("\\.bam", "", colnames(norm.counts)) 
colnames(norm.counts) <- gsub("_", "", colnames(norm.counts))
colnames(norm.counts)[2] <- "ANC"

mRNA.norm_counts <- norm.counts[order(row.names(norm.counts)), ]
mRNA_summary <- mRNA.norm_counts %>%
  group_by(genID) %>%         
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

mRNA_long_summary <- norm.counts %>%
  pivot_longer(
    cols = -c(genID),  # keep first three columns
    names_to = "sample",           # new column for sample names
    values_to = "counts"            # new column for the values
  )


merged_RNA <- inner_join(sRNA_long_summary, mRNA_long_summary, by = c("genID", "sample"))
merged_RNA <- as.data.frame(merged_RNA)

cor.test(merged_RNA$counts.x, merged_RNA$counts.y, method="pearson")
         
ggplot(merged_RNA, aes(x = log10(counts.x+1), y = log10(counts.y+1))) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  stat_cor(method = "pearson")+
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  labs(
    x = "Countsx",
    y = "Countsy",
  ) +
  theme_minimal(base_size = 14)

ALL.EPI <- readLines("all_epi.bed_geneid")
merged_epimut <- merged_RNA%>% dplyr::filter(genID%in%ALL.EPI)

ggplot(merged_epimut, aes(x = log10(counts.x+1), y = log10(counts.y+1))) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  stat_cor(method = "pearson")+
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  labs(
    x = "Countsx",
    y = "Countsy",
  ) +
  theme_minimal(base_size = 14)


load("/home/bioc1877-local/Documents/MA_pombe/mRNApombe/mRNABinaryTablesAll.RData")
head(BinaryTablesAll)
mat_mepi <- do.call(cbind, BinaryTablesAll)
head(mat_mepi)
filtered_mmat <- mat_mepi[apply(mat_mepi, 1, function(row) any(row != 0)), ]
mepimutations <- noquote(rownames(filtered_mmat))
merged_mepimut <- merged_RNA%>% dplyr::filter(genID%in%mepimutations)


ggplot(merged_mepimut, aes(x = log10(counts.x+1), y = log10(counts.y+1))) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  stat_cor(method = "pearson")+
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  labs(
    x = "Countsx",
    y = "Countsy",
  ) +
  theme_minimal(base_size = 14)

###############################################################################
figure_db$highlight <- "Normal"
figure_db$highlight[figure_db$Zscore > 3]  <- "Z > 3"
figure_db$highlight[figure_db$Zscore < -3] <- "Z < -3"
figure_db$highlight <- factor(figure_db$highlight,
                              levels = c("Z < -3", "Normal", "Z > 3"))


LX_TX.df <- norm_counts %>% select(count.ANC.windows.bed, count.L1T1.windows.bed)
LX_TX.df$zscore <- zscore.df$L1T1.windows.bed_norm

ggplot(LX_TX.df, aes(x=log2(count.ANC.windows.bed+1), y=log2(count.L1T1.windows.bed+1)))+
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "black")

method.plot <- ggplot(figure_db,
                      aes(x = log2(normalized_count.cont + 1),
                          y = log2(normalized_count.test + 1))) +
  geom_point(aes(color = highlight), size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  scale_color_manual(values = c(
    "Z < -3" = "orange",
    #"Normal" = "grey70",
    "Z > 3"  = "red"
  )) +
  theme_minimal() +
  labs(
    x = "log2(normalized control\n siRNAs counts + 1)",
    y = "log2(normalized L1T1\n siRNAs counts + 1)",
    color = "Z-score",
  )+
  theme(text = element_text(size = 15))


#make pdf of  residual plots:
cols <- names(residual.df)

# -------------------------
# Histograms in a grid PDF
# -------------------------
pdf("residual_histograms_grid.pdf", width = 12, height = 10)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (col in cols) {
  x <- residual.df[[col]]
  x <- x[!is.na(x)]
  
  hist(x, breaks = 40,
       main = col,
       xlab = "",
       col = "grey")
}

dev.off()

# -------------------------
# QQ plots in a grid PDF
# -------------------------
pdf("residual_qqplots_grid.pdf", width = 12, height = 10)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

for (col in cols) {
  x <- residual.df[[col]]
  x <- x[!is.na(x)]
  
  qqnorm(x, main = col)
  qqline(x, col = "red")
}

dev.off()