setwd("/home/bioc1877-local/Documents/MA_pombe/sRNApombe")
#Get the epimutation based on the log2FC
#create the function
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
  # Calculate log2FC: log2(test counts +1)/ log2(control counts+1)
  log2FC <- log2((merged$normalized_count.test + 1)/(merged$normalized_count.cont + 1))
  # Add results back
  merged$log2FC <- log2FC
  return(merged)
}
################################################################################
#Load the csv of the normalized counts
norm_counts <- read.csv("normalized_counts_windows.csv")
#########################################################################
#This is for centromere and telomeres
#heterochromatin <- read.table("~/Documents/MA_pombe/reference/centromeresandtelomeres.bed")
#norm_counts <- norm_counts %>%
#  rowwise() %>%
#  filter(
#    any(
#      chrom == heterochromatin$V1 &
#        start < heterochromatin$V3 &
#        end > heterochromatin$V2)) %>%
#  ungroup()
##########################################################################
#norm_counts <- norm_counts %>% filter(count.ANC.windows.bed < 5791)
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
################################################################################
sort_by_timepoint <- function(x) {
  timepoint <- as.numeric(sub(".*T([0-9]+).*", "\\1", x))
  x[order(timepoint)]}

all_norms <- ls(pattern = "_norm$")

#L1 and L2 =control
L1 <- grep("^L1.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L1 <- L1[!grepl("_R2", L1)]
L1 <- sort_by_timepoint(L1)

L2 <- grep("^L2.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L2 <- L2[!grepl("_R2", L2)]
L2 <- sort_by_timepoint(L2)
#L6 and L7 = 4% ethanol
L6 <- grep("^L6.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L6 <- L6[!grepl("_R2", L6)]
L6 <- sort_by_timepoint(L6)

L7 <- grep("^L7.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
L7 <- L7[!grepl("_R2", L7)]
L7 <- sort_by_timepoint(L7)

# Combine into list of all lineages

LineagesAll <- list(L1, L2) 
# Path to the control (ancestor) count file
Control <- grep("^ANC.*\\.windows\\.bed_norm$", all_norms, value = TRUE)
#STOP AT THE BINARY TABLE FOR SURVIVAL ANALYSIS
BinaryTablesAll <- list()
log2FCTablesAll <- list()
PairsAllOut <- list()

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
  #esto rea lo equivalnete a tener una lista de los Zscores
  log2FCTablesAll[[j]] <- OverallTab
  
  BinaryTab <- OverallTab 
  for (i in 1:ncol(OverallTab)) {
    x <- OverallTab[, i]
    Zmod <- (x - median(x)) / mad(x)
    
    BinaryTab[, i] <- 0
    BinaryTab[Zmod > 3, i] <- 1
    BinaryTab[Zmod < -3, i] <- -1
  }
  BinaryTablesAll[[j]] <- BinaryTab
  # Assign generation numbers for each sample
  Gen <- as.numeric(sub(".*T([0-9]+)\\..*", "\\1", colnames(BinaryTab)))
  
  # Compute pairwise Spearman correlations between samples
  OutPairs <- c()
  for (i in 1:ncol(BinaryTab)) {
    for (k in 1:ncol(BinaryTab)) {
      if (!(i == k)) {
        COR <- cor.test(BinaryTab[,i], BinaryTab[,k], method="spearman")$estimate
        OutPairs <- rbind(OutPairs, c(Gen[i], Gen[k], COR))}}
  }
  
  # Store output for this lineage
  PairsAllOut[[j]] <- OutPairs}

#save(BinaryTablesAll, file = "BinaryTablesAll.RData")
################################################################################
#Examine the Binary tables
load("BinaryTablesAll.RData")
mat_epi <- do.call(cbind, BinaryTablesAll)

colnames(mat_epi) <- gsub("\\.windows\\.bed_norm", "", colnames(mat_epi))
filtered_mat <- mat_epi[apply(mat_epi, 1, function(row) any(row != 0)), ]
filtered_mat.df <- as.data.frame(filtered_mat)

epi_windows.df <- filtered_mat.df %>%
  pivot_longer(
    cols = all_of(colnames(filtered_mat.df)),
    names_to = "sample",
    values_to = "epi_value")

epi_windows.df %>% filter(epi_value == -1)%>%
  summarise(n())
################################################################################
#Examine divergence plots
# Combine all pairwise correlation results
allpoints <- c()
for (i in 1:2) {
  allpoints <- rbind(allpoints, cbind(PairsAllOut[[i]][,2] - PairsAllOut[[i]][,1], 
                                      PairsAllOut[[i]][,3]))}

# Keep only positive generation differences
allpoints <- allpoints[which(allpoints[,1] > 0),]

df <- data.frame(
  generations = allpoints[,1],
  correlation = allpoints[,2])

lmfit <- lm(correlation ~ generations, data = df)
anova(lmfit)
pred <- predict(lmfit, newdata = df, se.fit = TRUE)

plotdf <- data.frame(
  generations = df$generations,
  fit   = pred$fit,
  upper = pred$fit + pred$se.fit,
  lower = pred$fit - pred$se.fit)

# GGplot
divergence_ethanol <- ggplot(df, aes(x = generations, y = correlation)) +
  geom_point() +
  geom_line(data = plotdf, aes(y = fit), color = "#ffbb6f", size = 1) +
  geom_ribbon(data = plotdf,
              aes(x = generations, ymin = lower, ymax = upper),
              fill = "#ffbb6f", alpha = 0.2, inherit.aes = FALSE) +
  #annotate("text", x = 8.2, y = 0.61, label = "p < 0.01", size =7)+
  labs(x = "Transfer of separation", y = "Correlation") +
  theme_minimal()+
  theme(text = element_text(size = 30))

ggsave(filename = "divergence_ethanol.svg", plot = divergence_ethanol,
       width = 15, height = 10, units = "cm")

divergence_control <- ggplot(df, aes(x = generations, y = correlation)) +
  geom_point() +
  geom_line(data = plotdf, aes(y = fit), color = "#298c8c", size = 1) +
  geom_ribbon(data = plotdf,
              aes(x = generations, ymin = lower, ymax = upper),
              fill = "#298c8c", alpha = 0.2, inherit.aes = FALSE) +
  labs(x = "Transfer of separation", y = "Correlation") +
  theme_minimal()+
  theme(text = element_text(size = 30))

ggsave(filename = "divergence_control_het.svg", plot = divergence_control,
       width = 15, height = 10, units = "cm")

    ################################################
     #Count the number of epimutations/transfers etc
    ################################################
load("/home/bioc1877-local/Documents/MA_pombe/sRNApombe/BinaryTablesAll.RData")
mat_epi <- do.call(cbind, BinaryTablesAll)
colnames(mat_epi) <- gsub("\\.windows\\.bed_norm", "", colnames(mat_epi))
filtered_mat <- mat_epi[apply(mat_epi, 1, function(row) any(row != 0)), ]
# Function to count epimutations (stable runs of the same value that last 
#>= 2 transfers)
count_stable_by_value <- function(x) {
  r <- rle(x)
  vals <- r$values
  lens <- r$lengths
  
  # true for a zero-run only if a nonzero run appeared before it
  zero_allowed <- cumsum(vals != 0) > 0
  
  c(
    epimutation_0  = sum(lens >= 2 & vals == 0 & zero_allowed),
    epimutation_1  = sum(lens >= 2 & vals == 1),
    epimutation_m1 = sum(lens >= 2 & vals == -1)
  )
}
#Function to count transitions
count_transitions <- function(x) {
  sum(diff(x) != 0)}

# Define which transfers belong to each line
line_transfers <- list(
  L1 = paste0("L1T", 1:10),
  L2 = paste0("L2T", 1:8),
  L6 = paste0("L6T", 1:10),
  L7 = paste0("L7T", 1:10))
sort_by_timepoint(line_transfers)

# Build the final dataframe
result_list <- lapply(names(line_transfers), function(line) {
  cols <- intersect(line_transfers[[line]], colnames(filtered_mat))
  if (length(cols) == 0) return(NULL)
  
  x <- filtered_mat[, cols, drop = FALSE]
  
  # keep only rows where at least one value in that line is nonzero
  keep <- rowSums(x != 0) > 0
  x <- x[keep, , drop = FALSE]
  
  if (nrow(x) == 0) return(NULL)
  
  stable_counts <- t(apply(x, 1, count_stable_by_value))
  
  data.frame(
    window = rownames(x),
    number_of_epimutation = rowSums(stable_counts),
    transitions = apply(x, 1, count_transitions),
    stable_counts,
    line = line,
    row.names = NULL
  )
})

result_df <- do.call(rbind, result_list)
#add a column indicating the condition
result_df <- result_df %>%
  dplyr::mutate(
    condition = dplyr::case_when(
      line %in% c("L1", "L2") ~ "control",
      line %in% c("L6", "L7") ~ "ethanol",
      TRUE ~ NA_character_
    ),
    number_transfers = dplyr::case_when(
      line %in% c("L1", "L2", "L7") ~ 10,
      line == "L6" ~ 9,
      TRUE ~ NA_real_ ))

#make summary with the total numer of transitions and epimutations. 
result_df %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    transition = sum(transitions),
    epimutation = sum(number_of_epimutation),
    epi_0 = sum(epimutation_0),
    epi_1 = sum(epimutation_1),
    epi_m1 =sum(epimutation_m1),
    .groups = "drop")
#model to test if there is difference in the number of epimutations

model_epimutations <- glmer( number_of_epimutation ~ condition + offset(log(number_transfers)) + (1 | line),
  data = result_df, family = poisson())

summary(model_epimutations)
anova(model_epimutations)
check_overdispersion(model_epimutations)
#model to test difference in the number of transitions
model_transitions <- glmer( transitions ~ condition + offset(log(number_transfers)) + (1 | line),
                            data = result_df, family = poisson)

summary(model_transitions)
check_overdispersion(model_transitions)

transitions.plot <- ggplot(result_df, aes(x = transitions, fill = condition, colour = condition)) +
  geom_bar(position = "identity", alpha = 0.3)+
  scale_color_manual(values=c("#298c8c","#ffbb6f"))+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  theme_minimal(base_size = 30) +
  labs(x = "Number of transitions", y = "Frequency")+
  theme(legend.position = "none")

ggsave(filename = "transitions.svg", plot = transitions.plot,
       width = 15, height = 15, units = "cm")

epimutations.plot <- ggplot(result_df, aes(x = number_of_epimutation, fill = condition, colour = condition)) +
  geom_bar(position = "identity", alpha = 0.3)+
  scale_color_manual(values=c("#298c8c","#ffbb6f"))+
  scale_fill_manual(values=c("#298c8c","#ffbb6f"))+
  theme_minimal(base_size = 30) +
  labs(x = "Number of epimutations", y = "Frequency")+
  theme(legend.position = "none")

ggsave(filename = "epimutations.svg", plot = epimutations.plot,
       width = 15, height = 15, units = "cm")

save(result_df, file = "epimutation_number_results.RData")
##count the initial amount of 0's to see if in control it takes on average more time 
#for the first epigenetic change to ocurr

count_initial_zeros <- function(x) {
  r <- rle(x)
  vals <- r$values
  lens <- r$lengths
  
  initial_zero <- cumsum(vals != 0) == 0
  
  sum(lens[initial_zero & vals == 0])
}

result_zerolist <- lapply(names(line_transfers), function(line) {
  cols <- intersect(line_transfers[[line]], colnames(filtered_mat))
  if (length(cols) == 0) return(NULL)
  
  x <- filtered_mat[, cols, drop = FALSE]
  
  # keep only rows where at least one value in that line is nonzero
  keep <- rowSums(x != 0) > 0
  x <- x[keep, , drop = FALSE]
  
  if (nrow(x) == 0) return(NULL)
  
  
  data.frame(
    window = rownames(x),
    initial_zero = apply(x, 1, count_initial_zeros),
    line = line,
    row.names = NULL
  )
})

result_zerodf <- do.call(rbind, result_zerolist)
#add a column indicating the condition
result_zerodf <- result_zerodf %>%
  dplyr::mutate(
    condition = dplyr::case_when(
      line %in% c("L1", "L2") ~ "control",
      line %in% c("L6", "L7") ~ "ethanol"))

result_zerodf %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    initial_zero_length = median(initial_zero),
    .groups = "drop")

model_initialzero <- glm(initial_zero ~ condition, data = result_zerodf, family = poisson())
summary(model_initialzero)

ggplot(result_zerodf, aes(x = factor(initial_zero), fill = condition, colour = condition)) +
  geom_bar(position = "identity", alpha = 0.3) +
  scale_color_manual(values=c("#298c8c","#ffbb6f")) +
  scale_fill_manual(values=c("#298c8c","#ffbb6f")) +
  theme_minimal(base_size = 30) +
  labs(x = "Count", y = "Initial zeros")
#there is no difference bewteen the length of inital zeros. 

    ##########################################
    #Where in the genome the epimutations are?
    ##########################################

epimutation_control <- result_df %>% filter(condition == "control") %>% 
  separate(window, into = c("chrom", "start", "end"), sep = ":")
#convert into genome ranges
control_ranges <- makeGRangesFromDataFrame(
  epimutation_control,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE)

epimutation_ethanol <- result_df %>% filter(condition == "ethanol") %>%
 separate(window, into = c("chrom", "start", "end"), sep = ":")
#convert into genome ranges
ethanol_ranges <- makeGRangesFromDataFrame(
  epimutation_ethanol,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE)

gtf.file <- import.gff("~/Documents/MA_pombe/reference/Schizosaccharomyces_pombe_all_chromosomes.gtf")
 
epi_genes_control <- subsetByOverlaps(gtf.file, control_ranges, minoverlap = 50)
epi_genes_control.df <- as.data.frame(epi_genes_control) %>% filter(type == "transcript")
control_geneID <- unique(epi_genes_control.df$gene_id)
  
epi_genes_ethanol <- subsetByOverlaps(gtf.file, ethanol_ranges, minoverlap = 50)
epi_genes_ethanol.df <- as.data.frame(epi_genes_ethanol) %>% filter(type == "transcript")
ethanol_geneID <- unique(epi_genes_ethanol.df$gene_id)              

save(control_geneID, ethanol_geneID, file = "epimutations_GeneID.RData")

##get the coordinates of the epimutations
#First the control
control_epimutations <- result_df %>% dplyr::filter(number_of_epimutation != 0, condition == "control")
control_epimutations[duplicated(control_epimutations$window), ]
sum(duplicated(control_epimutations$window))
control_epimutations <- control_epimutations[!duplicated(control_epimutations$window, fromLast = TRUE), ]
control.epicoord <- control_epimutations %>%
  separate(window, into = c("chrom", "start", "end"), sep = ":", convert = TRUE)


ethanol_epimutations <- result_df %>% dplyr::filter(number_of_epimutation != 0, condition == "ethanol")
ethanol_epimutations[duplicated(ethanol_epimutations$window), ]
sum(duplicated(ethanol_epimutations$window))
ethanol_epimutations <- ethanol_epimutations[!duplicated(ethanol_epimutations$window, fromLast = TRUE), ]
ethanol.epicoord <- ethanol_epimutations %>%
  separate(window, into = c("chrom", "start", "end"), sep = ":", convert = TRUE)

all.epimutations <- result_df[!duplicated(result_df$window, fromLast = TRUE), ]
all.epicoord <- all.epimutations %>%
  separate(window, into = c("chrom", "start", "end"), sep = ":", convert = TRUE)
  
save(control.epicoord, ethanol.epicoord, all.epicoord, file = "epimutations_coord.RData")


############################################################################################
#create the plots
filtered_mat <- as.data.frame(filtered_mat)
filtered_mat.df_L1 <- filtered_mat %>% dplyr::select("L1T1")
filtered_mat.df_L1$coord <- rownames(filtered_mat.df_L1)
rownames(filtered_mat.df_L1) <- NULL

norm_counts$coord <-  paste(norm_counts$chrom, norm_counts$start, norm_counts$end, sep=":")


# Match binary classification by coordinate
norm_counts$L1T1_highlight <- filtered_mat.df_L1$L1T1[
  match(norm_counts$coord, filtered_mat.df_L1$coord)]
#z mod score threshold 3
# Plot
example.plot <- ggplot(norm_counts, aes(x = log2(count.ANC.windows.bed + 1), y = log2(count.L1T1.windows.bed + 1))) +
  geom_point(size = 1, alpha = 0.4, colour = "grey") +
  geom_point(data = subset(norm_counts, L1T1_highlight != 0), size = 1, alpha = 0.9, colour = "red") +
  theme_minimal(base_size = 20)+
  labs(
    x = "log2(Anc count + 1)",
    y = "log2(L1T1 count + 1)")

ggsave(filename = "example_sup.svg", plot = example.plot,
       width = 10, height = 12, units = "cm")


      
             ###############################
             #Gene expression vs siRNA count
            ################################

mRNA_windows_norm_count <- read.csv("~/Documents/MA_pombe/mRNApombe/bed_files/normalized_counts_mRNAwindows.csv")

mRNA_windows_norm_count <- mRNA_windows_norm_count %>%
    rename("Anc" = "count.Anc",
           "L1_T10" = "count.L1_T10", 
           "L1_T5" = "count.L1_T5",  
           "L2_T10" = "count.L2_T10",
           "L2_T5" = "count.L2_T5",  
           "L6_T10" = "count.L6_T10", 
           "L6_T5" = "count.L6_T5",  
           "L7_T10" = "count.L7_T10", 
           "L7_T5" ="count.L7_T5")

mRNA_windows_log2FC  <- data.frame(chrom = mRNA_windows_norm_count$chrom,
                                   start = mRNA_windows_norm_count$start,
                                   end = mRNA_windows_norm_count$end,
                                   L1_T10 = log2((mRNA_windows_norm_count$L1_T10+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L1_T5 = log2((mRNA_windows_norm_count$L1_T5+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L2_T10 = log2((mRNA_windows_norm_count$L2_T10+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L2_T5 = log2((mRNA_windows_norm_count$L2_T5+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L6_T10 = log2((mRNA_windows_norm_count$L6_T10+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L6_T5 = log2((mRNA_windows_norm_count$L6_T5+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L7_T10 = log2((mRNA_windows_norm_count$L7_T10+1)/(mRNA_windows_norm_count$Anc+1)),
                                   L7_T5 = log2((mRNA_windows_norm_count$L7_T5+1)/(mRNA_windows_norm_count$Anc+1)))

#Load the siRNA data
norm_counts <- read.csv("/home/bioc1877-local/Documents/MA_pombe/sRNApombe/normalized_counts_windows.csv")
norm_counts_fil <- norm_counts %>% dplyr::select("chrom", "start", "end", "count.ANC.windows.bed",
                                                 "count.L1T10.windows.bed", "count.L1T5.windows.bed",
                                                 "count.L2T10.windows.bed", "count.L2T5.windows.bed", 
                                                 "count.L6T10.windows.bed", "count.L6T5.windows.bed", 
                                                 "count.L7T10.windows.bed", "count.L7T5.windows.bed")

norm_counts_fil <- norm_counts_fil %>%
  rename("Anc" = "count.ANC.windows.bed",
         "L1_T10" = "count.L1T10.windows.bed", 
         "L1_T5" = "count.L1T5.windows.bed",  
         "L2_T10" = "count.L2T10.windows.bed",
         "L2_T5" = "count.L2T5.windows.bed",  
         "L6_T10" = "count.L6T10.windows.bed", 
         "L6_T5" = "count.L6T5.windows.bed",  
         "L7_T10" = "count.L7T10.windows.bed", 
         "L7_T5" = "count.L7T5.windows.bed")

siRNA_log2FC  <- data.frame(chrom = norm_counts_fil$chrom,
                                start = norm_counts_fil$start,
                                end = norm_counts_fil$end,
                                L1_T10 = log2((norm_counts_fil$L1_T10+1)/(norm_counts_fil$Anc+1)),
                                L1_T5 = log2((norm_counts_fil$L1_T5+1)/(norm_counts_fil$Anc+1)),
                                L2_T10 = log2((norm_counts_fil$L2_T10+1)/(norm_counts_fil$Anc+1)),
                                L2_T5 = log2((norm_counts_fil$L2_T5+1)/(norm_counts_fil$Anc+1)),
                                L6_T10 = log2((norm_counts_fil$L6_T10+1)/(norm_counts_fil$Anc+1)),
                                L6_T5 = log2((norm_counts_fil$L6_T5+1)/(norm_counts_fil$Anc+1)),
                                L7_T10 = log2((norm_counts_fil$L7_T10+1)/(norm_counts_fil$Anc+1)),
                                L7_T5 = log2((norm_counts_fil$L7_T5+1)/(norm_counts_fil$Anc+1)))

line.names <- c("L1_T10", "L1_T5", "L2_T10", "L2_T5", "L6_T10", "L6_T5", "L7_T10", "L7_T5")

mRNA_log2FC.long <- mRNA_windows_log2FC%>%
  pivot_longer(
    cols = all_of(line.names),
    names_to = "sample",
    values_to = "mRNAFC")

siRNA_log2FC.long <- siRNA_log2FC %>%
  pivot_longer(
    cols = all_of(line.names),
    names_to = "sample",
    values_to = "siRNAFC")

merged_norm_counts <- inner_join(mRNA_log2FC.long, siRNA_log2FC.long, 
                                 by = c("chrom", "start", "end", "sample"))

merged_counts_ranges <- makeGRangesFromDataFrame(
  merged_norm_counts,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE)

########### Get the samples with epimutations
#filter the filtered mat to have only the samples that have mRNA
head(filtered_mat)
filtered_mat <- as.data.frame(filtered_mat)
mod_mat <- filtered_mat %>% dplyr::select("L1T10", "L1T5", "L2T10", "L2T5", "L6T10",
                                   "L6T5", "L7T10", "L7T5") %>%
                            dplyr::rename("L1_T10" = "L1T10", 
                                   "L1_T5" = "L1T5",  
                                   "L2_T10" = "L2T10",
                                   "L2_T5" = "L2T5",  
                                   "L6_T10" = "L6T10", 
                                   "L6_T5" = "L6T5",  
                                   "L7_T10" = "L7T10", 
                                   "L7_T5" = "L7T5")
mod_mat$coord <- rownames(filtered_mat)

mod_mat.long <- mod_mat %>%
  pivot_longer(
    cols = all_of(line.names),
    names_to = "sample",
    values_to = "epimutations")
#add a coord column to the "merged_norm_counts"               
merged_norm_counts$coord <- paste(merged_norm_counts$chrom, merged_norm_counts$start, merged_norm_counts$end, sep=":")                          

merged_epimutations <- inner_join(merged_norm_counts, mod_mat.long, 
                                 by = c("coord", "sample"))
#filter for only the windows/samples with an epimutation
merged_epimutations <- merged_epimutations  %>% dplyr::filter(epimutations!= 0)
#make it a Grange object
epi_gr <- GRanges(
  seqnames = merged_epimutations$chrom,
  ranges = IRanges(
    start = merged_epimutations$start,
    end = merged_epimutations$end)
)

#Load transposon coordinates
transposons <- read.table("~/Documents/MA_pombe/reference/transposons.bed")
transposons <- as.data.frame(transposons)
colnames(transposons) <- c("chrom", "start", "end")

transposons_ranges <- makeGRangesFromDataFrame(
  transposons,
  seqnames.field = "chrom",
  start.field = "start",
  end.field = "end")

# Find overlaps
hits <- findOverlaps(
  epi_gr,
  #gtf.file,
  transposons_ranges,
  minoverlap = 50,
  ignore.strand = TRUE)

#add true oe false column
merged_epimutations$transposon_overlap <- FALSE
merged_epimutations$transposon_overlap[unique(queryHits(hits))] <- TRUE

merged_epimutations <- merged_epimutations %>%
  mutate(
    quadrant = case_when(
      mRNAFC >= 0 & siRNAFC >= 0 ~ "High mRNA, high siRNA",
      mRNAFC <  0 & siRNAFC >= 0 ~ "Low mRNA, high siRNA",
      mRNAFC <  0 & siRNAFC <  0 ~ "Low mRNA, low siRNA",
      mRNAFC >= 0 & siRNAFC <  0 ~ "High mRNA, low siRNA"))

epimutation.plot <- ggplot(merged_epimutations, aes(x = mRNAFC, y = siRNAFC, color = quadrant)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_point(
    data = subset(merged_epimutations, transposon_overlap == TRUE),
    aes(x = mRNAFC, y = siRNAFC),
    color = "red",
    size = 1,
    alpha = 1,
    inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "High mRNA, high siRNA" = "#1b9e77",
    "High mRNA, low siRNA" = "#7570b3",
    "Low mRNA, low siRNA" = "#d95f02",
    "Low mRNA, high siRNA" = "#e7298a"))+
  theme_minimal(base_size = 25) +
  labs(
    x = "log2FC(mRNA)",
    y = "log2FC(siRNA)",
    color = NULL,
    title = "Epimutations")

control_genes.df <- merged_epimutations %>%
  filter(grepl("L1|L2", sample)) 

control.plot <- ggplot(control_genes.df, aes(x = mRNAFC, y = siRNAFC, color = quadrant)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_point(
    data = subset(control_genes.df, transposon_overlap == TRUE),
    aes(x = mRNAFC, y = siRNAFC),
    color = "red",
    size = 1,
    alpha = 1,
    inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "High mRNA, high siRNA" = "#1b9e77",
    "High mRNA, low siRNA" = "#7570b3",
    "Low mRNA, low siRNA" = "#d95f02",
    "Low mRNA, high siRNA" = "#e7298a"))+
  theme_minimal(base_size = 30) +
  labs(
    x = "log2(mRNA)",
    y = "log2(siRNA)",
    color = "Quadrant",
    title = "Control epimutations")


ethanol_genes.df <-merged_epimutations %>%
  filter(grepl("L6|L7", sample))

ethanol.plot <- ggplot(ethanol_genes.df, aes(x = mRNAFC, y = siRNAFC, color = quadrant)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_point(
    data = subset(ethanol_genes.df, transposon_overlap == TRUE),
    aes(x = mRNAFC, y = siRNAFC),
    color = "red",
    size = 1,
    alpha = 1,
    inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c(
    "High mRNA, high siRNA" = "#1b9e77",
    "High mRNA, low siRNA" = "#7570b3",
    "Low mRNA, low siRNA" = "#d95f02",
    "Low mRNA, high siRNA" = "#e7298a"))+
  theme_minimal(base_size = 30) +
  labs(
    x = "log2(mRNA)",
    y = "log2(siRNA)",
    color = "Quadrant",
    title = "Ethanol")
      
       #######################
        #permutation analysis
      ########################

# Observed overlap for each quadrant
obs <- merged_epimutations %>%
  group_by(quadrant) %>%
  summarise(observed = mean(transposon_overlap))

# Permutations
set.seed(123)

perm <- bind_rows(lapply(1:10000, function(i) {
  merged_epimutations %>%
    group_by(sample) %>%
    mutate(quadrant_perm = sample(quadrant)) %>%
    ungroup() %>%
    group_by(quadrant_perm) %>%
    summarise(expected_perm = mean(transposon_overlap)) %>%
    mutate(permutation = i)}))

# Summary + one-sided p-value for enrichment
results <- obs %>%
  left_join(
    perm %>%
      group_by(quadrant_perm) %>%
      summarise(
        expected = mean(expected_perm),
        p_value = (sum(expected_perm >=
                         obs$observed[match(first(quadrant_perm), obs$quadrant)]) + 1) /
          (n() + 1)
      ),
    by = c("quadrant" = "quadrant_perm")
  ) %>%
  mutate(enrichment = observed / expected)

results


