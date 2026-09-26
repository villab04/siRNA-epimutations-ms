setwd("/home/bioc1877-local/Documents/MA_pombe/mRNApombe/mutants")
#BiocManager::install("RUVSeq")
library(RUVSeq)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(vegan)
library(lmerTest)

trans.counts <- read.csv("all_transcounts.csv") %>% 
  dplyr::select(!c("Chr", "Start", "End", "Strand", "Length"))

trans.counts.etoh <- trans.counts %>% select(Geneid, L41:L60)
trans.counts.YES <- trans.counts %>% select(Geneid, L61:L80)

pombe.counts <- read.csv("pombe_counts.csv")
all.pombe.counts <- rbind(trans.counts.etoh, pombe.counts)
elegans.counts <- read.csv("elegans_counts.csv")

pombe.counts <- read.csv("pombe_countsYES.csv")
all.pombe.counts <- rbind(trans.counts.YES, pombe.counts)
elegans.counts <- read.csv("elegans_countsYES.csv")


all.counts <- rbind(all.pombe.counts, elegans.counts)

#Take all of the counts that have expression above 5
filter <- apply(all.counts[, -1], 1, function(x) length(x[x>5])>=2)
filtered <- all.counts[filter,]
row.names(filtered) <- filtered[,1]
filtered <- filtered %>% dplyr::select(!Geneid)
#Name the reads that are genes and spikes
genes <- rownames(filtered)[grep("^SP|Tf2-type", rownames(filtered))]
spikes <- rownames(filtered)[grep("^WB", rownames(filtered))]
 
x <- as.factor(rep(c("WT", "Ago"), each=10))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))


colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

set1 <- RUVg(set, spikes, k=1)
pData(set1)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)

#calculate a distance matrix 
counts.df <- as.data.frame(counts(set1))
counts.df <- na.omit(counts.df)
counts.genes_ethanol <- counts.df %>% filter(rownames(counts.df) %in% genes)

#save(counts.genes_control, file = "YES_count_genes.RData")
save(counts.genes_ethanol, file = "ethanol_count_genes.RData")


#####################################################
      #####################################
        #Coefficient of variaton per gene 
      #####################################
load("ethanol_count_genes.RData")
load("YES_count_genes.RData")

WT_ethanol.df <- counts.genes_ethanol %>% dplyr::select(L41, L42, L43, L44, L45, L46, L47, L48, L49, L50)
ago_ethanol.df <- counts.genes_ethanol %>% dplyr::select(L51, L52, L53, L54, L55, L56, L57, L58, L59, L60)

WT_control.df <- counts.genes_control %>% dplyr::select(L61, L62, L63, L64, L65, L66, L67, L68, L69, L70)
ago_control.df <- counts.genes_control %>% dplyr::select(L71, L72, L73, L74, L75, L76, L77, L78, L79, L80)


counts.genes_ethanol$cv.wt <- apply(WT_ethanol.df, 1, function(x) sd(x) / mean(x))
counts.genes_ethanol$cv.ago <- apply(ago_ethanol.df, 1, function(x) sd(x) / mean(x))
counts.genes_ethanol$mean.all <- apply(counts.genes_ethanol, 1, function(x) mean(x))

counts.genes_control$cv.wt <- apply(WT_control.df, 1, function(x) sd(x) / mean(x))
counts.genes_control$cv.ago <- apply(ago_control.df, 1, function(x) sd(x) / mean(x))
counts.genes_control$mean.all <- apply(counts.genes_control, 1, function(x) mean(x))

fit.wt <- loess(counts.genes_ethanol$cv.wt ~ counts.genes_ethanol$mean.all, span = 0.3)
fit.ago <- loess(counts.genes_ethanol$cv.ago ~ counts.genes_ethanol$mean.all, span = 0.3)
pred.wt <- predict(fit.wt, se =TRUE)
pred.ago <- predict(fit.ago, se =TRUE)

pred.wt$fit <- na.omit(pred.wt$fit)
counts.genes_ethanol <- na.omit(counts.genes_ethanol)
counts.genes_ethanol$fit.wt <- pred.wt$fit
counts.genes_ethanol$lower.wt = pred.wt$fit - 1.96 * pred.wt$se.fit
counts.genes$upper.wt = pred.wt$fit + 1.96 * pred.wt$se.fit

counts.genes$fit.ago <- pred.ago$fit
counts.genes$lower.ago = pred.ago$fit - 1.96 * pred.ago$se.fit
counts.genes$upper.ago = pred.ago$fit + 1.96 * pred.ago$se.fit

counts.genes <- counts.genes %>%
  mutate(number.wt = cv.wt > lower.wt & fit.wt < upper.wt)

counts.genes <- counts.genes %>%
  mutate(number.ago = cv.ago > lower.ago & fit.ago < upper.ago)

sum(counts.genes$number.wt == TRUE, na.rm = TRUE)
sum(counts.genes$number.wt == FALSE, na.rm = TRUE)

sum(counts.genes$number.ago == TRUE, na.rm = TRUE)
sum(counts.genes$number.ago == FALSE, na.rm = TRUE)


control.noiseplot <- ggplot(counts.genes_control) +
  geom_point(
    aes(x = log10(mean.all), y = log10(cv.wt)), color = "#F8766D", alpha = 0.1) +
  geom_point(
    aes(x = log10(mean.all), y = log10(cv.ago)), color = "#00BA38", alpha = 0.1) +
  geom_smooth(
    aes(x = log10(mean.all), y = log10(cv.wt)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#F8766D",
    fill = "#F8766D",
    alpha = 0.2) +
  geom_smooth(
    aes(x = log10(mean.all), y = log10(cv.ago)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#00BA38",
    fill = "#00BA38",
    alpha = 0.2) +
  theme_minimal() +
  theme(text = element_text(size = 30)) +
  labs(
    x = "Log10(mean)",
    y = "Log10(Coefficient Variation)") +
  coord_cartesian(
    xlim = c( 0, 3.5),
    ylim = c(-1, 0.5))

ethanol.noiseplot <- ggplot(counts.genes_ethanol) +
  geom_point(
    aes(x = log10(mean.all), y = log10(cv.wt)), color = "#F8766D", alpha = 0.1) +
  geom_point(
    aes(x = log10(mean.all), y = log10(cv.ago)), color = "#00BA38", alpha = 0.1) +
  geom_smooth(
    aes(x = log10(mean.all), y = log10(cv.wt)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#F8766D",
    fill = "#F8766D",
    alpha = 0.2) +
  geom_smooth(
    aes(x = log10(mean.all), y = log10(cv.ago)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#00BA38",
    fill = "#00BA38",
    alpha = 0.2) +
  theme_minimal() +
  theme(text = element_text(size = 30)) +
  labs(
    x = "Log10(mean)",
    y = "Log10(Coefficient Variation)") +
  coord_cartesian(
    xlim = c( 0, 3.5),
    ylim = c(-1, 0.5))

ggsave("ethanol.noiseplot.svg", plot = ethanol.noiseplot, width = 15, height = 15, units = "cm", dpi = 300)
ggsave("control.noiseplot.svg", plot = control.noiseplot, width = 15, height = 15, units = "cm", dpi = 300)

#####################################
#MODEL 
counts.genes_control$treatment <- rep("control", nrow(counts.genes_control))
counts.genes_control$geneid <- rownames(counts.genes_control)
counts.genes_control <- dplyr::filter(counts.genes_control, geneid != "Tf2-type")

counts.genes_ethanol$treatment <- rep("ethanol", nrow(counts.genes_ethanol))
counts.genes_ethanol$geneid <- rownames(counts.genes_ethanol)
counts.genes_ethanol <- dplyr::filter(counts.genes_ethanol, geneid != "Tf2-type")

counts.genes_control.long <- counts.genes_control %>%
  dplyr::select(cv.wt, cv.ago, treatment, mean.all, geneid) %>%
  pivot_longer(
    cols = c(cv.wt, cv.ago),
    names_to = "strain",
    values_to = "cv",
    names_prefix = "cv.")


counts.genes_ethanol.long <- counts.genes_ethanol %>%
  dplyr::select(cv.wt, cv.ago, treatment, mean.all, geneid) %>%
  pivot_longer(
    cols = c(cv.wt, cv.ago),
    names_to = "strain",
    values_to = "cv",
    names_prefix = "cv.")

model_control <- lmer(log(cv) ~ strain + (1 | geneid), data = counts.genes_control.long)
summary(model_control)
anova(model_control)

par(mfrow = c(2, 2))
plot(model_control)

model_ethanol <- lmer(log(cv) ~ strain + (1 | geneid), data = counts.genes_ethanol.long)
summary(model_ethanol)
anova(model_ethanol)

counts.all.long <- rbind(counts.genes_control.long, counts.genes_ethanol.long)

model_all <- lmer(log(cv) ~ strain * treatment + (1 | geneid), data = counts.all.long)
summary(model_all)
anova(model_all)
par(mfrow = c(2, 2))
plot(model_all)

ggplot(
  counts.all.long,
  aes(
    x = log10(mean.all),
    y = log10(cv),
    color = strain,
    linetype = treatment
  )
) +
 # geom_point(alpha = 0.1) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    linewidth = 1.5
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 25)
  ) +
  labs(
    x = "Log10(mean)",
    y = "Log10(Coefficient Variation)",
    color = "Strain",
    linetype = "Treatment"
  ) +
  coord_cartesian(xlim = c(NA, 3.5))

## raw data
raw.df <- as.data.frame(counts(set))
raw.df <- na.omit(raw.df)
raw.genes <- raw.df %>% filter(rownames(raw.df) %in% genes)
raw.spikes <- raw.df %>% filter(rownames(raw.df) %in% spikes)

WT.raw <- raw.genes %>% select(L41, L42, L43, L44, L45, L46, L47, L48, L49, L50)
ago.raw <- raw.genes %>% select(L51, L52, L53, L54, L55, L56, L57, L58, L59, L60)

raw.genes$cv.wt <- apply(WT.raw, 1, function(x) sd(x) / mean(x))
raw.genes$cv.ago <- apply(ago.raw, 1, function(x) sd(x) / mean(x))
raw.genes$mean.all <- apply(raw.genes, 1, function(x) mean(x))

raw.spikes$cv.spike <- apply(raw.spikes, 1, function(x) sd(x) / mean(x))
raw.spikes$mean.all <- apply(raw.spikes, 1, function(x) mean(x))


ggplot() +
  geom_point(data =raw.genes, aes(x = log10(mean.all), y = log10(cv.wt)),
             color = "#00688B", alpha=0.1) +
  geom_point(data =raw.genes, aes(x = log10(mean.all), y = log10(cv.ago)),
             color = "#FFA500", alpha=0.1) +
  geom_point(data =raw.spikes, aes(x = log10(mean.all), y = log10(cv.spike)),
             color = "#747474", alpha=0.1) +
  geom_smooth(data =raw.genes, aes(x = log10(mean.all), y = log10(cv.wt)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#00688B",
    fill = "#00688B",
    alpha = 0.2)+
  geom_smooth(data =raw.genes, aes(x = log10(mean.all), y = log10(cv.ago)),
    method = "loess",
    se = TRUE,
    level = 0.95,
    color = "#FFA500",
    fill = "#FFA500",
    alpha = 0.2)+
  geom_smooth(data =raw.spikes, aes(x = log10(mean.all), y = log10(cv.spike)),
              method = "loess",
              se = TRUE,
              level = 0.95,
              color = "#747474",
              fill = "#747474",
              alpha = 0.2)


       
                    ###############
                    #transposons
                    ##################

load("YES_count_genes.RData")
load("ethanol_count_genes.RData")

trans.YES <- counts.genes_control %>% 
  filter(row.names(counts.genes_control) == "Tf2-type")

trans.YES.long <- trans.YES %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "count")


strain <- (rep(c("WT", "Ago"), each=10))
condition <- (rep(c("control"), each= 20))
trans.YES.long$strain <- strain
trans.YES.long$condition <- condition

trans.etoh <- counts.genes_ethanol %>% 
  filter(row.names(counts.genes_ethanol) == "Tf2-type")

trans.etoh.long <- trans.etoh %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "count")

trans.etoh.long$strain <- strain
condition_e <- (rep(c("ethanol"), each= 20))
trans.etoh.long$condition <- condition_e

all_trans <- as.data.frame(rbind(trans.YES.long, trans.etoh.long))

CV_summary <- all_trans %>%
  summarise(
    .by = c(strain, condition),
    CV = sd(count) / mean(count))

boxplot_trans <- ggplot(all_trans, aes(x= condition, y=count, fill = strain))+
  geom_boxplot()+
  scale_fill_manual(values=c("#F8766D","#00BA38")) +
  theme_minimal(base_size = 25) +
  labs(x = NULL, y = "Tf2 read counts")

pointplot_trans <- ggplot(CV_summary, aes(x= condition, y=CV, color = strain))+
  geom_point()+
  scale_color_manual(values=c("#F8766D","#00BA38")) +
  theme_minimal(base_size = 25) +
  labs(x = "condition", y = "Coeff Variation")

trans.plot <- ggdraw() +
  draw_plot(boxplot_trans, x = 0, y = 0.4, width = 1, height = 0.55) +
  draw_plot(pointplot_trans, x = 0, y = 0, width = 1, height = 0.37) 


