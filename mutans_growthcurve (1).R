library(growthcurver)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(broom)

#setwd("~/Documents/Mariana/")
setwd("/home/bioc1877-local/Documents/MA_pombe/growth_curve/")
data <- read.csv("mutants_growthcurve.csv")
head(data)

gc_mutantout <- SummarizeGrowthByPlate(data, plot_fit = TRUE, 
                                 plot_file = "gc_mutantplots.pdf")
head(gc_mutantout)

write.csv(gc_mutantout, "gc_mutantout.csv")
############################################################### 
data_gc <- read.csv("gc_mutantout_edit.csv")
head(data_gc)
data_gc$repetition <- factor(data_gc$repetition)
data_gc$treatment <- factor(data_gc$treatment)
data_gc$treatment <- relevel(data_gc$treatment, ref = "YES")

ggplot(data_gc, aes(x = treatment, y = r, fill = strain)) +
  geom_boxplot() +
  facet_wrap(~ strain) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(data_gc, aes(x = treatment, y = r, fill = strain)) +
  geom_boxplot() +
  #facet_wrap(~ strain) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(data_gc, aes(x = treatment, y = k, fill = strain)) +
  geom_boxplot() +
  facet_wrap(~ strain)

#statistical model
model1 <-lmer(r ~ treatment + strain + (1 | repetition), data = data_gc)
summary(model1)
anova(model1)
 
rdp1_data <- data_gc %>% filter(strain == "Δ rdp1")
chp2_data <- data_gc %>% filter(strain == "Δ chp2")
ago1_data <- data_gc %>% filter(strain == "Δ ago1")

model_rpd1 <-glm(r ~ treatment, data = rdp1_data)
summary(model_rpd1)

model1_chp2 <-glm(r ~ treatment, data = chp2_data)
summary(model1_chp2)

model1_ago1 <-glm(r ~ treatment, data = ago1_data)
summary(model1_ago1)

##################################################################
#T1 05.02.2026. 
data_t1 <- read.csv("OD600_T1.csv")
head(data_t1)

gc_t1_out <- SummarizeGrowthByPlate(data_t1, plot_fit = TRUE, 
                                       plot_file = "gc_mutantplots.pdf")
head(gc_t1_out)

write.csv(gc_t1_out, "gc_t1_out.csv")
#T2 07.02.2026
data_t2 <- read.csv("OD600_T2.csv")
head(data_t2)

gc_t2_out <- SummarizeGrowthByPlate(data_t2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t2.pdf")

write.csv(gc_t2_out, "gc_t2_out.csv")
#T3 09.02.206
data_t3 <- read.csv("OD600_T3.csv")
head(data_t3)

gc_t3_out <- SummarizeGrowthByPlate(data_t3, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t3.pdf")

write.csv(gc_t3_out, "gc_t3_out.csv")

#T4 11.02.206
data_t4 <- read.csv("OD600_T4.csv")
head(data_t4)

gc_t4_out <- SummarizeGrowthByPlate(data_t4, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t4.pdf")

write.csv(gc_t4_out, "gc_t4_out.csv")
#T5 13.02.206
data_t5 <- read.csv("OD600_T5.csv")
head(data_t5)

gc_t5_out <- SummarizeGrowthByPlate(data_t5, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t5.pdf")

write.csv(gc_t5_out, "gc_t5_out.csv")

#Load growth curve parameters
data_edit_T1 <- read.csv("gc_t1_out_edit.csv")
data_edit_T2 <- read.csv("gc_t2_out_edit.csv")
data_edit_T2 <- data_edit_T2[-37,]
data_edit_T3 <- read.csv("gc_t3_out_edit.csv")
data_edit_T4 <- read.csv("gc_t4_out_edit.csv")
data_edit_T5 <- read.csv("gc_t5_out_edit.csv")

all_transfers <- rbind(data_edit_T1, data_edit_T2, data_edit_T3, data_edit_T4, data_edit_T5)
WT_alltransfers <- filter(all_transfers, strain == "WT")
rdp1_alltransfers <- filter(all_transfers, strain == "Δ rdp1")
ago1_alltransfers <- filter(all_transfers, strain == "Δ ago1")
chp2_alltransfers <- filter(all_transfers, strain == "Δ chp2")

models <- all_transfers %>%
  group_by(sample, treatment, strain) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(r ~ transfer, data = .x)),
    results = map(model, tidy)
  ) %>% unnest(results)

models2 <- all_transfers %>%
  group_by(treatment, strain) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(r ~ transfer, data = .x)),
    results = map(model, tidy)
  ) %>% unnest(results)

models1<- filter(models, term == "transfer")

model1 <- lm(estimate ~ treatment * strain, data = models)
summary(model1)

ggplot(all_transfers , aes(x = transfer, y = r, color = strain)) +
 #geom_point() +
  geom_smooth(method = "lm", se = TRUE)+
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())



ggplot(data_edit_T1, aes(x = treatment, y = r, fill = strain)) +
  geom_boxplot() +
  facet_wrap(~ strain) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(data_edit_T2, aes(x = treatment, y = r, fill = strain)) +
  geom_point() +
  facet_wrap(~ strain) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(data_edit_T3, aes(x = treatment, y = r, fill = strain)) +
  geom_boxplot() +
  facet_wrap(~ strain) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(WT_alltransfers, aes(x = transfer, y = r)) +
  geom_point() +
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")


ggplot(rdp1_alltransfers, aes(x = transfer, y = r)) +
  geom_point() +
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(ago1_alltransfers, aes(x = transfer, y = r)) +
  geom_point() +
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

ggplot(chp2_alltransfers, aes(x = transfer, y = r)) +
  geom_point() +
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none")

###############################################################################
# EXP2 DATA

##################################################################
#T1 22.02.2026. 
data_t1.2 <- read.csv("OD600_T1_EXP2.csv")
head(data_t1.2)

gc_t1.2_out <- SummarizeGrowthByPlate(data_t1.2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_EXP2.pdf")
head(gc_t1.2_out)

write.csv(gc_t1.2_out, "gc_t1.2_out.csv")
#T2 24.02.2026
data_t2.2 <- read.csv("OD600_T2_EXP2.csv")
head(data_t2.2)

gc_t2.2_out <- SummarizeGrowthByPlate(data_t2.2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t2.2.pdf")

write.csv(gc_t2.2_out, "gc_t2.2_out.csv")
#T3 26.02.206
data_t3.2 <- read.csv("OD600_T3_EXP2.csv")
head(data_t3.2)

gc_t3.2_out <- SummarizeGrowthByPlate(data_t3.2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t3.2.pdf")

write.csv(gc_t3.2_out, "gc_t3.2_out.csv")

#T4 28.02.206
data_t4.2 <- read.csv("OD600_T4_EXP2.csv")
head(data_t4.2)

gc_t4.2_out <- SummarizeGrowthByPlate(data_t4.2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t4.2.pdf")

write.csv(gc_t4.2_out, "gc_t4.2_out.csv")
#T5 02.03.206
data_t5.2 <- read.csv("OD600_T5_EXP2.csv")
head(data_t5.2)

gc_t5.2_out <- SummarizeGrowthByPlate(data_t5.2, plot_fit = TRUE, 
                                    plot_file = "gc_mutantplots_t5.2.pdf")

write.csv(gc_t5.2_out, "gc_t5.2_out.csv")

#Load growth curve parameters exp2
data_edit_T1.2 <- read.csv("gc_t1.2_out_edit.csv")
#data_edit_T1.2 <- data_edit_T2.2[-c(25,10),]
data_edit_T2.2 <- read.csv("gc_t2.2_out_edit.csv")
data_edit_T3.2 <- read.csv("gc_t3.2_out_edit.csv")
data_edit_T4.2 <- read.csv("gc_t4.2_out_edit.csv")
data_edit_T5.2 <- read.csv("gc_t5.2_out_edit.csv")

all_transfers_exp2 <- rbind(data_edit_T1.2, data_edit_T2.2, 
                            data_edit_T3.2, data_edit_T4.2, data_edit_T5.2)

models_exp2 <- all_transfers_exp2 %>%
  group_by(sample, treatment, strain) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(r ~ transfer, data = .x)),
    results = map(model, tidy)
  ) %>% unnest(results)

ggplot(all_transfers_exp2, aes(x = as.factor(transfer), y = r, color = strain)) +
  geom_boxplot() +
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())

ggplot(all_transfers_exp2 , aes(x = transfer, y = r, color = strain)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)+
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())


#Bind the two databases
all_transfers$experiment <- rep_len("exp1", nrow(all_transfers))
all_transfers <- all_transfers %>% select(!note) %>% filter(treatment != "YES – 6%")
all_transfers_exp2$experiment <- rep_len("exp2", nrow(all_transfers_exp2))

all_exp <- rbind(all_transfers_exp2, all_transfers)

CoeffVar <- all_exp %>% group_by(strain, treatment) %>% 
  mutate(mean_strain = mean(r)) %>% ungroup() %>% group_by(sample) %>% 
  mutate(sd_line= sd(r)) %>% ungroup() %>% mutate(cv = (sd_line/mean_strain))

as.data.frame(CoeffVar)


ggplot(all_exp , aes(x = transfer, y = r, color = strain)) +
  geom_point() +
  geom_smooth(aes(group = sample), method = "lm", se = FALSE)+
  facet_wrap(~ treatment + experiment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())

ggplot(all_exp , aes(x = transfer, y = r, color = strain)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)+
  facet_wrap(~ treatment + experiment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())

ggplot(all_exp , aes(x = transfer, y = r, color = strain)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)+
  facet_wrap(~ treatment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())

ggplot(all_exp, aes(x = as.factor(transfer), y = r, color = strain)) +
  geom_boxplot() +
  facet_wrap(~ treatment + experiment) +
  labs(y = "Growth rate (s⁻¹)") +
  theme_light() +
  theme(
    axis.title.x = element_blank())

ggplot(CoeffVar , aes(x = strain, y = cv, color = strain)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ treatment ) +
  labs(y = "Coefficient ofa variation") +
  theme_light() 

  theme(
    axis.title.x = element_blank())