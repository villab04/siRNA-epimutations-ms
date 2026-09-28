library(growthcurver)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(broom)

#setwd("~/Documents/Mariana/")
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
