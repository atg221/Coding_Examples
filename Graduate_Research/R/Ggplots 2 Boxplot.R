library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Mep1a/Mep1a Immunostaining/Mep1aKO/')
#IF Cell Counting + Lesion Area
Mep_B6_IF_Cell <- read_csv("Mep1aKO_B6_IF_Aorta_Cell_Counting_SMC_Mac3.csv")
Mep_B6_IF_Cell_ly6g <- read_csv("Mep1aKO B6 Ly6g Lesion Area Comparison.csv")
Mep_B6_IF_Cell_cxcl5 <- read_csv("Mep1aKO B6 CXCL5 Lesion Area Comparison.csv")
#HE Plaque Stability
Mep_B6_HE_necrotic <- read_csv("Mep1aKO B6 HE Nectrotic.csv")
Mep_B6_tri_nonstained <- read_csv("Mep1aKO B6 Trichrome Non-Stained.csv")
Mep_B6_HE_fibrous_mid_out <- read_csv("Mep1aKO B6 HE Fibrous Cap Thickness Middle Outer.csv")
Mep_B6_HE_fibrous_min_max <- read_csv("Mep1aKO B6 Max Min Fibrous Cap Thickness.csv")

#plots
#HE
ggplot(Mep_B6_HE_fibrous_min_max, aes(Genotype, y = `Max Fibrous Cap Thickness`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20) +ylab("Max Fibrous Cap Thickness(um)")
#IF Cell Counting
ggplot(Mep_B6_IF_Cell, aes(Strain, y = `Percent Cells Mac3+`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20)+ aes(ymax = 100)
ggplot(Mep_B6_IF_Cell, aes(Strain, y = `Percent Cells SMC-A+`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20)+ aes(ymax = 100)
ggplot(Mep_B6_IF_Cell, aes(Strain, y = `Media Percent Cells SMC-A+`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20)+ aes(ymax = 100)

#IF Percent Area
ggplot(Mep_B6_IF_Cell_ly6g, aes(Strain, y = `Percent Ly6g+`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20)+ aes(ymax = 25)
ggplot(Mep_B6_IF_Cell_cxcl5, aes(Strain, y = `Percent CXCL5+`)) + geom_boxplot(fill = "grey") + geom_jitter( width = 0, alpha = .3) + theme_stata(base_size = 20)+ aes(ymax = 100)

