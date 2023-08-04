library(ggplot2)
library(ggthemes)
setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr9/Chr9 Congenics/B6 vs Chr9 Congenic Fat MRI')
visc <- read_csv('Manual_WF_Normal_Visceral_Comparison.csv')
v<-ggplot(visc, aes(Manual, Normal))
vp <- v + geom_point(shape=1, size = 3)
vp + geom_abline(intercept = 0, slope = 1, color = "black", alpha = .7) + ylim(0,600)+ xlim(0,600) + ylab("Normal") +theme_bw(base_size = 20)


sub <- read_csv('Manual_WF_Normal_Subcutaneous_Comparison.csv')
s<-ggplot(sub, aes(Manual, Water_Filtered))
sp <- s + geom_point(shape=1, size = 3)
sp + geom_abline(intercept = 0, slope = 1, color = "black", alpha = .7) + ylim(0,600)+ xlim(0,600) + ylab("Water-Filtered") +theme_bw(base_size = 20)
