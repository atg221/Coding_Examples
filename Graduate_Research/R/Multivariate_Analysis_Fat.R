library(readr)
library(dplyr)
setwd('/Users/owc/Downloads')
auto_manual <- read_csv('MSE_Test_S_V.csv')
lm_v <- lm(V_Manual~V_WF, auto_manual)
sm_v <- summary(lm_v)
mse_v <- mean(sm_v$residuals^2)
rmse_v <- sqrt(mean(sm_v$residuals^2))
lm_s <- lm(S_Manual~S_WF, auto_manual)
sm_s <- summary(lm_s)
mse_s <- mean(sm_s$residuals^2)
rmse_s <- sqrt(mean(sm_s$residuals^2))



setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr9/Chr9 Congenics/B6 vs Chr9 Congenic Fat MRI')
fat<-read_csv('B6_Chr9_Autoseg_Manova_Subcutaneous_Visceral.csv')
sig_fat<-read_csv('B6_Chr9_Autoseg_Manova_sigsub.csv')
man_visc_sum<-summary(manova(cbind(Slice,Visceral)~Strain, data = fat))
man_visc<-manova(cbind(Slice,Visceral)~Strain, data = fat)
man_sub_sum<-summary(manova(cbind(Slice,Subcutaneous)~Strain, data = fat))
man_sub<-manova(cbind(Slice,Subcutaneous)~Strain, data = fat)

man_sigsub<-summary(manova(cbind(Slice,Sub)~Strain, data = sig_fat))

