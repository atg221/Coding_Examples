setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Human CT/UNET_CT_Backup/Data/TestingData')
m_a <- read.csv('Human_ct_training_validation_fat_volume_manual_auto_subq.csv')
#linear regression
m_a_lm_s <- lm(Manual_S~Auto_S, m_a)
m_a_lm_v <- lm(Manual_V~Auto_V, m_a)
#Get RMSE
rmse_s<- sqrt(mean(m_a_lm_s$residuals^2))
rmse_v<- sqrt(mean(m_a_lm_v$residuals^2))
#Print Residuals
m_a_lm_s_residual <- as.data.frame(m_a_lm_s$residuals)
write_csv(m_a_lm_s_residual, "subq_residual.csv")
m_a_lm_v_residual <-as.data.frame(m_a_lm_v$residuals)
write_csv(m_a_lm_v_residual, "visceral_residual.csv")
#to get R-squared
summary(m_a_lm_s)$adj.r.squared
summary(m_a_lm_v)$adj.r.squared
#to get p-value
summary(m_a_lm_s)$coefficients[,4]
summary(m_a_lm_v)$coefficients[,4]


#plot two to show correlation, r-squared, and p
v<-ggplot(m_a, aes(Manual_V, Auto_V))
vp <- v + geom_point(shape=1, size = 3)
vp + geom_smooth(method = 'lm', formula = y ~ x, color = "black", size = .5, se = F) +
annotate("text", x=50000, y=10000, label = "R^2 == 0.79", parse = T) +
 annotate("text", x=50000, y=8000, label = "p = 7.929e-31 ")






