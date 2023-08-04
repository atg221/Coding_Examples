setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Male_Trei/')
library(plyr)
aorta <- read.csv('BalbxSMJ_Male_F2_Aortic_Trey.csv', header = T)
allmouseid <- subset(aorta, Mouse.ID != "")
allmouseid$Slide <- NULL
allmouseid$Section <- NULL
allmouseid$Count <- NULL
allmouseid$X <- NULL

write.csv(allmouseid, 'Balc_X_SMJ_F2_Heart_4x_Average_Max.csv')

#import geno/pheno data
geno_pheno <- read.csv('BALB_SMJ_male_F2_geno_phenotype-11-2014.csv', header = T)

#find the difference in column sizes so I can cbind the new lesion data to the rest of the 
#geno/pheno data; add that number of NAs to the data so it has the same number of nrows as other data.frame
x<-data.frame(matrix(NA, ncol = 2, nrow = 19))
colnames(x) <- colnames(aorta)
allmouseid_altered <- rbind(x, aorta)
colnames(allmouseid_altered)[1] <- "SampleID"

y<-data.frame(matrix(NA, ncol = 21, nrow = 140))
colnames(x) <- colnames(geno_pheno)
allmouseid_altered_2 <- cbind(y, geno_pheno)
#cbind that shit
fusion_hah <- merge(aorta, geno_pheno, by.x="Sample", by.y="Sample", all.x=TRUE, all.y=TRUE)
yar<- rbind(aorta, geno_pheno)
by_our_powers_combined <- join(allmouseid_altered, geno_pheno)

#write that shit
write.csv(fusion_hah, 'BALB_SMJ_male_F2_geno_phenotype_aortic_combined.csv')
write.csv(by_our_powers_combined, 'BALBxSM-F2-female-genopheno_joined_lesion_added.csv')
