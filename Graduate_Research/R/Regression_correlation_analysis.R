setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female/Aorta/New Samples/')
setwd('/Users/owc/Downloads/')
pheno <- read.csv('pheno_non_q_redo_checked.csv', header = T)
pheno2<-read.csv('all_data_only_mice_with_carotid.csv', header = T)
pheno2<- pheno2[3:nrow(pheno2),]
genotypes<-pheno2[,c(1,32:ncol(pheno2))]
p <- read.csv('hdl_vs_aorta_b6_balb.csv')
#Dem Tests Tho
#basic correlation test
cor.test(pheno2$Log2_Carotid_Lesion_Total_Area,pheno2$NonHDL_W)
cor.test(p$HDL.Western, p$Aortic.lesion_Aveg)
#linear regression + Anova
regression <- lm(p$HDL.Western~p$Aortic.lesion_Aveg)
all_pheno<-lm(pheno$Aorta_Sum_Top_8~pheno$Totalchol_W+pheno$HDL_W+pheno$NonHDL_W+pheno$Triglyceride_W+pheno$Glucose_W+pheno$Log2_Carotid_Lesion_Total_Area)
re<-lm(pheno$Aorta_Sum_Top_8~pheno$HDL_W)
sum <- summary(regression)
#isolate the R-squared and p values from the regression summary
sum_r<-sum$adj.r.squared
sum_r<-format(sum_r, digits = 3)
sum_p<-sum$coefficients[2,4]
sum_p<-format(sum_p, digits=3)
#check with anova
anova_pheno<-anova(all_pheno)
#plot scatterplot and add sample names and R/p values
regress<- qplot(p$Aortic.lesion_Aveg,p$HDL.Western, xlab = "Aortic Root Lesion Area (um²)", ylab = "HDL (mg/dl)") + geom_abline(slope = -1.810e-04, intercept = 9.226e+01)
text(y=12.5, x=37, labels =paste("R²=",sum_r), cex=.75);text(y=11.5, x=37, labels=paste("p=",sum_p), cex=.75)
text(pheno$Carotid.Lesion.Total.Area, pheno$Aorta_Lesion_Avg_Top_5, labels = (pheno2$SampleID), cex=0.6, pos=4, col="red")
#Hits (With 'BALBxSMJ_F2_Heart_Female_Lesion_Geno_All_Mice.csv'):
  #Negative corelation between HDL_W and Average..Top.5
  #

###
#do cluster analysis
###
#make a matrix of the data you want to cluster
matrix<-data.matrix(pheno2$Log2_Carotid_Lesion_Total_Area)
matrix2<-data.matrix(pheno2$Carotid_Lesion_Total_Area)
#Figure out number of clusters
cnum <- (nrow(matrix2)-1)*sum(apply(matrix2,2,var))
for (i in 2:15){
  cnum[i] <- sum(kmeans(matrix2, centers=i)$withinss)}
plot(1:15, cnum, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
#use kmeans to cluster, specifying how many clusters you want(in my case, I can see there is probably 2 for Log10_Carotid;3 for Carotid Total Area)
c<-kmeans(matrix,2)
c2<-kmeans(matrix2,3)
#look at clusters to see what their means are
aggregate(matrix,by=list(c$cluster),FUN=mean)
aggregate(matrix2,by=list(c2$cluster),FUN=mean)
#add cluster number to original data as a data.frame
cluster<-data.frame(matrix,c$cluster)
cluster2<-data.frame(matrix2,c2$cluster)
rownames(cluster)=pheno2$SampleID
rownames(cluster2)=pheno2$SampleID
#make data.frame with SampleIDs/caroti data/cluster
cluster_sampleid<-data.frame(pheno2$SampleID,matrix,c$cluster)
cluster_sampleid2<-data.frame(pheno2$SampleID,matrix2,c2$cluster)
colnames(cluster_sampleid2)<-c("SampleID", "Carotid_Lesion_Total_Area","Cluster")
colnames(cluster_sampleid)<-c("SampleID", "Log2_Carotid_Lesion_Total_Area","Cluster")
# Ward Hierarchical Clustering
d <- dist(cluster, method = "euclidean") # distance matrix
d2 <- dist(cluster2, method = "euclidean")
fit <- hclust(d, method="ward") 
fit2 <- hclust(d2, method="ward") 
plot(fit, xlab="SampleID", main = "Carotid Lesion (Total Area)") # display dendogram
plot(fit2, xlab="SampleID", main = "Carotid Lesion (Total Area)")
groups <- cutree(fit, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=2, border="red")
rect.hclust(fit2, k=3, border="red")

#see what alleles have bigger average lesion
cluster_genotypes_2<-merge(cluster_genotypes, pheno[,c(1,23)], by.x = "SampleID", by.y="SampleID")
cluster_3_13_AA<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 3 & cluster_genotypes$gnf13.092.499 == "AA")
cluster_3_13_AA<-cluster_3_13_AA[,c(1:3,24,25,26,148)]
cluster_3_13_AB<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 3 & cluster_genotypes$gnf13.092.499 == "AB")
cluster_3_13_AB<-cluster_3_13_AB[,c(1:3,24,25,26,148)]
cluster_3_13_BB<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 3 & cluster_genotypes$gnf13.092.499 == "BB")
cluster_3_13_BB<-cluster_3_13_BB[,c(1:3,24,25,26,148)]
cluster_2_13_AA<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 2 & cluster_genotypes$gnf13.092.499 == "AA")
cluster_2_13_AA<-cluster_2_13_AA[,c(1:3,24,25,26,148)]
cluster_2_13_AB<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 2 & cluster_genotypes$gnf13.092.499 == "AB")
cluster_2_13_AB<-cluster_2_13_AB[,c(1:3,24,25,26,148)]
cluster_2_13_BB<-subset(cluster_genotypes_2, cluster_genotypes$Cluster == 2 & cluster_genotypes$gnf13.092.499 == "BB")
cluster_2_13_BB<-cluster_2_13_BB[,c(1:3,24,25,26,148)]
cluster_3_13<-rbind(cluster_3_13_AA,cluster_3_13_AB,cluster_3_13_BB)
cluster_2_13<-rbind(cluster_2_13_AA,cluster_2_13_AB,cluster_2_13_BB)
cluster_13<-rbind(cluster_3_13,cluster_2_13)
write.csv(cluster_3_13, "Cluster_3_Chr_13_120cM_Carotid_Genotype.csv")
write.csv(cluster_2_13, "Cluster_2_Chr_13_120cM_Carotid_Genotype.csv")
write.csv(cluster_13, "Chr2_97cM_111cM_127cM_Carotid_Aortic_Genotype.csv")
write.csv(cluster_genotypes_2, "Cluster_Carotid_Aortic_Genotype.csv")
###
#look at genotypes of samples in each cluster at Chr13 peak
###
cluster_phenotypes<-merge(cluster_sampleid, pheno2, by.x="SampleID", by.y="SampleID")
cluster_pheno_geno <- merge(cluster_phenotypes, genotypes, by.x="SampleID", by.y="SampleID")
write.csv(cluster_pheno_geno, "Log2_Carotid_Lesion_Total_Area_Clustered_With_All_Phenotypes_And_Geno.csv")

####
#Sort Log2 Carotid clustered data by cluster, then carotid size, then aortic side, then other
#phenos to find healthiest/sickest
###
cluster_pheno_sorted <- sort(cluster_genotypes$Cluster)
                               , cluster_genotypes$Log2_Carotid_Lesion_Total_Area.x, cluster_genotypes$Aorta_Average_Top_5, cluster_genotypes$Glucose_W))

