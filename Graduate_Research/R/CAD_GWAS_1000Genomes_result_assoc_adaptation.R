setwd('/Users/owc/Downloads/cad.additive.Oct2015.pub/')
p<- cad.add.160614.website[1:10,]
p1 <- cad.add.160614.website[,c(1:3,11)]
p2<- p1[,c(2,1,3:4)]
colnames(p2) <- c("CHR", "SNP", "BP", "P")
p1$P <- replace(p1$P, p1$P == 0, 1)
write.table(p2,"cad.add.160614.altered.2.txt", quote = FALSE, col.names = T, row.names = FALSE, sep="\t")

w1<- subset(cad.add.160614.altered.2, CHR <= 12)
w2<- subset(cad.add.160614.altered.2, CHR >12)
w <- w1[,c(2,4)]
ww <- w2[,c(2,4)]
write.table(w,"cad.add.160614.altered.locuszoom.1-12.txt", quote = FALSE, col.names = T, row.names = FALSE, sep="\t")
write.table(ww,"cad.add.160614.altered.locuszoom.12-22.txt", quote = FALSE, col.names = T, row.names = FALSE, sep="\t")

z<- subset(cad.add.160614.altered.locuszoom, P.value <= 1e-5)
write.table(z, "cad.add.160614.altered.locuszoom.main.effect.snps.txt", quote = FALSE, col.names = T, row.names = FALSE, sep="\t")

y <- subset(cad.add.160614.altered.2, CHR == 19)
oasl <- subset(y, BP > 121458095 & BP < 121477045 & P <= 1e-4); dim(oasl)
cxcl5 <- subset(y, BP > 74700000 & BP < 74900000 & P <= 1e-3); dim(cxcl5)
fgf5_prdm8 <- subset(y, BP > 81000000 & BP < 81400000 & P <= 1e-4); dim(fgf5_prdm8)
fam114a1_tlr1_tlr6 <- subset(y, BP > 38500000 & BP < 39500000 & P <= 1e-4)
cenpw_rspo3 <- subset(y, BP > 126600000 & BP < 127200000 & P <= 1e-3)
syt1_pawr <- subset(y, BP > 79500000 & BP < 80300000 & P <= 1e-3)
apoe <- subset(y, BP > 45300000 & BP < 45450000 & P<= 1e-4)

write.table(apoe,"apoe_cad_suggestive_snps.txt", quote = FALSE, col.names = T, row.names = FALSE, sep="\t")
