setwd('/Users/normangarrett/Desktop')

#load in
BxS.M.F2 <- read.cross(format = "csv", file="BALB_SMJ_male_F2_geno_phenotype_aortic_combinedORGINAL.csv", genotypes=c("S","H","B","D","C"), na.strings=c("-"),convertXdata = T)
BxS.M.F2_lesion_par <- scanone(cross=BxS.M.F2, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), pheno.col=c(6), model="normal")

#simulate the genotypes between markers given the observed marker data
BxS.M.F2.sim <- sim.geno(BxS.M.F2)

#make a "fake" qtl at the nearest marker to the peak 
#can do separately or combined (separate seems to give more accurate LOD scores than grouping them all together)
BxS.M.F2.chr10 <- makeqtl(BxS.M.F2.sim, chr = 10, pos = 9.158)
BxS.M.F2.chr9 <- makeqtl(BxS.M.F2.sim, chr = 9, pos = 26.842)
BxS.M.F2.qtl <-  makeqtl(BxS.M.F2.sim, chr = c(9,10), pos = c(26.842, 9.158))

#fit the QTL back to the simulated genotype cross to calculate the LOD score and effect size of the QTL (%var)
effect.size <- fitqtl(BxS.M.F2.sim, pheno.col = 6, qtl = BxS.M.F2.qtl)
effect.size.10 <- fitqtl(BxS.M.F2.sim, pheno.col = 6, qtl = BxS.M.F2.chr10)
effect.size.9 <- fitqtl(BxS.M.F2.sim, pheno.col = 6, qtl = BxS.M.F2.chr9)
