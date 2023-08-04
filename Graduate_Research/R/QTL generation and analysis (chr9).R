setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr9/Glucose/Whole Genome CC/')
library(qtl)
library(qtlcharts)
data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAkAAAAJCAYAAADgkQYQAAAAMElEQVR42mNgIAXY2Nj8x8cHC8AwMl9XVxe3QqwKcJmIVwFWhehW4LQSXQCnm3ABAHD6MDrmRgfrAAAAAElFTkSuQmCC
#Load in cross(104 mice with both aortic and carotid lesion data)

a16_22_23_cc <- read.cross(format="csv", file="16387874 x 22294616_bb x 23938286 trial.csv", genotypes=c("High", "H", "Low", "D", "C"), na.strings=c("-"), convertXdata=TRUE)

a1_22_23 <- read_csv("cc_transformations.csv")

## REVIEW CC DATA ##
#Run for A_Max5
# read in cross
BxH_F2_M <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BALBxB6(F)_C3HxB6(F)_C3HxB6(M)/BxH F2 Male Rowlan-2013.csv", genotypes=c("C", "H", "B", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BxH_F2_F <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BALBxB6(F)_C3HxB6(F)_C3HxB6(M)/BxH F2 female Su-2005.csv", genotypes=c("A", "H", "B", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BB_F2_F <- read.cross(format="csv", file="BALB_B6_F2_Female .csv", genotypes=c("C", "H", "B", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BxS_F2_F <- read.cross(format="csv", file="BALB_SMJ_F2_Female.csv", genotypes=c("S", "H", "B", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BxS_F2_M <- read.cross(format="csv", file="BALB_SMJ_F2_Male.csv", genotypes=c("S", "H", "B", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BB.F_BxH.F_BxH.M <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BALBxB6(F)_C3HxB6(F)_C3HxB6(M)/balbxb6(F)_c3hxb6(F)_c3hxb6(M)_whole_genome_cc_sig_sug_marker_modded.csv", genotypes=c("High", "H", "Low", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
BB.F_BxH.F_BxH.M <- jittermap(BB.F_BxH.F_BxH.M)
# running scanone (GLUCOSE)
BB_F2_F_GlucoseW <- scanone(cross=BB_F2_F, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), pheno.col=c(23), model="normal", method="em")
BxS_F2_F_GlucoseW <- scanone(cross=BxS_F2_F, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), pheno.col=c(16), model="normal", method="em")
BxS_F2_M_GlucoseW<- scanone(cross=BxS_F2_M, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), pheno.col=c(23), model="normal", method='em')
BB.F_BxH.F_BxH.M_amax5_nonpar <- scanone(cross=BB.F_BxH.F_BxH.M, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), pheno.col=c(2), model="np")
#mqm augmentation + scan
a16_22_23_r_aug <- mqmaugment(a16_22_23_r)
a16_22_23_r.a16_22_23_mqm <- mqmscan(a16_22_23_r_aug, pheno.col = 2)
#plot permutation data
plot(a16_22_23_r.a16_22_23.permutations)
#visualize mqm permutations calculations
mqmplot.permutations(a16_22_23_r.a16_22_23_perm)
#Visualize heatmap (Use qtlcharts)
iplotScanone(BB_F2_F_GlucoseW_nonpar)
iplotScanone(BxS_F2_F_GlucoseW)
iplotScanone(a22.a22)
iplotScanone(a23.a23)
plot(a16.a16)
#visualize all crosses + cc together
plot(BB.F_BxH.F_BxH.M_amax5_nonpar) 
plot(BxH_F2_M_amax5_nonpar,BxH_F2_F_amax5_nonpar, BB_F2_F_amax5_nonpar, col = c("violetred", "blue", "green"), ylim = c(0, 9.161), chr = c(1:19)) 
plot(BB.F_BxH.F_BxH.M_amax5_nonpar, BxH_F2_M_amax5_nonpar, BxH_F2_F_amax5_nonpar) + plot(BB_F2_F_amax5_nonpar,col = "green", add=T)
plot(BB.F_BxH.F_BxH.M_amax5_nonpar, BxH_F2_M_amax5_nonpar, BB_F2_F_amax5_nonpar) 
plot(BB.F_BxH.F_BxH.M_amax5_nonpar, BxH_F2_M_amax5_nonpar, BxH_F2_F_amax5_nonpar, chr=9, col = c("black","violetred", "blue")) + plot(BB_F2_F_amax5_nonpar, chr=9, col = "green", add=T)
plot(BB.F_BxH.F_BxH.M_amax5_nonpar, BxH_F2_M_amax5_nonpar, BxH_F2_F_amax5_nonpar, chr=2) + plot(BB_F2_F_amax5_nonpar, chr=2, col = "green", add=T)
plot(BxH_F2_F_amax5_nonpar)
#find all suggestive + markers
a23.a23_df<- as.data.frame(a23.a23)
a23.a23_sug<- subset(a23.a23_df, lod >= 2)
a23_markers<- rownames(a23.a23_sug)
#find high/low allele at each of these markers
plotPXG(a23, "D15Mit171", pheno.col = 9)
## ***USE THIS TO SEE DISTRIBUTION OF LESION SIZE FOR EACH GENOTYPE AT A MARKER **** ##
plotPXG(a16_22_23_r, "D3Mit42", pheno.col = 3)


## BALBxSM F2 Carotid CC DATA ##
#read in
balbxs <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/BALB_SMJ_Female_F2_JQTL_Raw_Data.csv", genotypes=c("S", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
balbxb6 <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/Zhang2012_BL6xBALB_B37_Data.csv", genotypes=c("C", "H", "B", "D"), na.strings=c("NA"), convertXdata=TRUE)
c3hxb6 <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/Li2008_C3HxB6_B37_Data.csv", genotypes=c("C", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
balbxs_balbxb6_c3hxb6 <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/bxs_balbxb6_b6xc3h_combined.csv", genotypes=c("High", "H", "Low", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
#scanone
balbxs_balbxb6_c3hxb6.balbxs_balbxb6_c3hxb6 <- scanone(cross=balbxs_balbxb6_c3hxb6, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(2), model="np", method = 'em')
balbxs.balbxs <- scanone(cross=balbxs, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(16), model="np", method = 'em')
balbxb6.balbxb6 <- scanone(cross=balbxb6, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(12), model="np", method = 'em')
c3hxb6.c3hxb6 <- scanone(cross=c3hxb6, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(12), model="np", method = 'em')
#plot
plot(balbxs_balbxb6_c3hxb6.balbxs_balbxb6_c3hxb6, chr=12) + plot(balbxs.balbxs, balbxb6.balbxb6, c3hxb6.c3hxb6, chr=12, col = "black", lty = c(6,2,3), add=T)
plot(balbxs_balbxb6_c3hxb6.balbxs_balbxb6_c3hxb6, chr=13) + plot(balbxs.balbxs, balbxb6.balbxb6, c3hxb6.c3hxb6, chr=13,col = "black", lty = c(6,2,3), add=T)


## Combined Cross Analysis (Whole Genome Aorta) ##
#read in
balbxs_aorta <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BALBxSMJ F2 Female Aorta.csv", genotypes=c("S", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
balbxb6_aorta <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BB-F2 genotype and phenotype for JQTL CH12 updated 06_20_2012.csv", genotypes=c("C", "H", "B", "D"), na.strings=c("NA"), convertXdata=TRUE)
c3hxb6_aorta <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/BxH F2 female Su-2005.csv", genotypes=c("A", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
cc_whole_genome_multigeno_aorta <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/balbxsm_balbxb6_c3hxb6_whole_geno_cc_combined_multigeno.csv", genotypes=c("S", "H", "Ba", "C", "B"), na.strings=c("-"), convertXdata=TRUE)
cc_whole_genome_high_low_modded_aorta <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/Review/Combined cross (whole genome aorta)/balbxsm_balbxb6_c3hxb6_whole_geno_cc_combined_high_low_modded.csv", genotypes=c("High", "H", "Low", "D", "C"), na.strings=c("-"), convertXdata=TRUE)
jittermap(cc_whole_genome_high_low_modded_aorta)
#scanone
balbxs_aorta.avgmax5_nonpar <- scanone(cross=balbxs_aorta, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(16), model="np", method = 'em')
balbxb6_aorta.avgmax5_nonpar <- scanone(cross=balbxb6_aorta, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(23), model="np", method = 'em')
c3hxb6_aorta.avgmax5_nonpar <- scanone(cross=c3hxb6_aorta, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(6), model="np", method = 'em')
cc_whole_genome_multigeno_aorta.avgmax5_nonpar <- scanone(cross=cc_whole_genome_multigeno_aorta, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(2), model="np", method = 'em')
cc_whole_genome_high_low_modded_aorta.avgmax5_nonpar <- scanone(cross=cc_whole_genome_high_low_modded_aorta, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(2), model="np", method = 'em')
#plot
plot(balbxs_aorta.avgmax5_nonpar)
plot(balbxb6_aorta.avgmax5_nonpar)
plot(c3hxb6_aorta.avgmax5_nonpar)
plot(cc_whole_genome_multigeno_aorta.avgmax5_nonpar)
plot(cc_whole_genome_high_low_modded_aorta.avgmax5_nonpar)
#find all suggestive + markers
balbxs_aorta.avgmax5_nonpar_df<- as.data.frame(balbxs_aorta.avgmax5_nonpar)
balbxs_aorta.avgmax5_nonpar_sug<- subset(balbxs_aorta.avgmax5_nonpar_df, lod >= 2)
balbxs_markers<- rownames(balbxs_aorta.avgmax5_nonpar_sug)

balbxb6_aorta.avgmax5_nonpar_df<- as.data.frame(balbxb6_aorta.avgmax5_nonpar)
balbxb6_aorta.avgmax5_nonpar_sug<- subset(balbxb6_aorta.avgmax5_nonpar_df, lod >= 2)
balbxb6_markers <- rownames(balbxb6_aorta.avgmax5_nonpar_sug)

c3hxb6_aorta.avgmax5_nonpar_df<- as.data.frame(c3hxb6_aorta.avgmax5_nonpar)
c3hxb6_aorta.avgmax5_nonpar_sug<- subset(c3hxb6_aorta.avgmax5_nonpar_df, lod >= 2)
c3hxb6_markers <- rownames(c3hxb6_aorta.avgmax5_nonpar_sug)
#find high/low allele at each of these markers
plotPXG(balbxs_aorta, "mCV22888090", pheno.col = 16)
plotPXG(balbxb6_aorta, "D2Mit263", pheno.col = 23)
plotPXG(c3hxb6_aorta, "D19M103", pheno.col = 6)
plotPXG(cc_whole_genome_multigeno_aorta, "D11M214", pheno.col = 2)
plotPXG(cc_whole_genome_high_low_modded_aorta, "D1M512", pheno.col = 2)


## Combined Cross Analysis (Whole Genome Carotid) ##
#read in
balbxs.carotid <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/BALB_SMJ_Female_F2_JQTL_Raw_Data.csv", genotypes=c("S", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
balbxb6.carotid <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/Zhang2012_BL6xBALB_B37_Data.csv", genotypes=c("C", "H", "B", "D"), na.strings=c("NA"), convertXdata=TRUE)
c3hxb6.carotid <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/Li2008_C3HxB6_B37_Data.csv", genotypes=c("C", "H", "B", "D"), na.strings=c("-"), convertXdata=TRUE)
whole.genome.cc.female.carotid <- read.cross(format="csv", file="/Users/owc/Documents/UVA Grad Labs/Weibin/BALB x SMJ/Female Carotid Paper/Whole Genome CC/bxs_balbxb6_b6xc3h_whole_gneome_cc_sig_sug_marker_modded.csv", genotypes=c("High", "H", "Low", "D"), na.strings=c("-"), convertXdata=TRUE)
#scanone
balbxs.carotid.amax5.nonpar <- scanone(cross=balbxs.carotid, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(16), model="np", method = 'em')
balbxb6.carotid.amax5.nonpar <- scanone(cross=balbxb6.carotid, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(12), model="np", method = 'em')
c3hxb6.carotid.amax5.nonpar <- scanone(cross=c3hxb6.carotid, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(12), model="np", method = 'em')
whole.genome.cc.female.carotid.amax5.nonpar <- scanone(cross=whole.genome.cc.female.carotid, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(2), model="np", method = 'em')
whole.genome.cc.female.carotid.amax5.nonpar.perm <- scanone(cross=whole.genome.cc.female.carotid, chr=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"), pheno.col=c(2), model="np",n.perm = 1000, method = 'em')

#find all suggestive + markers
balbxs.carotid.amax5.nonpar.df<- as.data.frame(balbxs.carotid.amax5.nonpar)
balbxs.carotid.amax5.nonpar.sug<- subset(balbxs.carotid.amax5.nonpar.df, lod >= 2)
balbxs_markers<- rownames(balbxs.carotid.amax5.nonpar.sug)

balbxb6.carotid.amax5.nonpar.df<- as.data.frame(balbxb6.carotid.amax5.nonpar)
balbxb6.carotid.amax5.nonpar.sug<- subset(balbxb6.carotid.amax5.nonpar.df, lod >= 2)
balbxb6_markers <- rownames(balbxb6.carotid.amax5.nonpar.sug)

c3hxb6.carotid.amax5.nonpar.df<- as.data.frame(c3hxb6.carotid.amax5.nonpar)
c3hxb6.carotid.amax5.nonpar.sug<- subset(c3hxb6.carotid.amax5.nonpar.df, lod >= 2)
c3hxb6_markers <- rownames(c3hxb6.carotid.amax5.nonpar.sug)
#find high/low allele at each of these markers
plotPXG(balbxs.carotid, "rs3726547", pheno.col = 16)
plotPXG(balbxb6.carotid, "D5Mit95", pheno.col = 12)
plotPXG(c3hxb6.carotid, "D5Mit95", pheno.col = 12)
plotPXG(whole.genome.cc.female.carotid, "D5Mit309", pheno.col = 2)

