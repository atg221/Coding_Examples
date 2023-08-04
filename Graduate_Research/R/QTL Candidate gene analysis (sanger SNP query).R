setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr12 carotid/')
setwd('/Users/owc/Downloads/')
chr5.qtl <- read.csv('sanger_snps_cath1_22-25cM_51887089-54649249mb.csv', header = T)

#make rows characters so can remove ones with no data
chr5.qtl$C57BL_6J <- as.character(chr5.qtl$C57BL_6J)
chr5.qtl$BALB_cJ <- as.character(chr5.qtl$BALB_cJ)
chr5.qtl$C3H_HeJ <- as.character(chr5.qtl$C3H_HeJ)
chr5.qtl$SM_J <- as.character(chr5.qtl$SM_J)
chr5.qtl$FVB.NJ <- as.character(chr5.qtl$FVB.NJ)
chr5.qtl$X129S1.SvImJ <- as.character(chr5.qtl$X129S1.SvImJ)
chr5.qtl$CAST_EiJ <- as.character(chr5.qtl$CAST_EiJ)
chr5.qtl$NZB.BlNJ <- as.character(chr5.qtl$NZB.BlNJ)
chr5.qtl$RF_J <- as.character(chr5.qtl$RF_J)
chr5.qtl$NZO_HlLtJ <- as.character(chr5.qtl$NZO_HlLtJ)
chr5.qtl$DBA_2J <- as.character(chr5.qtl$DBA_2J)
chr5.qtl$AKR_J <- as.character(chr5.qtl$AKR_J)
#subset for similar high and low alleles
s129.c3h.sim <- subset(chr5.qtl, SM_J == C57BL_6J & C3H_HeJ == BALB_cJ & C57BL_6J != C3H_HeJ); dim(s129.c3h.sim)
s129.c3h.sim <- subset(chr5.qtl, C3H_HeJ == BALB_cJ & C57BL_6J != C3H_HeJ); dim(s129.c3h.sim)
s129.c3h.sim <- subset(chr5.qtl, AKR_J == C3H_HeJ &  BALB_cJ == "-" & BALB_cJ == DBA_2J & C3H_HeJ != "-" & BALB_cJ != AKR_J & BALB_cJ != C3H_HeJ ); dim(s129.c3h.sim)
s129.c3h.bleh <- subset(s129.c3h.sim, C57BL_6J != " " &  BALB_cJ != " " &  C3H_HeJ != " " & SM_J != " " ); dim(s129.c3h.sim)

#look for instances where high and low differ
high.low.diff <- subset(s129.c3h.sim, S129S1_SvImJ != B6); dim(high.low.diff)
missense <- subset(high.low.diff, Csq =="missense_variant"); dim(missense)

#Export for Viewing
write.csv(s129.c3h.sim, "sanger_cath1_22-25cM_haplotype_analysis_no_SM.csv")

