setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/Chr17 Aorta QTL (BalbxSM)/HDMP aorta microarray correlation (manual)/')

strain_name <- read_csv('HDMP_Name_Strain.csv');dim(strain_name)
strain_lesion <- read_csv('HDMP_Strain_Lesionavg.csv'); dim(strain_lesion)
gene <- read_csv('mep1a aorta expression (all mice).csv'); dim(gene)


strain_name_gene <- merge(strain_name, gene, by.x = "Name", by.y = "Name"); dim(strain_name_gene)
strain_name_gene_lesion <- merge(strain_name_gene, strain_lesion, by.x = "Strain", by.y = 'Strain', all.x = TRUE);dim(strain_name_gene_lesion)

#From WGCNA Package (Bicor)
bicor <- bicorAndPvalue(strain_name_gene_lesion$Value, strain_name_gene_lesion$lesion_avg)


qplot(con_mac_lesion$Hpcal1,con_mac_lesion$lesion_avg, xlab= 'Snx6', ylab = "Average Lesions Area (um²)") 
+ geom_abline(slope = -72350, intercept = 924802, col = "blue")