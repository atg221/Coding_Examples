setwd('/Users/owc/Documents/UVA Grad Labs/Weibin/HOMER/C3H vs B6 Aortic Wall/')

rcn2<- read_csv("rcn2kd_vs_control_endothelial_cells_1.5_fold.csv")
rcn2_affy_d<- read.csv("b6_vs_c3h_western_aortic_wall_1.5_fold_up_affy.csv", header = F)
rcn2_affy_u<- read.csv("sw_vs_c3h_bone_marrow_1.5_fold_up_affy.csv", header = F)
rcn3_affy<- read.csv("rcn3kd_vs_control_endothelial_cells_1.5_fold_affy.csv", header = F)
sdf4_affy<- read.csv("sdf4kd_vs_control_endothelial_cells_1.5_fold_affy.csv", header = F)
chr12_affy<- read.csv("chr12_vs_b6_carotid_1.5_fold_affy_down.csv", header = F)
rcn3<- read_csv("rcn3kd_vs_control_endothelial_cells_1.5_fold.csv")
sdf4<- read_csv("sdf4kd_vs_control_endothelial_cells_1.5_fold.csv")



write_delim(rcn2_affy_d, 'b6_vs_c3h_western_aortic_wall_1.5_fold_up_affy.txt', delim = "/n")
write_delim(rcn2_affy_u, 'sw_vs_c3h_bone_marrow_1.5_fold_up_affy.txt', delim = "/n")
write_delim(sdf4, 'sdf4kd_vs_control_endothelial_cells_1.5_fold.txt', delim = "/n")

rcn2_rcn3 <- merge(rcn2, rcn3, by.x = "ID", by.y= "ID"); dim(rcn2_rcn3)
rcn2_rcn3_affy <- merge(rcn2_affy, rcn3_affy); dim(rcn2_rcn3_affy)
rcn2_rcn3_sdf4_affy <- merge(rcn2_rcn3_affy, sdf4_affy); dim(rcn2_rcn3_sdf4_affy)
rcn2_rcn3_sdf4 <-merge(rcn2_rcn3, sdf4, by.x = "ID", by.y="ID"); dim(rcn2_rcn3_sdf4)
write_delim(rcn2_rcn3_sdf4_affy, "rcn2_rcn3_sdf4_kd_shared_endothelial_cell_1.5_fold_affy.txt", delim = "/n")

shared<-read.csv("rcn2_rcn3_sdf4_kd_shared_endothelial_cell_1.5_fold.csv", header=F)
write_delim(shared, 'rcn2_rcn3_sdf4_kd_shared_endothelial_cell_1.5_fold.txt', delim = "/n")
