library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(dbplyr)
library(BiocFileCache)
library(rtracklayer)
library(AnnotationHub)

## Question 1: Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome (hg19)
## retrieve the cpg islands
ah <- AnnotationHub()
ah_human <- subset(ah, species == "Homo sapiens")
ah_human_cpg <- query(ah_human, "CpG Islands")
ah_human_cpg_record <- ah_human_cpg[[1]]
## extract autosomes
autosome_filter <- c(paste("chr", 1:22, sep=""))
split_record <- split(ah_human_cpg_record, seqnames(ah_human_cpg_record))
autosomes <- split_record[autosome_filter]
## check the number of autosomes
unlist(autosomes)


##Question 2: How many CpG Islands exists on chromosome 4.
autosomes[4]


## Question 3: Question 3
## Obtain the data for the H3K4me3 histone modification for the H1 cell line from Epigenomics Roadmap, 
## using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
## Question: How many bases does these regions cover?
#Query to find the right file (only want narrow peaks data)
ah_human_h3k4me3 <- query(ah_human, c("H3K4me3", "H1", "narrowpeak"))
#download the file
ah_human_h3k4me3_record <- ah_human_h3k4me3[[1]]
## extract autosomes and check the number of regions cover
split_H3K4me3 <- split(ah_human_h3k4me3_record, seqnames(ah_human_h3k4me3_record))
H3K4me3_autosomes <- split_H3K4me3[autosome_filter]
sum(width(unlist(H3K4me3_autosomes)))


##Question 4: Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, 
##using the AnnotationHub package. Subset these regions to only keep regions mapped to the autosomes. 
##In the return data, each region has an associated "signalValue".
##Question: What is the mean signalValue across all regions on the standard chromosomes?
#Query to find the right file (only want narrow peaks data)
ah_human_h3k27me3 <- query(ah_human, c("H3K27me3", "H1", "narrowpeak"))
#download the file
ah_human_h3k27me3_record <- ah_human_h3k27me3[[1]]
## extract autosomes
split_H3K27me3 <- split(ah_human_h3k27me3_record, seqnames(ah_human_h3k27me3_record))
H3K27me3_autosomes <- split_H3K27me3[autosome_filter]
## create a subset of extracted autosomes
ah_H3K27me3_autosomes <- subset(ah_human_h3k27me3_record, seqnames %in% autosome_filter)
#check mean signalValue
mean_signal_H3K27me3 <- mean(ah_H3K27me3_autosomes$signalValue)
mean_signal_H3K27me3


## Question 5: Bivalent regions are bound by both H3K4me3 and H3K27me3.
## Question: Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?
## intersect between two records
bivalent <- intersect(unlist(H3K4me3_autosomes), unlist(H3K27me3_autosomes))
sum(width(bivalent))


## Question 6: We will examine the extent to which bivalent regions overlap CpG Islands.
## Question: how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands?
# find bivalent regions overlap CpG Islands
cpg_autosomes <- autosomes
cpg_bivalent <- findOverlaps(bivalent, unlist(cpg_autosomes))
# calculate the fraction of the bivalent regions overlap CpG Islands
fraction_bivalent <- length(unique(queryHits(cpg_bivalent)))/length(bivalent)
#queryHits is used as that is the number of hits in the bivalent that overlap with cpg_autosomes
fraction_bivalent


## Question 7:Question: How big a fraction (expressed as a number between 0 and 1) 
## of the bases which are part of CpG Islands, are also bivalent marked.
cpg_bivalent_intersect <- intersect(bivalent, unlist(cpg_autosomes))
# calculate the fration of the bases intersected between CpG Islands and bivalent
fraction_bivalent_intersect <- sum(width(reduce(cpg_bivalent_intersect)))/sum(width(unlist(cpg_autosomes)))
fraction_bivalent_intersect


## Question 8: How many bases are bivalently marked within 10kb of CpG Islands?
## Tip: consider using the “resize()”" function.
# extract CpG Islands within 10kb
cpg_10kb <- resize(unlist(cpg_autosomes), width = 20000 + width(unlist(cpg_autosomes)), fix = "center")
cpg_10kb_bivalent <- intersect(cpg_10kb, bivalent)
sum(width(cpg_10kb_bivalent))


## Question 9: How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?
## Tip 1: the object returned by AnnotationHub contains "seqlengths".
## Tip 2: you may encounter an integer overflow. As described in the session on R Basic Types, 
## you can address this by converting integers to numeric before summing them, "as.numeric()".
# calculate human genome size
chr_list <- c(paste("chr", 1:22, sep=""))
genome <- keepSeqlevels(ah_human_cpg_record, chr_list, pruning.mode = "coarse")
genome_size <- sum(as.numeric(seqlengths(genome)))
# calculate the fraction of human genome which contained a CpG Island
cpg_autosomes_size <- sum(as.numeric(width(unlist(cpg_autosomes))))
cpg_autosomes_size / genome_size

## Question 10: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
## calculate InOut matrix
overlapMat <- matrix(0, ncol = 2, nrow = 2)
colnames(overlapMat) <- c("cpg_in", "cpg_out")
rownames(overlapMat) <- c("bivalent_in", "bivalent_out")
overlapMat[1,1] <- sum(width(cpg_bivalent_intersect))
overlapMat[1,2] <- sum(width(setdiff(bivalent, unlist(cpg_autosomes))))
overlapMat[2,1] <- sum(width(setdiff(unlist(cpg_autosomes), bivalent)))
overlapMat[2,2] <- genome_size - sum(overlapMat)
## calculate odds-ratio
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio









