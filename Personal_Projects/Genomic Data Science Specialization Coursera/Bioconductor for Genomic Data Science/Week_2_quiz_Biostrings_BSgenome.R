#Load in relevant libraries
library(AnnotationHub)
library(Biostrings)
library(BSgenome)
  library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
#install genomes if needed
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
#load genomes (can use these with Views later on to look at regions that overlap ranges from relevant tracks)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#Question 1:
# Question: What is the GC content of “chr22” in the “hg19” build of the human genome?
# Tip: The reference genome includes “N” bases; you will need to exclude those.
# count total bases on chr22
alphabet_frequency <- alphabetFrequency(Hsapiens$chr22)
total_bases <- sum(alphabet_frequency[c('A','G','C','T')])
# count GC bases on chr22  
GC_bases <- sum(alphabet_frequency[c('G','C')])
#calculate GC ratio
GC_content <- GC_bases/total_bases
GC_content


# Question 2: 
# Background: In the previous assessment we studied H3K27me3 “narrowPeak” regions from the H1 cell line 
# (recall that the Roadmap ID for this cell line is “E003”). We want to examine whether the GC content of the regions influence the signal; 
# in other words wether the reported results appear biased by GC content.
# Question: What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics Roadmap from the H1 stem cell line on chr 22.
# Clarification: Compute the GC content for each peak region as a percentage and then average those percentages to compute a number between 0 and 1.
# retrieve record
ah <- AnnotationHub()
H3K27me3_qh <- query(ah, c("H3K27me3", "E003", "narrowPeak"))
H3K27me3_record <- H3K27me3_qh[["AH29892"]]
# extract chr 22
H3K27me3_chr22 <- subset(H3K27me3_record, seqnames == "chr22")
#use Views to grab the subset that only contains peak sequences
H3K27me3_chr22_views <- Views(Hsapiens, H3K27me3_chr22)
# calculate mean GC content
GC_contents <- letterFrequency(H3K27me3_chr22_views, "GC", as.prob = TRUE)
mean_GC <- mean(GC_contents)
mean_GC

# Question 3:
# The “narrowPeak” regions includes information on a value they call “signalValue”.
# Question: What is the correlation between GC content and “signalValue” of these regions (on chr22)?
signal_value <- mcols(H3K27me3_chr22_views)$signalValue
#use "cor" to find the correlation
cor(signal_value, GC_contents)



# Question 4
# The “narrowPeak” regions are presumably reflective of a ChIP signal in these regions. To confirm this, we want to obtain the “fc.signal” data 
# from AnnotationHub package on the same cell line and histone modification. This data represents a vector of fold-change enrichment of ChIP signal over input.
# Question: what is the correlation between the “signalValue” of the “narrowPeak” regions and the average “fc.signal” across the same regions?
# Clarification: First compute the average “fc.signal” for across each region, for example using “Views”; this yields a single number of each region. 
# Next correlate these numbers with the “signalValue” of the “narrowPeaks”.
# retrieve record
H3K27me3_fc <- query(ah, c("H3K27me3", "E003", "fc.signal"))
H3K27me3_fc_record <- H3K27me3_fc[["AH32033"]]
# get subset data on chr22
gr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
#import the file with the specific GRange you want as an Rle
H3K27me3_fc_gr <- import(H3K27me3_fc_record, which = gr22, as = "Rle")
H3K27me3_fc_gr22 <- H3K27me3_fc_gr$chr22
# view fc.signal data
fc.signal <- Views(H3K27me3_fc_gr22, start = start(H3K27me3_chr22), end = end(H3K27me3_chr22))
# calculate the correlation between the average of fc.signal and signalValue
fc.signal_mean <- mean(fc.signal)
cor(fc.signal_mean, signal_value)

# Question 5
# Referring to the objects made and defined in the previous question.
# Question: How many bases on chr22 have an fc.signal greater than or equal to 1?
sum(H3K27me3_fc_gr22 >= 1)


# Question 6
# The H1 stem cell line is an embryonic stem cell line, a so-called pluripotent cell. Many epigenetic marks change upon differentiation. 
# We will examine this. We choose the cell type with Roadmap ID “E055” which is foreskin fibroblast primary cells.
# We will use the “fc.signal” for this cell type for the H3K27me3 mark, on chr22. We now have a signal track for E003 and a signal track for 
# E055. We want to identify regions of the genome which gain H3K27me3 upon differentiation. These are regions which have a higher signal 
# in E055 than in E003. To do this properly, we would need to standardize (normalize) the signal across the two samples; we will ignore this for now.
# Question: Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.
# Tip: If you end up with having to intersect two different Views, note that you will need to convert the Views to IRanges or GRanges first 
# with ir <- as(vi, "IRanges")ir <- as(vi, "IRanges").
# get record
H3K27me3_E055 <- query(ah, c("H3K27me3", "E055"))
H3K27me3_E055_record <- H3K27me3_E055[["AH32470"]]
# get subset data on chr22
gr_chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
H3K27me3_fc_gr_E055 <- import(H3K27me3_E055_record, which = gr_chr22, as = "Rle")
H3K27me3_fc_gr22_E055 <- H3K27me3_fc_gr_E055$chr22
# identify region (all regions in which E003 has a signal < 0.5 and E055 has a signal >2)
region_E003 <- as(slice(H3K27me3_fc_gr22, upper = 0.5), "IRanges")
region_E055 <- as(slice(H3K27me3_fc_gr22_E055, lower = 2), "IRanges")
inter_region <- intersect(region_E003, region_E055)
sum(width(inter_region))


# Question 7
# CpG Islands are dense clusters of CpGs. The classic definition of a CpG Island compares the observed to the expected frequencies of CpG 
# dinucleotides as well as the GC content.
# Specifically, the observed CpG frequency is just the number of “CG” dinucleotides in a region. The expected CpG frequency is defined as 
# the frequency of C multiplied by the frequency of G divided by the length of the region.
# Question: What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?
## retrieve the cpg
ah_human <- subset(ah, species == "Homo sapiens")
ah_human_cpg <- query(ah_human, "CpG Islands")
ah_human_cpg_record <- ah_human_cpg[["AH5086"]]
# get subset data on chr22
ah_human_cpg_chr22 <- subset(ah_human_cpg_record, seqnames == "chr22")
ah_human_cpg_chr22_views <- Views(Hsapiens, ah_human_cpg_chr22)
# calculate observed GC bases
observed_GC <- dinucleotideFrequency(ah_human_cpg_chr22_views)[,7]/width(ah_human_cpg_chr22_views)
# calculate expected GC bases
freq_C <- letterFrequency(ah_human_cpg_chr22_views, "C")
freq_G <- letterFrequency(ah_human_cpg_chr22_views, "G")
expected_GC <- (freq_C/width(ah_human_cpg_chr22_views))*(freq_G/width(ah_human_cpg_chr22_views))
# calculate the average observed-to-expected ratio of CpG dinucleotides
mean(observed_GC/expected_GC)


# Question 8:
# A TATA box is a DNA element of the form “TATAAA”. Around 25% of genes should have a TATA box in their promoter. 
# We will examine this statement.
# Question: How many TATA boxes are there on chr 22 of build hg19 of the human genome?
# Clarification: You need to remember to search both forward and reverse strands.
TATA_boxes <- countPattern("TATAAA", Hsapiens$chr22) + countPattern("TATAAA", reverseComplement(Hsapiens$chr22))
TATA_boxes

# Question 9:
# Question: How many promoters of transcripts on chromosome 22 containing a coding sequence, contains a TATA box on the same strand as 
# the transcript?
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. 
# Here, we defined a promoter to be 900bp upstream and 100bp downstream of the transcription start site.
#load up the transcript database and set the GRanges
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
# find promoters of transcripts on chr 22
gr_trans_chr22 <- subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
proms <- promoters(gr_trans_chr22, upstream = 900, downstream = 100)
# find coding sequences on chr 22
gr_cds_chr22 <- subsetByOverlaps(cds(txdb), gr, ignore.strand = TRUE)
# find overlaps between promoters of transcripts and coding sequences
gr_proms_cds_chr22 <- subsetByOverlaps(proms, gr_cds_chr22, ignore.strand = TRUE)
# calculate TATA box on overlaps
gr_proms_cds_chr22_views <- Views(Hsapiens, gr_proms_cds_chr22)
TATA_boxes_prom_cds_chr22 <- vcountPattern("TATAAA", DNAStringSet(gr_proms_cds_chr22_views))
sum(TATA_boxes_prom_cds_chr22)

# Question 10:
# It is possible for two promoters from different transcripts to overlap, in which case the regulatory features inside the overlap might 
# affect both transcripts. This happens frequently in bacteria.
# Question: How many bases on chr22 are part of more than one promoter of a coding sequence?
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. Here, we define a promoter to 
# be 900bp upstream and 100bp downstream of the transcription start site. In this case, ignore strand in the analysis.
# calculate transcript lengths
trans_len_chr22 <- transcriptLengths(txdb, with.cds_len = TRUE)
trans_len_chr22 <- trans_len_chr22[trans_len_chr22$cds_len > 0,]
# find promoters from different transcripts to overlap
trans_eval <- proms[mcols(proms)$tx_id %in% trans_len_chr22$tx_id]
result = sum(coverage(trans_eval) > 1)
result["chr22"]





