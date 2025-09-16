library(data.table)
data2 <- fread("remap2022_EP300_nr_macs2_hg38_v1_0.bed.gz")
####################
library(GenomicFeatures)
library(Homo.sapiens)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(clusterProfiler)
library(dplyr)
genome <- BSgenome.Hsapiens.UCSC.hg38
data <- data2[data2$V5 > 6,c(1,7)]
data$end <- data$V7 + 100
data$start <- data$V7 - 100
data <- data[,c(1,3,4)]
colnames(data) <- c("seqnames","end","start")
t1 <- as.data.frame(table(data2$V5))
write.table(data[,c(1,3,2)], file = "ep300.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
bed_file <- import("ep300.bed", format = "BED")
sequences <- getSeq(genome, bed_file)
write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  fasta_lines <- character()
  for (i in seq_along(sequences)) {
    chr <- bed_file[i, "seqnames"]
    start <- bed_file[i, "start"]
    end <- bed_file[i, "end"]
    seq_name <- paste0(chr, ":", start, "-", end)
    fasta_lines <- c(fasta_lines, paste0(">", seq_name), as.character(sequences[[i]]))
  }
  writeLines(fasta_lines, con = filepath)
}
write_sequences_to_fasta(sequences,data, filepath = "EP300.fa")
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gkmSVM)
neg2 <- read.table("ep300_neg.fa")
nchar(neg2$V1[4]) 
genNullSeqs('ep300.bed',nMaxTrials=180,xfold=30,genome=BSgenome.Hsapiens.UCSC.hg38.masked,outputPosFastaFN='n1.fa', outputBedFN='neg1x_n3.bed', outputNegFastaFN='neg1x_n3.fa')
####################

