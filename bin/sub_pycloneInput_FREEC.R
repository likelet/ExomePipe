#!/bin/Rscript
options(warn = -1) 
args <- commandArgs(trailingOnly = TRUE)

OutDir=args[1]
patient=args[2]

suppressMessages(library(data.table))
suppressMessages(library(reshape2))

tsv <- read.table(paste(patient, ".tsv", sep=""), header=F, check.names=F, stringsAsFactors = FALSE)
colnames(tsv) <- c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "sample")
datasets <- split(x=tsv, f=tsv$sample)

for (dat in datasets){
	sampleID <- unique(dat[,7])
	absent_mut <- setdiff(tsv[,1], dat[,1])
	absent_df <- data.frame()
	for (mut in absent_mut){
		absent_df <- rbind(absent_df, data.frame(mutation_id=mut, ref_counts=0, var_counts=0, ormal_cn=2, minor_cn=0, major_cn=2, sample=sampleID))
	}
	dat <- rbind(dat,absent_df)
	write.table(dat[,1:6], paste(OutDir, "/", sampleID, ".tsv", sep = ""), quote=F, sep="\t", row.names=F)
}
