args = commandArgs(T)
library("facets")
set.seed(1234)

# compute and fit models 
rcmat = readSnpMatrix(args[1])
xx = preProcSample(rcmat)
oo=procSample(xx,cval=150)
fit=emcncf(oo)

# get purity 
pur.str <- paste(args[2], " Purity:", fit$purity, sep = "\t")

# write purity into file 
filecon=file(paste(args[2], ".purity.txt", sep = ""))
writeLines(pur.str,filecon)
close(filecon)


# write cnv and LOH result 
pdf(file = paste(args[2], ".pdf", sep = ""))
plotSample(x=oo,emfit=fit)
out.loh = fit$cncf[which(fit$cncf$tcn>0 & fit$cncf$lcn==0 & abs(fit$cncf$mafR) > fit$purity*fit$purity),]
output = data.frame("chrom" = out.loh$chrom, "start" = out.loh$start, "end" = out.loh$end, 
                    "nhet" = out.loh$nhet, "mafR" = out.loh$mafR, "mafR.clust" = out.loh$mafR.clust, 
                    "cnlr.media" = out.loh$cnlr.median, "cnlr.median.clust" = out.loh$cnlr.median.clust,
                    "tcn.em" = out.loh$tcn.em, "lcn.em" = out.loh$lcn.em)
if (length(rownames(output)) > 0){
  output$chrom = paste("chr", output$chrom, sep = "")
}
write.table(output, file = paste(args[2], ".loh.strict.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
out.loh = fit$cncf[which(fit$cncf$tcn>0 & fit$cncf$lcn==0),]
output = data.frame("chrom" = out.loh$chrom, "start" = out.loh$start, "end" = out.loh$end, 
                    "nhet" = out.loh$nhet, "mafR" = out.loh$mafR, "mafR.clust" = out.loh$mafR.clust, 
                    "cnlr.media" = out.loh$cnlr.median, "cnlr.median.clust" = out.loh$cnlr.median.clust,
                    "tcn.em" = out.loh$tcn.em, "lcn.em" = out.loh$lcn.em)
if (length(rownames(output)) > 0){
  output$chrom = paste("chr", output$chrom, sep = "")
}
write.table(output, file = paste(args[2], ".loh.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
dev.off()
