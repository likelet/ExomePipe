# useage: Rscript Maf_analysis.R <Mut.maf> <Sample_info.txt>
args = commandArgs(T)
library("maftools")
library("ggsci")
#args[1] = "../2.MAF/1.ALL_sample_filt/Merged_SNV_INDEL.maf"
#args[2] = "Clinical_info/Sample_info.txt"
if (!is.na(args[2])) {
  laml = read.maf(maf = args[1], clinicalData = args[2])
}else{
  laml = read.maf(maf = args[1])
}
Name = sub(".maf", "", args[1])

flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
            "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
            "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

# Variants summary
print("##### Variants summary #####")
pdf(file = paste(Name, ".summary.pdf", sep=""))
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()
pdf(file = paste(Name, ".tcga_compare.pdf", sep = ""), height=3, width=6)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'Our_Data')
dev.off()
pdf(file = paste(Name, ".titv.pdf", sep = ""))
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
dev.off()
pdf(file = paste(Name, ".gene_cloud.pdf", sep = ""))
geneCloud(input = laml, top = 50)
dev.off()
pdf(file = paste(Name, ".vaf_plot.pdf", sep = ""))
vafPlot = plotVaf(maf = laml, vafCol = 'VAF', flip = TRUE)
dev.off()

# Oncoplot
print("##### Oncoplot #####")
pdf(file = paste(Name, ".oncoplot.pdf", sep = ""), height=8, width=12)
col = ggsci::pal_npg("nrc")(10)
names(col) = c('Frame_Shift_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', 
               'Nonsense_Mutation', 'In_Frame_Del', 'Splice_Site', 'Frame_Shift_Ins',
               'Nonstop_Mutation', 'Translation_Start_Site')
oncoplot(maf = laml, top = 50, removeNonMutated = F, colors = col, fontSize = 8, genesToIgnore = flags)
dev.off()

#Mutational Signatures
print("##### Mutational Signatures #####")
pdf(file = paste(Name, ".Each_signature.pdf", sep = ""))
laml.tnm = trinucleotideMatrix(maf = laml, ref_genome = '/disk/database/human/hg19/genome.fa', 
                               prefix = 'chr', add = TRUE, ignoreChr = 'chr23', useSyn = TRUE)
MS_single = laml.tnm$nmf_matrix
MS_single_fraction = as.matrix(MS_single / rowSums(MS_single))
write.table(t(MS_single_fraction), file = paste(Name, ".MS_fraction.txt", sep = ""), 
            sep = "\t", quote = FALSE)
plotData = MS_single_fraction
nsigs = nrow(plotData)
par(mfrow = c(5, 2), oma = c(5, 4, 0, 0) + 0.1, 
    mar = c(0, 0, 1, 1) + 0.1)
color = c("coral4", "lightcyan4", "deeppink3", 
          "lightsalmon1", "forestgreen", "cornflowerblue")
colors = rep(color, each = 16)
for (i in 1:nsigs) {
  d = as.matrix(t(plotData[i, ]))
  barplot(d, xaxt = "n", yaxt = "n", border = FALSE, 
          col = colors, beside = TRUE, ylim = c(-0.1, 0.2), main = rownames(plotData)[i])
  axis(side = 2, at = c(0, 0.05, 0.1, 0.15, 0.2), 
       labels = c(0, 0.05, 0.1, 0.15, 0.2), pos = c(0, 0.05, 0.1, 0.15, 0.2), las = 2)
  abline(h = c(0.05, 0.1, 0.15, 0.2, 0.25), lty = 2, 
         lwd = 0.3, col = "gray70")
  rect(xleft = seq(0, 192, 32), ybottom = -0.05, 
       xright = 192, ytop = -0.02, col = color, border = "gray70")
  text(labels = c("C>A", "C>G", "C>T", "T>A", "T>C", 
                  "T>G"), y = rep(-0.08, 6), x = seq(0, 192, 32)[2:7] - 16, cex = 0.8)
}
dev.off()

#Signature analysis
print("##### Cismic signature analysis #####")
require('NMF')
pdf(file = paste(Name, ".signature.pdf", sep = ""))
laml.sign = extractSignatures(mat = laml.tnm, nTry = 6, plotBestFitRes = FALSE)
plotSignatures(laml.sign)
plotSignatures(nmfRes = laml.sign, contributions = TRUE)
nsig = ncol(laml.sign$signatures) #number of signatures identified
contrib = t(laml.sign$contributions) #contribution of signatures in every sample
require('ggfortify')
contrib.gg = autoplot(kmeans(contrib, nsig), data = contrib, frame = TRUE, frame.type = 'norm', 
                      loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4, 
                      loadings.colour = 'blue')+cowplot::theme_cowplot(line_size = 1)+
  cowplot::background_grid(major = 'xy')
print(contrib.gg)
require('corrplot')
corrplot::corrplot(corr = laml.sign$coSineSimMat, col = RColorBrewer::brewer.pal(n = 9, name = 'Oranges'), 
                   is.corr = FALSE, tl.cex = 0.6, tl.col = 'black', cl.cex = 0.6)
pheatmap::pheatmap(mat = laml.sign$coSineSimMat, cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures", display_numbers = TRUE)
dev.off()

# Prepare for MutSigCV analysis
print("##### Prepare for MutSigCV analysis #####")
laml.mutsig.corrected = prepareMutSig(maf = laml)
write.table(laml.mutsig.corrected, file = paste(Name, ".for_MuSigCV.maf", sep = ""), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Driver_gene detection (Oncodrive)
print("##### Driver_gene detection (Oncodrive) #####")
pdf(file = paste(Name, ".oncodrive_bubble.pdf", sep = ""))
laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()

# Oncoplot of driver genes
print("##### Oncoplot of driver genes (Oncodrive) #####")
driver_gene = as.character(laml.sig[fdr<0.1]$Hugo_Symbol)
pdf(file = paste(Name, ".oncodrive_oncoplot.pdf", sep = ""), height=8, width=12)
col = ggsci::pal_npg("nrc")(10)
names(col) = c('Frame_Shift_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', 
               'Nonsense_Mutation', 'In_Frame_Del', 'Splice_Site', 'Frame_Shift_Ins',
               'Nonstop_Mutation', 'Translation_Start_Site')
oncoplot(maf = laml, removeNonMutated = F, colors = col, fontSize = 8, genes = driver_gene, genesToIgnore = flags)
dev.off()

# Somatic Interactions of driver genes
print("##### Somatic Interactions of driver genes #####")
pdf(file = paste(Name, ".oncodrive_interaction.pdf", sep = ""))
somaticInteractions(maf = laml, top = 50, pvalue = 0.1, genes = driver_gene)
dev.off()