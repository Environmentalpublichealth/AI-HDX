setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")

library(pheatmap)
library(RColorBrewer)
library(WGCNA)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
freqData <- read.csv("../results/AA freq unnormalized.csv", header = T, row.names = 1)
pheatmap(freqData, cluster_rows = F, cluster_cols = F, color = rev(colors))
xlabel = c("0-20%","20-40%","40-60%","60-80%","80-100%")

labeledHeatmap(Matrix = freqData,
               xLabels = xlabel,
               yLabels = rownames(freqData),
               ySymbols = rownames(freqData),
               colorLabels = FALSE,
               colors = blues9,
               textMatrix = round(freqData,3),
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(0,2),
               main = paste("Amino Acid Frequency"))
