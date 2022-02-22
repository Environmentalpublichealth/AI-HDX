setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")

data <- read.csv("../results/1xyn.csv", header = T)

cor.test(data$Measured, data$Predicted, method = "spearman", exact = F)

eg <- read.csv("../results/eg.csv", header = T)

cor.test(eg$Measured, eg$Predicted, method = "pearson", exact = F)
