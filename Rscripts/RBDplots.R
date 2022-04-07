setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/")

library("ggplot2")
library("reshape2")
library("plyr")

# import data file
predData <- read.csv("results/sars/RBD_pred0324data.csv", header = T)[1:30,]

data <- melt(predData[,c(1,2,4,6,8)])
data$SD <- c(predData$SD, predData$SD.1, predData$SD.2, predData$SD.3)/10

ggplot(data, aes(x= peptide, y=value, group = variable, color=variable))+
  geom_line(size = 1.1, alpha = 0.9) + 
  geom_errorbar(aes(ymin = value-SD, ymax=value+SD), width=.2)+
  theme_classic(base_size = 16)+
  labs(x = "Peptide fragment", y = "HDX rate", color = "")+
  scale_color_manual(labels=c("Wuhan","Omicron","Wuhan+STE90-C11","Omicron+S309"),
                     values = c("#e6b89c","#ead2ac","#9cafb7","#4281a4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8))

ggsave("images/RBDprediction.pdf", width = 10, height = 5)
