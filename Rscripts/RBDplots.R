setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/")

library("ggplot2")
library("reshape2")
library("plyr")

# import data file
predData <- read.csv("results/sars/RBD_pred1206data.csv", header = T)[1:30,]

data <- melt(predData[,c(1,2,4,6,8,10,12)])
data$SD <- c(predData$SD, predData$SD.1, predData$SD.2, predData$SD.3, predData$SD.4)/10
#colorscale = c("#e6b89c","#ead2ac","#9cafb7","#4281a4")
colorscale <- c("#554cc8","#7bd38c","#e5c54b","#F19417","#545151")

ggplot(data)+
  geom_line(aes(x= peptide, y=value, group = variable, color=variable, linetype = variable), size = 1.3, alpha = 0.8) + 
  geom_errorbar(aes(x = peptide, y = value, color=variable, ymin = value-SD, ymax=value+SD), width=.2)+
  geom_point(data = predData, aes(x = peptide, y = 0.8, size = 1, shape = SS))+
  theme_classic(base_size = 16)+
  labs(x = "Peptide fragment", y = "HDX rate", color = "", linetype="")+
  scale_color_manual(labels=c("Wuhan","Omicron","Wuhan+STE90-C11","Omicron+S309", "Wuhan+S309"),
  values = colorscale)+
  scale_linetype_manual(labels=c("Wuhan","Omicron","Wuhan+STE90-C11","Omicron+S309", "Wuhan+S309"),
                        values = c("solid","solid","solid","solid","dashed"))+
  scale_shape_manual(name = "Secondary Structure", labels = c("C","C,E","C,H","E","H"), values = c(14:16, 18, 17))+
  scale_size(guide = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.8), legend.position="top", legend.box = "vertical",
        legend.title=element_text(size=12))

ggsave("images/RBDprediction1206.pdf", width = 9, height = 6)
ggsave("images/RBDprediction1206.png", width = 9, height = 6)
