setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")

# read data file
pred = read.csv("../results/stepRes XYN pred.xlsx - Sheet1.csv", header = T)

# calculate mean and SE for each fragment
library(Rmisc)
library(reshape2)

meltdf <- melt(pred[,c(1,6:15)], id.vars = "X")
datasum <- summarySE(meltdf, measurevar = "value", groupvars = "X")

datasum$col <- 1
datasum$X <- datasum$X + 1

library(ggplot2)
ggplot(datasum, aes(x = col, y = X, fill= value))+
  geom_tile()+
  geom_text(aes(label=round(value, 2)))+
  scale_fill_gradient(low="white",high="#0072B2")+
  theme_classic()+
  scale_y_continuous(trans = "reverse", breaks = datasum$X)+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank())+
  labs(x= "", y = "")
ggsave("../images/stepRes heatmap.pdf", width = 2, height = 8)  
