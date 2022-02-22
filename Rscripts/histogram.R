setwd("~/Desktop/Jiali/TAMU/Dynamics/")

data <- read.csv("XYN1_pred.csv", header = T)

data$diff <- data$Predicted - data$Measured

library(ggplot2)

ggplot(data, aes(x=diff)) +
  geom_histogram(aes(y=..density..), color = "black", fill = "#3776ab", bins = 8) +
  geom_density(alpha=.4, fill = "#3776ab")+
  theme_classic(base_size = 16)+
  labs(x= "Distance between predicted and measured HDX (%)", y = "Frequency")
ggsave("Seq2HDX/images/Xyn_pred_distance.pdf", width = 6, height = 4)
