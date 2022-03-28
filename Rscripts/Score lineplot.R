setwd("~/Desktop/Jiali/TAMu/Dynamics/Seq2HDX/Rscripts/")

# read data file
df_data <- read.csv("../results/1xyn_alphaFold_conf.csv", header = T)

library(reshape2)
library(ggplot2)

df_conf <- df_data[,c(1,2,5,6)]

df_plot <- melt(df_conf, id = "peptides")

# plot line with sorted confidence score
x_order <- df_conf$peptides[order(df_conf$error, decreasing = FALSE)]
df_plot$peptides <- factor(df_plot$peptides, levels = x_order)
ggplot(df_plot, aes(x = peptides, y = value, group = variable, color = variable))+
  geom_line(stat = "identity", aes(color = variable))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = c("#0072B2",'#66CC99', 'black'), name = "", labels = c("AI-HDX","AlphaFold","Error"))+
  labs(y = "Score")
ggsave("../images/AlphaFold Conf.pdf", height = 2.5, width = 9)

# confidence index vs error

