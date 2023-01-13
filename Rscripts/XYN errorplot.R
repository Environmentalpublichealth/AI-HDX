setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")
data_df <- read.csv("../results/XYN error.csv", header = T)
data_df <- data_df[order(data_df$length),]
data_df$loc <- paste0(data_df$X0, "-", data_df$X1)
data_df <- data_df[!duplicated(data_df$loc),]
data_df <- data_df[,c(4:9)]
names(data_df)[1] <- "Measured"
data_df$Measured <- data_df$Measured/100
data_df$loc <- factor(data_df$loc, levels = data_df$loc)


ggplot(data=data_df, aes(x=loc, y = error))+
  #geom_point(shape = 1)+ geom_line(size =0.6, color="black") + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1)+
  #geom_point(data = data_df, aes(x=loc, y=Measured), color = "red")+
  geom_bar( stat = "identity", alpha = 0.8)+
  geom_point(data = data_df, aes(x = loc, y = error+ 10, size = 0.7, shape = SS)) +
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_x_discrete(labels=data_df$length)+
  scale_shape_manual(labels = c("C,E", "C,H","C,H,E","E"),values = c(15:18))+
  labs(x="Peptide lengths",y="Error %")+
  guides(color = "none", shape = guide_legend(order=1), size = "none")

ggsave("../images/length error XYN.pdf", width = 8, height = 4)
