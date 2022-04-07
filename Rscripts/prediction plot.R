setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")

# read the prediction results
data_df <- read.csv("../results/Q92731 pred3.xlsx - Sheet1.csv", header = T, row.names = 1)
data_df$loc <- paste0(data_df$X0, "-", data_df$X1)
data_df <- data_df[,c(4:14)] # XYN use 4:19, ER use 4:14
names(data_df)[1] <- "Measured"
data_df$Measured <- data_df$Measured/100

# XYN data has duplicated peptides, we will remove the duplication to reduce the confusion
data_df <- data_df[!duplicated(data_df$loc),]
# reshape data table 
library(reshape2)
library(Rmisc)
melt_df <- melt(data_df[,-c(1,9)], id.vars = c("loc", "confidence","CV.","SS"))
datasum <- summarySE(melt_df, measurevar = "value", groupvars = c("loc","confidence","CV.","SS"))

# draw plot
library(ggplot2)
ggplot(data = melt_df, aes(x = loc, y = value, group=variable))+
  geom_line(aes(color = variable))+
  geom_point(aes(color = variable, alpha = confidence, size = confidence, shape=SS))+ theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  labs(x="Peptides", y="HDX rate", color="Models",size="Confidence") +
  guides(alpha = "none")+
  scale_color_manual(values = c("red","#999999", "#E69F00", "#56B4E9","#D55E00","#CC79A7","#0072B2"))

ggsave("../images/HDX_plot.pdf", height = 4.1, width = 8)  

# draw plot with error bar
datasum$loc <- factor(datasum$loc, levels = data_df$loc)
data_df$loc <- factor(data_df$loc, levels = data_df$loc)
ggplot(data=datasum, aes(x=loc, y = value, group=1))+
  geom_point(aes(alpha = confidence, size = 2, shape = SS))+ geom_line(size =0.6, color="black") + geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.1)+
  geom_point(data = data_df, aes(x=loc, y=Measured), color = "red")+
  geom_line(data = data_df, aes(x=loc, y=Measured), color = "red", size = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_shape_manual(values = c(15:18))+
  labs(x="Peptides",y="HDX rate",alpha = "Confidence")+
  guides(color = "none", shape = guide_legend(order=1), size = "none", alpha = guide_legend(order=2))

ggsave("../images/XYN_plot_errorbar.pdf", height = 4.1, width = 8)
# plot histogram of the prediction diff
data_df$dff <- data_df$Measured - data_df$average
data_filter <- data_df[which(data_df$Measured > 0.2 & data_df$Measured < 0.7),]
ggplot(data_filter, aes(x=dff)) +
  geom_histogram(aes(y=..density..), color = "black", fill = "#0072B2", bins = 8) +
  geom_density(alpha=.4, fill = "#0072B2")+
  theme_classic(base_size = 16)+
  labs(x= "Distance between predicted and measured HDX", y = "Frequency")
ggsave("../images/ER_pred_distance0405.pdf", width = 6, height = 5)



# calculate the correlation coeffient
cor.test(data_df$Measured, data_df$average, method = "spearman", exact = F)
# learn where are these high confidence peptides from!
ggplot(data = data_df, aes(x=Measured, y=average))+
  geom_point(aes(color = CV., shape = SS), size = 4)+
  scale_color_gradient(low = "#0072B2",high = "black")+
  scale_shape_manual(values = c(15:18))+
  theme_classic()+
  geom_abline(intercept = 0, slope = 1)+
  guides(color = guide_colorbar(order=1),
         shape = guide_legend(order=2))
ggsave("../images/XYN_dots.pdf", height = 4, width = 5)
