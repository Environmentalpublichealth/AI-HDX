setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")
library(ggplot2)

# dataset 
CI_data <- read.csv("../results/validationCI.csv", header = T)
ggplot(CI_data, aes(x = mean.CI, y = RMSE, label= ID))+
  geom_point()+
  geom_text(hjust = 0, vjust=1.3)+
  # geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept = 0.6, linetype="dotted")+
  # ylim(-0.25, 0.1)+
  xlim(0.42,0.7)+
  labs(x = "mean CI")+
  theme_classic()

ggsave("../images/validateCI.pdf", height = 6, width = 5.5)

# barplot
ggplot(CI_data)+
  geom_bar(aes(x = ID, y = RMSE), stat = "identity")+
  geom_line(aes(x = ID, y = mean.CI), group=1)+
  geom_point(aes(x= ID, y = mean.CI))+
  # geom_hline(yintercept=0, linetype="dashed")+
  geom_hline(yintercept = 0.6, linetype="dotted")+
  # ylim(-0.25, 0.1)+
  theme_classic()

# CI for two testing proteins peptide

