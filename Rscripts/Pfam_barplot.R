setwd("~/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/Rscripts/")

library(ggplot2)
library(reshape2)
library(plyr)
metadata <- read.table("../all_data.txt", header = F)

# reshape data for a bar plot
table_data <- data.frame(table(metadata$V4))
names(metadata)[4] <- "Var1"
table_data <- join(table_data, metadata[,3:4], by= "Var1", match="first")

fix_data <- edit(table_data) # mannual edit the Pfam have more than one species

# bar plot
fix_data$V3[35] <- "Bacillus"
fix_data$V3 <- factor(fix_data$V3, levels = c("Virus","Ecoli","Rhipicephalus","Salmonella","Drosophila","Klebsiella","Bacillus","Yeast","Human","Myxococcus","Mouse"))
ggplot(fix_data, aes(x=Var1, y=Freq, fill=V3)) + 
  geom_bar(stat="identity") +
  labs(x="Pfam ID", y = "Number of proteins", fill="Organisms")+
  theme_classic()+
  scale_fill_manual(values = c('#713e5a', '#63a375', '#edc79b', '#d57a66', '#ca6680', '#395B50', '#92AFD7', '#b0413e', '#4381c1', '#736ced', '#631a86'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("../images/Barplot_pfam.pdf", height = 3.8, width = 6.8)
