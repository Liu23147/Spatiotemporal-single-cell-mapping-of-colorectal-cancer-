##############
Cellratio <- prop.table(table(Idents(scRNA1), scRNA1$tissue), margin = 2)
Cellratio

Cellratio <- as.data.frame(Cellratio)

library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+ scale_color_manual(
    values = zzm8colors
  )+ scale_fill_manual(
    values = zzm8colors 
  )

Cellratio <- prop.table(table(scRNA1$tissue, scRNA1$sample), margin = 2)
Cellratio1 <- prop.table(table(Idents(scRNA1), scRNA1$sample), margin = 2)

Cellratio <- as.data.frame(Cellratio)
Cellratio1 <- as.data.frame(Cellratio1)

merged_data <- left_join(Cellratio1, Cellratio, by = c("Var2" = "Var2"))
merged_data  <- filter(merged_data, Freq.y != 0)
merged_data$Freq <- as.numeric(merged_data$Freq)
ggplot(merged_data, aes(x = Var1.x, y = Freq.x, color =  Var1.y)) +
  geom_boxplot(aes(color = Var1.y), alpha = 0.3, outlier.shape  = NA, width = 0.6) + 
  geom_jitter(position = position_jitterdodge(jitter.width  = 0.1, dodge.width  = 0.6), 
              size = 1.5, alpha = 0.6) + 
  theme_minimal(base_size = 15) + 
  theme(panel.grid  = element_blank(), 
        axis.line  = element_line(color = "black"), 
        axis.ticks  = element_line(color = "black")) + 
  labs(x = "Cell Type", y = "Cell Proportion", color = "Sample") + 
  ggtitle("Cell Proportion") + 
  stat_compare_means(aes(label = ..p.signif..),  label.x = 1.5) + 
  scale_color_manual(values = c("HI" = "#4991C1", "ANT" = "#79B99D", "Tumor" = "#C6307C"))
