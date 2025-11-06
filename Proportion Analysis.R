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

library(ggplot2)
library(ggcorrplot)
library(reshape2)
Cellratio <- prop.table(table( scRNA1$sample,Idents(scRNA1)), margin = 2)
Cellratio <- as.data.frame(Cellratio)
matrix_data <- reshape2::dcast(Cellratio , Var1 ~ Var2, value.var = "Freq", fill = 0)
rownames(matrix_data) <- matrix_data$Var2
matrix_data <- matrix_data[, -1]

cor_res <- Hmisc::rcorr(as.matrix(matrix_data))

correlation_matrix <- cor_res$r
p_value_matrix <- cor_res$P

cor_long <- melt(correlation_matrix)
p_long <- melt(p_value_matrix)
colnames(cor_long) <- c("Var1", "Var2", "Correlation")
colnames(p_long) <- c("Var1", "Var2", "P_value")
cor_p_data <- merge(cor_long, p_long, by = c("Var1", "Var2"))
cor_p_data$Significance <- ifelse(cor_p_data$P_value < 0.001, "***",
                                  ifelse(cor_p_data$P_value < 0.01, "**",
                                         ifelse(cor_p_data$P_value < 0.05, "*", "")))
ggplot(cor_p_data, aes(Var1, Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  geom_text(aes(label = Significance), color = "black", size = 4) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 14), 
        axis.text.y = element_text(face = "bold", size = 14), 
        axis.title = element_text(face = "bold", size = 16)) + 
  labs(x = "Cell Types", y = "Cell Types", title = "Correlation Heatmap with Significance") +
  scale_fill_gradient2(low = "#1f77b4", high = "#d62728", mid = "#ffffff", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation")















