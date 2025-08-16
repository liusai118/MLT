##########################
library(tidyverse) #
library(data.table) 
a1  <- fread('read.count', header = T, data.table = F)
a1 <- na.omit(a1)
counts <- a1[,7:ncol(a1)] 
rownames(counts) <- a1$Geneid 
geneid_efflen <- subset(a1,select = c("Geneid","Length"))
colnames(geneid_efflen) <- c("geneid","efflen")  
geneid_efflen_fc <- geneid_efflen #用于之后比较
dim(geneid_efflen)
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)       
  PMSC_rpk <- sum(RPK)/1e6       
  RPK/PMSC_rpk                    
}  
colnames(counts) <- gsub("_","",colnames(counts))
counts <- as.matrix(counts)

tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)
write.csv(tpm, file = './tpm.csv',row.names = T)
########################
library(ggplot2)
library(ggExtra)
library(Rtsne)
result <- tpm 
group <- c(rep("P7",3),rep("P12",3),rep("P12 + MLT",3))
colnames <- paste0(group,"_",c(1:3,1:6,1:6))
colnames(result) <- colnames
tsne_result <- Rtsne(t(result), 
                     dims = 3,           
                     perplexity = 3,    
                     theta = 0.01,       
                     max_iter = 2000)    
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
rownames(tsne_df) <- colnames


p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group , shape = group )) +
  geom_point(size = 2) +
  labs(title = "t-SNE Projection with Custom Shapes for Groups") +
  scale_shape_manual(values = c(16, 17, 15, 18)) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    # 移除背景线
    panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
p
p <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group , shape = group)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_color_manual(values = c("#91D1C2", "#E64B35", "#4DBBD5")) +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),    # 移除背景线
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  )
p_with_marginal <- ggMarginal(p, type = "density",groupColour = T, groupFill = TRUE)
p_with_marginal

#############
library(limma)
DES <- result
DES <- log(DES +1)
group <- as.matrix(c(rep("p12",3),rep("p12_mlt",3)))
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))

rownames(design) = colnames(DES)
cc<- paste0("p12_mlt", " - ", "p12")
contrast.matrix <- makeContrasts(cc, levels = design) 


fit <- lmFit(DES,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG_ot <- topTable(fit2, adjust.method="fdr", coef=1,sort.by="logFC",n = Inf)
DEG_ot$ID <- rownames(DEG_ot)
library(ggplot2)
library(ggrepel)
DEG_ot <- as.data.frame(DEG_ot)
DEG_ot$log2fc <- DEG_ot$logFC * log2(10)
ggplot(DEG_ot, aes(logFC, -log10(P.Value))) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(P.Value), color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  # geom_text_repel(data = DEG_ot_label , aes(label = DEG_ot_label [,7]), color = "black", size = 3) +
  xlab("Log2FC") +
  ylab("-Log10(p-value)")
#############
library(clusterProfiler)
library(dplyr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
df5 <- DEG_ot %>%
  dplyr::filter(P.Value < 0.05)
ids <- bitr(rownames(df5),'SYMBOL','ENTREZID','org.Hs.eg.db')
ego_ALL <- enrichGO(gene = ids$ENTREZID,#我们上面定义???
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#富集的GO类型
                    pAdjustMethod = "BH",#这个不用管，???般都用的BH
                    minGSSize = 1,
                    pvalueCutoff = 1,#P值可以取0.05
                    qvalueCutoff = 1,
                    readable = TRUE) 
write.csv(ego_ALL,"GO.csv")
library(stringr)
library(ggalluvial)
p5 <- ggplot(ego_ALL@result[1:15,],
             aes(x = Count, y = Description, color = -log(pvalue), size = Count)) +
  geom_point() +
  scale_color_gradientn(
    name = "-log10(pvalue)",
    values = seq(0, 1, 0.5),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  #scale_x_continuous(limits = c(5, 50)) +
  geom_point(shape = 21, color = "black", stroke = 1, alpha = 1) +  # 添加黑色边框
  scale_size_continuous(name = "Count", range = c(1, 8), breaks = c(10, 40, 80, 120, 160)) +  # 调整点的大小范围和刻度，添加图例标题
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
    axis.line = element_line(colour = "black"),  # 添加坐标轴线
    panel.background = element_rect(fill = "white"),  # 添加面板背景
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    axis.title.x = element_blank(),  # 移除 x 轴标题
    axis.title.y = element_blank(),  # 移除 y 轴标题
    axis.text = element_text(family = "Arial", size = 12),  # 设置坐标轴文本字体为 Arial 12 号
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank()  # 去除次要网格线
  )
p5
library(yulab.utils)

