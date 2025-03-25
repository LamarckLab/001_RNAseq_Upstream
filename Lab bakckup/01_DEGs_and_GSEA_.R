library(DESeq2)#计算差异基因
library(reshape2)#转数据格式用得到
library(dplyr)#处理数据用的
library(ggplot2)#画图的
library(clusterProfiler) # GSEA富集/基因集读取
library(org.Rn.eg.db)#大鼠基因集，其他物种记得下载了换
library(pheatmap)#画热图

#切换工作路径，读取和保存默认改路径（也可写补充上其他的完整路径）
setwd('G:/转录组3/result/01_DEGs')

load("all_data.Rdata")#重新画图就加载一下保存的数据，第一次跑忽视


###########################计算差异基因并出图######################################
#读入基因表达矩阵与分组信息表
row_count = read.csv("gene_tpm_matrix.csv", sep=",", header=F)#不需要处理行名可以直接T
group = read.table("G_JE_3m.txt", header=T, row.names = 1)

##################处理数据

########先搞行名，保留样本名字（数据有一堆前缀，xxxx--样本名字）
countset_col <- row_count #先赋值一下
countset_col[1, ] <- sub(".*\\--", "", row_count[1, ])#处理一下行名，处理好的行名（样本名字）是要能对应分组样本名字的
colnames(countset_col) <- countset_col[1, ]#第一行作为列名
countset_col <- countset_col[-1, ]#删除第一行

########现在搞列名
countset_row <- countset_col
###去除没有注释到基因ID的以及处理数值类型，
countset1 <- countset_row[grepl("\\|", countset_row[, 1]), ]#去除没有注释到基因id的
countset1[, 2:ncol(countset1)] <- sapply(countset1[, 2:ncol(countset1)], as.numeric)# 将 countset1 中除了第一列之外的所有列转换为数值型

###处理过数值类型才好去除全为0 的行，没注释到基因ID的一般为xxxxx,注释到的为xxxxx|XXXX
countset2 <- countset1[rowSums(countset1[, 2:ncol(countset1)]) > 0, ]#去除全为0的行
countset2[, 1] <- sub(".*\\|", "", countset2[, 1])#将注释到基因id的只保留基因id,去除上游stringtie软件分配的id
countset3 <- aggregate(. ~ gene_id, data = countset2, FUN = sum)#ID一样的合并一下

###一些ID为ENSEMBL，将其统一为SYMBOL（注意有一些太新的还改不了）
#获取的是能改的
count_id <- bitr(countset3[, 1],    
                fromType = "ENSEMBL",#现有的ID类型
                toType = "SYMBOL",#需转换的ID类型
                OrgDb = "org.Rn.eg.db")#大鼠的物种注释包

#可以瞅一眼这几个
cols <- countset3$gene_id[countset3$gene_id %in% count_id$ENSEMBL]
cols #这样 [1] "ENSRNOG00000064192" "ENSRNOG00000067883" "ENSRNOG00000071049" "ENSRNOG00000066643" 
#这一行进行修改

#统一ID格式
matched_symbols <- count_id$SYMBOL[match(countset3$gene_id, count_id$ENSEMBL)]
countset3$gene_id <- ifelse(is.na(matched_symbols), countset3$gene_id, matched_symbols)


#改完就啥都没有了，空
cols <- countset3$gene_id[countset3$gene_id %in% count_id$ENSEMBL]
cols #character(0)


# 将 countset2 数据框保存为 CSV 文件
write.csv(countset3, "gene_tpm_py.csv", row.names = FALSE)



#############################
countset3 = read.csv("gene_count.csv", sep=",", header=T,row.names = 1)#做其他组做就直接读这个
group = read.table("G_JE_6m_vs_3m.txt", header=T, row.names = 1)


countset3<-round(countset3, digits=0)#将输入数据取整,之后将 count 转换为整数类型

########对数据进行分组处理啦！！！！
###将我们count表中的样本删一下，只保留分组表里面的
group_samples <- rownames(group)
cols <- group_samples[rownames(group) %in% colnames(countset3)]
print(cols)#看一眼对不对，打印出来是样本名字
count <- countset3[, cols, drop = FALSE]#去除完成


##################各种计算
coldata=data.frame(state = group$Group, row.names = rownames(group))#选择你要比较的分组
coldata$state <- as.factor(coldata$state)# 将字符型变量转换为因子型
rownames(coldata)=colnames(count)#将样本名与组名对应

#分析差异表达基因
dds=DESeqDataSetFromMatrix(countData = count, colData = coldata, design = ~state)
dds<-dds[rowSums(counts(dds))>1,]
dds<-DESeq(dds)
# 查看结果的名称
resultsNames(dds)

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c(contrast = c('state', 'JE_RC_6m', 'JE_RC_3m')))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
DEG_deseq2 <- na.omit(DEG)

# 输出差异表达基因结果的前几行
head(DEG_deseq2)


######################### 画图 —— 火山图和热图
# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 1
padj = 0.05
k1 <- (DEG_deseq2$padj <= padj) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$padj <=padj) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_deseq2$change)
# down stable up ：看上调下调有几个


# 火山图
p1 <- ggplot(data = DEG_deseq2, 
            aes(x = log2FoldChange, 
                y = -log10(padj))) +
  geom_point(alpha = 0.6, size = 2.5, 
             aes(color = change),stroke = 0) +
  ylab("-log10(Padj)")+
  scale_color_manual(values = c("#386EAA", "grey", "#BD242A"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(padj), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()+
  theme(
    text = element_text(size = 14)  # 设置所有文本的字体大小为12
  )
p1
# 保存图形为PDF文件
ggsave("1_Volcano_Plot.pdf", p1, width = 6, height = 5, dpi = 300)


# 差异基因热图
deg_opt <- DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
count_heatmap <- count %>% filter(rownames(count) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(count_heatmap) 
my_colors <- colorRampPalette(c("#386EAA", "white", "#BD242A"))(100)
#my_colors <- colorRampPalette(c("#2983b5", "white", "#DB311A"))(50)
group_colors <- list(Group=c('JE_RC_6m'="#E69191", 'JE_RC_3m'="#4D8CC3"))#annotation_col中
#count_heatmap <- as.matrix(count_heatmap)

p2 <- pheatmap(count_heatmap, annotation = annotation_col,
               annotation_colors = group_colors,#annotation_col会与这一行重合
               show_colnames = T, show_rownames = F,fontsize_row = 10,
               color = my_colors,
               scale = "row",
               cluster_cols = T,
               breaks = seq(-2.5, 2.5, length.out =100),
               border_color = NA,
               angle_col = 45) 
p2
ggsave("2_Heatmap_Plot_1.pdf", p2, width = 6, height = 5, dpi = 300)



p3 <- pheatmap(count_heatmap, annotation = annotation_col,
               annotation_colors = group_colors,#annotation_col会与这一行重合
               show_colnames = T, show_rownames = T,fontsize_row = 10,
               color = my_colors,
               scale = "row",
               cluster_cols = T,
               breaks = seq(-2.5, 2.5, length.out =100),
               border_color = NA,
               angle_col = 45) 
p3
ggsave("2_Heatmap_Plot_2.pdf", p3, width = 6, height = 350,limitsize = FALSE, dpi = 300)
###############################################原师兄的代码，一次保存差异基因结果##########

#基因过滤函数
DEG_filter<-function(logfc,paj){
  folder_name3 <- paste(folder_name2,"/","logFC",logfc,"_pvalue",paj,sep = "")
  dir.create(folder_name3)
  DEG_AvsB <- subset(resAvsB,abs(resAvsB$log2FoldChange) >= logfc & resAvsB$padj <= paj)
  DEG_AvsB_up <- subset(DEG_AvsB,DEG_AvsB$log2FoldChange > 0 & DEG_AvsB$padj <= paj)
  DEG_AvsB_down <- subset(DEG_AvsB,DEG_AvsB$log2FoldChange < 0 & DEG_AvsB$padj <= paj)
  count_AvsB <- c(nrow(subset(DEG_AvsB,DEG_AvsB$log2FoldChange > 0)),nrow(subset(DEG_AvsB,DEG_AvsB$log2FoldChange < 0)))
  geneID_AvsB <- DEG_AvsB@rownames
  geneID_AvsB_up <- DEG_AvsB_up@rownames
  geneID_AvsB_down <- DEG_AvsB_down@rownames
  
  write.csv(DEG_AvsB,paste(folder_name3, "/", "DEG_", sam1, "_VS_",sam2 ,".csv", sep = ""))
  write.table(geneID_AvsB,paste(folder_name3, "/", "geneID_", sam1, "_VS_",sam2 ,".csv", sep = ""),quote = F,row.names = F,col.names = F)
  write.table(geneID_AvsB_up,paste(folder_name3, "/", "geneID_", sam1, "_VS_",sam2 , "_up", ".csv", sep = ""),quote = F,row.names = F,col.names = F)
  write.table(geneID_AvsB_down,paste(folder_name3, "/", "geneID_", sam1, "_VS_",sam2 , "_down", ".csv", sep = ""),quote = F,row.names = F,col.names = F)
  write.table(count_AvsB, paste(folder_name1, "/", sam1, "_VS_others_counts", ".csv", sep = ""),quote = F,row.names = F,col.names = F, append = TRUE)
  write.table(geneID_AvsB_up, paste(folder_name1, "/", sam1, "_VS_others_all_geneID_up", ".csv", sep = ""),quote = F,row.names = F,col.names = F, append = TRUE)
  write.table(geneID_AvsB_down, paste(folder_name1, "/", sam1, "_VS_others_all_geneID_down", ".csv", sep = ""),quote = F,row.names = F,col.names = F, append = TRUE)
  
}

#设置阈值，可以根据自己的需要进行修改
DEG_compare <- function(sampleA,sampleB){
  pset <- c(0.05)#校正后的p值，一般选取0.05和0.01
  lfset <- c(1,2)#logFC，一般选取1和2
  for(i in pset){
    for(j in lfset){
      DEG_filter(j,i)
    }
  }
}

samples <- levels(dds$state)

#输出差异基因文件
for(sam1 in samples){
  folder_name1 <- paste(sam1,"_DEG",sep = "")
  dir.create(folder_name1)
  for(sam2 in samples){
    if(sam1 != sam2){
      folder_name2 <- paste(folder_name1,"/",sam1,"_VS_",sam2,sep = "")
      dir.create(folder_name2)
      resAvsB<-results(dds, contrast=c("state",sam1,sam2))
      DEG_compare(sam1,sam2)
    }
  }
}
save(list = ls(), file = "all_data.Rdata")
########################GSEA从这里开始###########

# 导入差异分析后的数据，以便后续使用logFC进行基因排序
load("all_data.Rdata")#之前存了一下的数据，从新打开就可以从这里开始跑
head(DEG_deseq2)
length(rownames(DEG_deseq2)) #共13332个基因


#添加entrez ID列：
##symbol转entrez ID：
symbol <- rownames(DEG_deseq2)
entrez <- bitr(symbol,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Rn.eg.db")#大鼠的物种注释
head(entrez)

#准备genelist文件(使用全部基因而非仅差异基因)：
genelist <- DEG_deseq2$log2FoldChange
names(genelist) <- rownames(DEG_deseq2)

#genelist缺失基因过滤(ID转换中缺失的部分)：
genelist <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
head(genelist)

#将genelist按照log2FC值从高到低进行排序：
genelist <- sort(genelist, decreasing = T)
head(genelist)

KEGG_ges <- gseKEGG(
  geneList = genelist,
  organism = "rno",
  minGSSize = 0,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#将ENTREZID重转为symbol：
KEGG_ges <- setReadable(KEGG_ges,
                        OrgDb = org.Rn.eg.db,
                        keyType = "ENTREZID")

#提取结果表并将富集结果保存到本地：
KEGG_ges_result <- KEGG_ges@result
write.csv(KEGG_ges_result, file = c('gsea(KEGG).csv'))


GO_ges <- gseGO(geneList = genelist,
                OrgDb = org.Rn.eg.db,
                ont = "ALL", #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                minGSSize = 0,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                eps = 0,
                verbose = FALSE)
GO_ges <- setReadable(GO_ges, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
GO_ges_result <- GO_ges@result
write.csv(GO_ges_result, file = c('gsea(GO).csv'))

save(list = ls(), file = "all_data.Rdata")
#########################################
# 经典的gsea结果图
library(enrichplot)

gseaplot2(KEGG_ges,1,color="black",pvalue_table = T,title="",base_size=10,ES_geom="line") # 按第一个条目做二维码图，并显示p值#color是enrichment score线的颜#base_sizexy轴标题字的大小

#ES_geom是enrichment score线用线或点显示
gseaplot2(GO_ges,1,pvalue_table = F,title="",base_size=10,ES_geom="line") 

plot_list<-list()
for (i in 1:100) {
  plot_list[[i]]<-gseaplot2(KEGG_ges,i,pvalue_table = F,title = KEGG_ges_result$Description[i],base_size = 10,ES_geom = 'line')
  ggsave(plot_list[[i]],filename=paste0("KEGG_gsea_",KEGG_ges_result$Description[i],".pdf"), width = 10, height = 8)
}

plot_list<-list()
for (i in 1:100) {
  plot_list[[i]]<-gseaplot2(GO_ges,i,pvalue_table = F,title = GO_ges_result$Description[i],base_size = 10,ES_geom = 'line')
  ggsave(plot_list[[i]],filename=paste0("GO_gsea_",GO_ges_result$Description[i],".pdf"), width = 10, height = 8)
}

