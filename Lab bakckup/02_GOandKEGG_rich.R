#### 加载包#######
###1.富集相关包
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(org.Rn.eg.db)
library(patchwork)
library(enrichplot)
library(ggrepel)
library(aPEAR)
library(DOSE)
library(tidyverse)
library(ggsankey)
library(ggplot2)
library(cols4all)

####切换目录
setwd('G:/转录组3/result/02_KEGG_GO_Rich(FC=1.5)/RM_3m_VS_CC_3m/')

#### 加载数据
data <- read.csv(file = 'DEG_RC_9m_VS_RC_3m.csv', header = T)#这个文件在上一个脚本存的文件中找

#### 差异基因提取
#上调
sig_gene_up <- data %>%
  dplyr::filter(padj<0.05 & log2FoldChange > 0)#其实padj都是小于0.05的
ids_up <- bitr(sig_gene_up$X, 'SYMBOL', 'ENTREZID', 'org.Rn.eg.db')   
sig_gene_up <- merge(sig_gene_up, ids_up, by.x='X', by.y='SYMBOL')
gene_diff_up <- sig_gene_up$ENTREZID


#下调
sig_gene_down <- data %>%
  dplyr::filter(padj<0.05 & log2FoldChange < 0)
ids_down <- bitr(sig_gene_down$X, 'SYMBOL', 'ENTREZID', 'org.Rn.eg.db')   
sig_gene_down <- merge(sig_gene_down, ids_down, by.x='X', by.y='SYMBOL')
gene_diff_down <- sig_gene_down$ENTREZID







######### GO富集分析
###UP
GO_UP_ego <- enrichGO(gene_diff_up, OrgDb = "org.Rn.eg.db", qvalueCutoff = 0.05,
                      pvalueCutoff = 0.05, ont="all")
GO_UP_ego <- setReadable(GO_UP_ego, OrgDb = org.Rn.eg.db, keyType="ENTREZID")#转换ID

###DOWN
GO_DOWN_ego <- enrichGO(gene_diff_down, OrgDb = "org.Rn.eg.db",qvalueCutoff = 0.05,
                        pvalueCutoff = 0.05,ont="all")
GO_DOWN_ego <- setReadable(GO_DOWN_ego, OrgDb = org.Rn.eg.db, keyType="ENTREZID")#转换ID

#######KEGG分析
###UP
KEGG_UP_ego <- enrichKEGG(gene_diff_up, organism = 'rno', qvalueCutoff = 0.05)
KEGG_UP_ego <- setReadable(KEGG_UP_ego, OrgDb = org.Rn.eg.db, keyType="ENTREZID")#转换ID
###DOWN
KEGG_DOWN_ego <- enrichKEGG(gene_diff_down, organism = 'rno', qvalueCutoff = 0.05)
KEGG_DOWN_ego <- setReadable(KEGG_DOWN_ego, OrgDb = org.Rn.eg.db, keyType="ENTREZID")#转换ID


#####保存结果
write.csv(GO_UP_ego, file = './GO_UP.csv')
write.csv(GO_DOWN_ego, file = './GO_DOWN.csv')
write.csv(KEGG_UP_ego, file = './KEGG_UP.csv')
write.csv(KEGG_DOWN_ego, file = './KEGG_DOWN.csv')

save(list = ls(), file = "rich_data.Rdata")


##########画图#######

###第二次跑读取
load("rich_data.Rdata")


####循环画图必不可少
objects <- list(GO_UP = GO_UP_ego, GO_DOWN = GO_DOWN_ego, 
                KEGG_UP = KEGG_UP_ego, KEGG_DOWN = KEGG_DOWN_ego)

########### 1. 气泡图
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  result <- ego@result[ego@result$p.adjust < 0.05, ]
  
  if (nrow(result) == 0) {
    cat("Skipping object", f, "because result is empty\n")
    next
  }
  
  #merge_data1 <- subset(merge_data1, ONTOLOGY == "BP")# subset 函数用于选择符合条件的行,BP,CC,MF
  result$GeneRatio <- sapply(strsplit(result$GeneRatio, "/"),# 确保merge_data2$GeneRatio中的所有值都是字符类型 
                           function(x) as.numeric(x[1])/as.numeric(x[2]))
  result$Description <- factor(result$Description, # 将 y 列转换为有序的因子变量
                             levels = result$Description)
  result <- head(result[order(result$Count, decreasing = TRUE), ], 15)# 排序并选择前x行
  multiplier1 <- 0.95
  multiplier2 <- 1.05
  p1 <- ggplot(data = result, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                                size = Count, color = p.adjust)) +
        geom_point() +
        theme_bw() +
        theme(axis.text = element_text(color = "black", size = 13),
              axis.title.y = element_blank(),
              aspect.ratio = 2.5,# 设置绘图区域的宽高比例为2.5
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +  
        scale_x_continuous(labels = scales::number_format()) +
        labs(x = "GeneRatio", y = "") +# 设置y轴标签为空
        guides(size = guide_legend(order = 1))+#设置一下图例顺序
        scale_size(range = c(2, 8))+# 设置点的大小范围在2到8之间
        coord_cartesian(xlim = c(min(result$GeneRatio) * multiplier1, #调整坐标的，不然两边的点会显示不全
                                 max(result$GeneRatio) * multiplier2)) 
        #labs(size = "Count", color = "Q Value")# 如果需要修改图例标签，可以使用labs函数
        #scale_color_continuous_c4a_seq(palette = 'yl_gn_bu')#一个配色方案，好看但是有些点颜色太浅了看不起，弃用

  ggsave(filename = paste0("01_Bubble_", f, ".pdf"),  p1, width = 10, height = 4.4, dpi = 300)
  print(paste0("01_Bubble_", f, ".pdf")) # 检查文件名
}

########## 2.关联图（出的图不好看）且懒得优化（clusterProfiler包自带的）###弃用
geneList <- data$log2FoldChange
names(geneList) <- data$X
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  p2 <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 10)
  p3 <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 20, circular = TRUE)
  p4 <- cnetplot(ego, color.params = list(foldChange = geneList, edge = FALSE), showCategory = 3, node_label = "gene")
  p5 <- cnetplot(ego, color.params = list(foldChange = geneList), showCategory = 3, layout = "gem")
  #scale_fill_gradient2(low = "#386EAA", mid="white",high = "#BD242A") # 手动指定离散变量对应的颜色
  ggsave(filename = paste0("Cir_", f, ".pdf"), (p2 + p3) / (p4 + p5), width = 20, height = 20, dpi = 300)
  print(paste0("Cir_", f, ".pdf")) # 检查文件名
}


########## 3.热图（clusterProfiler包自带的）
geneList <- data$log2FoldChange
names(geneList) <- data$X
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  result <- ego@result[ego@result$p.adjust < 0.05, ]
  if (nrow(result) == 0) {
    cat("Skipping object", f, "because result is empty\n")
    next
  }
  
  p6 <- heatplot(ego, foldChange=geneList, showCategory = 10)+#展示前10项
        scale_y_discrete(labels=function(x) str_wrap(x, width = 150))+
        theme(axis.text = element_text(color = "black"))
  #scale_fill_gradient2(low = "#386EAA", mid="white",high = "#BD242A") # 手动指定离散变量对应的颜色
  ggsave(filename = paste0("02_Heat_", f, ".pdf"), p6, width = 10, height = 3, dpi = 300)
  print(paste0("02_Heat_", f, ".pdf")) # 检查文件名
}


######### 4.富集通路之间的网络关系图（clusterProfiler包自带的）#不是很好看弃用
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  result <- ego@result[ego@result$p.adjust < 0.05, ]
  if (nrow(result) == 0) {
    cat("Skipping object", f, "because result is empty\n")
    next
  }
  res1 <- enrichplot::pairwise_termsim(ego)
  p7 <- emapplot(res1, split = "ONTOLOGY", showCategory = 30, 
                 layout.params = list(layout = "kk"), #circle，kk:布局方式
                 edge.params = list(min = 0.7), #富集到两个通路的基因交集除并集为0.7
                 cex.params = list(category_node = 1.5)) +
        theme_minimal() +#改一些主题样式的
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                panel.grid = element_blank()) # 隐藏网格线
  ggsave(filename = paste0("Net_", f, ".pdf"), p7, width = 6, height = 6, dpi = 300)
  print(paste0("Net_", f, ".pdf")) # 检查文件名
}


######### 4.富集通路之间的网络关系图_2（aPEAR包，具体可看相应文章说明）
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  result <- ego@result[ego@result$p.adjust < 0.05, ]
  if (nrow(result) == 0) {
    cat("Skipping object", f, "because result is empty\n")
    next
  }

  set.seed(123)
  p8 <- enrichmentNetwork(result,
                          colorBy = 'p.adjust',
                          colorType = c("pval"),
                          nodeSize = 'Count',
                          fontSize = 3,#设置节点标签的字体大小
                          pCutoff = -10,
                          drawEllipses = T, #是否绘制置信区间椭圆
                          verbose = TRUE)+
        guides(size = guide_legend(order = 1))+#将圆圈图例作为第一图例（就是放上方）
        labs(color = "-log10(p,adjust)")+#调整颜色图例的名称
        theme(legend.key = element_blank())+#去除圆圈图例上的方框主题
        scale_color_continuous_c4a_seq('yl_gn_bu')#配色
  ggsave(filename = paste0("03_Net_", f, ".pdf"), p8, width = 8, height = 8, dpi = 300)
  print(paste0("03_Net_", f, ".pdf")) # 检查文件名
}


###########4.桑基图
for (i in seq_along(objects)) {
  f <- names(objects)[i]
  ego <- objects[[i]]
  merge_data1 <- ego@result[ego@result$p.adjust < 0.05, ]
  if (nrow(merge_data1) == 0) {
    cat("Skipping object", f, "because result is empty\n")
    next
  }
  merge_data1 <- head(merge_data1[order(merge_data1$Count, decreasing = TRUE), ], 10)# 排序并选择前x行
  #提取数据，整理为画图格式
 test<-list()
  for (i in 1:nrow(merge_data1)) {
    test[[merge_data1$Description[i]]]<-data.frame(Description=rep(merge_data1$Description[i],length(strsplit(merge_data1$geneID[i],'/')[[1]])),
                                                gene=rep('null',length(strsplit(merge_data1$geneID[i],'/')[[1]])))
    for (j in 1:length(strsplit(merge_data1$geneID[i],'/')[[1]])) {
      test[[merge_data1$Description[i]]]$gene[j]<-strsplit(merge_data1$geneID[i],'/')[[1]][j]
    }
  }
  test_2<-data.frame(Description=c('null'),gene=c('null'))
  for (i in names(test)) {
    test_2<-rbind(test_2,test[[i]])
  }
  test_2 <- test_2[-1, ]# 删除第一行
  merge_data1<- merge(test_2, data, by.x = "gene", by.y = "X", all.x = TRUE)# 使用merge函数将两个数据框按照gene和X列进行合并
# 根据Description列匹配合并数据框，并只保留 result 数据框中的 BgRatio 列
 merge_data1<- merge(merge_data1,ego@result[, c("Description", "GeneRatio")], by = "Description", all.x = TRUE)
# 确保merge_data2$GeneRatio中的所有值都是字符类型
  merge_data1$GeneRatio <- sapply(strsplit(merge_data1$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  dt <- merge_data1[, c("gene", "Description", "log2FoldChange")]#本地数据读入：
  df <- dt %>%
  make_long(gene, Description)#将数据转换为绘图所需格式:
  df$node <- factor(df$node,levels = c(dt$gene %>% unique()%>% rev(),#指定因子，调整显示顺序：
                                     dt$Description %>% unique() %>% rev()))
  unique_count <- length(unique(merge_data1$gene)) + length(unique(merge_data1$Description))# 计算需要的颜色个数
  mycol <- c4a('rainbow_wh_rd',unique_count)
  p9 <- ggplot(df, aes(x = x,next_x = next_x,
                      node = node,next_node = next_node,
                      fill = node,label = node)) +
  scale_fill_manual(values = mycol) + #更改配色
  geom_sankey(flow.alpha = 0.5, #条带不透明度
              smooth = 8, #条带弯曲度
              width = 0.12) + #节点宽度
  geom_sankey_text(size = 3.2,
                   color = "black") +
  theme_void() +
  theme(legend.position = 'none')

  ggsave(filename = paste0("04_Sankey_", f, ".pdf"), p9, width = 7, height = 7, dpi = 300)
  print(paste0("04_Sankey_", f, ".pdf")) # 检查文件名
}

###########5.弦图 
library(circlize)
#一个个跑，循环会报错，没有得到解决
ego <- GO_UP_ego
ego <- GO_DOWN_ego
ego <- KEGG_UP_ego
ego <- KEGG_DOWN_ego
  merge_data1 <- ego@result[ego@result$p.adjust < 0.05, ]
  merge_data1 <- head(merge_data1[order(merge_data1$Count, decreasing = TRUE), ], 10)
  #提取数据，整理为画图格式
  test<-list()
  for (i in 1:nrow(merge_data1)) {
    test[[merge_data1$Description[i]]]<-data.frame(Description=rep(merge_data1$Description[i],length(strsplit(merge_data1$geneID[i],'/')[[1]])),
                                                   gene=rep('null',length(strsplit(merge_data1$geneID[i],'/')[[1]])))
    for (j in 1:length(strsplit(merge_data1$geneID[i],'/')[[1]])) {
      test[[merge_data1$Description[i]]]$gene[j]<-strsplit(merge_data1$geneID[i],'/')[[1]][j]
    }
  }
  test_2<-data.frame(Description=c('null'),gene=c('null'))
  for (i in names(test)) {
    test_2<-rbind(test_2,test[[i]])
  }
  test_2 <- test_2[-1, ]# 删除第一行
  merge_data1<- merge(test_2, data, by.x = "gene", by.y = "X", all.x = TRUE)# 使用merge函数将两个数据框按照gene和X列进行合并
  
  all_genes <- unique(merge_data1$gene)
  all_pathway <- unique(merge_data1$Description)
  
  # 创建一个空的数据框，行数为全集基因ID的数量,列数为全集pathway的数量
  sets <- data.frame(matrix(ncol = length(all_pathway), nrow = length(all_genes)))
  row.names(sets) <- all_genes
  colnames(sets) <- all_pathway
  
  # 将数据框中的所有值初始化为0
  sets[] <- 0
  # 填充数据
  for (gene in all_genes) {
    for (pathway in all_pathway) {
      matching_row <- merge_data1[merge_data1$gene == gene & merge_data1$Description == pathway, ]# 查找merge_data1中基因和描述匹配的行
      if (nrow(matching_row) > 0) { # 检查是否找到匹配的行
        value <- matching_row$log2FoldChange
        sets[gene, pathway] <- value# 如果找到匹配的行，获取数据并填充到sets中
      }
    }
  }
  sets <- sets[rowSums(sets[, 1:ncol(sets)]) != 0, ]
  set1 <- as.matrix(sets)
  
  #设定弦图颜色
  
  unique_count1 <- length(unique(merge_data1$Description))
  unique_count2 <- length(unique(merge_data1$gene)) 
  grid.col = NULL#设定颜色
  grid.col[colnames(sets)] = c4a('brewer.set3',unique_count1) 
  grid.col[rownames(sets)] = c4a('rainbow_wh_rd',unique_count2)
  
  # 设置文本颜色
  grid.col2 = NULL#设定颜色
  grid.col2[colnames(sets)] = c("transparent")#把列名改为无色。后面加图例了
  grid.col2[rownames(sets)] = c("black")


pdf("05_Chordal_GO_DOWN.pdf", width =8, height = 4)####挨个存，写循环后面会报错
 
par(mar = c(4, 1, 1, 16))  # 下、左、上、右的边距值

chordDiagram(
  set1,
  diffHeight = 0.06,
  annotationTrack = c("grid"),
  annotationTrackHeight = c(0.1, 10),
  grid.col = grid.col,
  link.lwd = 0.02,
  transparency = 0.5
)
#文本注释颜色
for(si in get.all.sector.index()) {
  myCol <- grid.col2[si]
  xlim = get.cell.meta.data("xlim",sector.index = si,track.index = 1)
  ylim = get.cell.meta.data("ylim",sector.index = si,track.index = 1)
  circos.text(mean(xlim), ylim[2],labels = si,sector.index = si,
              track.index = 1, facing = "clockwise", col = myCol,
              cex=0.8,adj=c(0,.5),niceFacing = T,xpd=T)#cex控制字体大小
}

#设定图例
legend(x = 1.2, y=0.6,                      # 调整图例位置，可能需要根据实际情况微调y值
       pch = 22,                               # 图例中点的形状（方形）
       title="Pathway ",title.adj=0 ,         #标题名字与位置，0为左对齐1为右对齐
       bty='n',                                #n为不加外边框，o为添加外边框
       legend = colnames(set1),                # 图例中的标签
       col = grid.col[colnames(set1)],         # 图例中点的颜色
       pt.bg = grid.col[colnames(set1)],       # 图例中点的填充颜色，与颜色相同
       cex = 0.8,                                # 图例文本的大小
       pt.cex = 2.5,                             # 图例中点的大小
       border = "black",                       # 图例中点的边框颜色
       ncol = 1,xpd=T) ###一列，xpd表示图例是否可以超过图本身边界绘制，与画布par(mar=c())同用。

dev.off()















