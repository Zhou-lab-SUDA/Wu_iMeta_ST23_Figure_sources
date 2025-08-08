
#KEgg富集分析
library(clusterProfiler)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(ggrepe1)

enrich = read.csv(file = "to_enrich.csv",header = TRUE)

id_list <- mapIds(org.Hs.eg.db,enrich$gene,"ENTREZID","SYMBOL")

#GO富集函数
go <- enrichGO(gene = id_list, # Entrez ID列表
               OrgDb = org.Hs.eg.db, # 指定物种数据库
               keyType = "ENTREZID", # 指定给定的名称类型
               ont = "ALL", # 可选,BP(生物学过程)/CC(细胞组分)/MF(分子功能)/ALL(同时指定)
               pAdjustMethod = "none", # P值校正方法,还可以是fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q值阈值
               readable = T # 将ID转换为symbol
)

go.res <- as.data.frame(go)
write.csv(go.res,file = "go.csv")

#绘制GO富集分析条形图，结果默认按qvalue排序，分别筛选出前十的term进行绘图即可
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:10,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
#绘图
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
go_bar
ggsave("GO.pdf",go_bar,wi=10,he=7)


#KEGG富集分析
kegg <- enrichKEGG(gene = id_list, 
                   organism = "hsa",keyType = "kegg", 
                   pAdjustMethod = "none",pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                   minGSSize = 10,maxGSSize = 500,use_internal_data = F)

kegg = setReadable(kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")

kegg = data.frame(kegg)
kegg$Description = factor(kegg$Description,levels = rev(kegg$Description))


kegg$ratio = parse_ratio(kegg$GeneRatio)

write.csv(kegg,file = "kegg.csv")
kegg <- kegg[1:10,]

#绘图
p = ggplot(data = kegg,aes(x = ratio, y = reorder(Description,Count)))+
  geom_point(aes(size = Count,color = -log10(p.adjust)))+
  theme_bw()+
  scale_colour_gradient(low = "blue",high = "red")+
  scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
  labs(x = "GeneRatio",y = "",title = "Dotplot",
       color = expression(-log10(p.adjust)),size = "Count")+
  theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 11),legend.text = element_text(size = 10))

p
ggsave("KEGG.pdf",p,wi=10,he=7)


