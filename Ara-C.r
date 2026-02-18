library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

gene_exp_with_symbol <- gene_exp_with_symbol %>% 
  mutate(sig=case_when(G2_vs_G1_log2FoldChange <=  -1 & G2_vs_G1_FDR < 0.05 ~ 'down',
                       G2_vs_G1_log2FoldChange >=  1 & G2_vs_G1_FDR < 0.05 ~ 'up',
                       G2_vs_G1_FDR > 0.05 ~ 'non',
                       G2_vs_G1_log2FoldChange > -1 | G2_vs_G1_log2FoldChange < 1 ~ 'non'))




gene_exp_with_symbol_sort <- gene_exp_with_symbol %>% 
  arrange(desc(G2_vs_G1_log2FoldChange))

gene_list <- gene_exp_with_symbol_sort$G2_vs_G1_log2FoldChange 
names(gene_list) <- gene_exp_with_symbol_sort$Symbol
head(gene_list)

ego <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  #nPerm        = 1000, #置换检验的次数，默认为1000
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.2,
  verbose      = FALSE,
  seed = FALSE, 
  by = "fgsea")


KEGG_Flor_entrez <- gseKEGG(geneList = gene_list,
                            organism = 'mmu', 
                            pvalueCutoff = 0.25)


kegg_diff_Flor <- gene_exp_with_symbol_sort %>% 
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(Symbol, name=G2_vs_G1_log2FoldChange) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db) %>% 
  .$ENTREZID %>%
  enrichKEGG(gene = .,
             organism = 'mmu',  #例如，oas 代表绵羊，其它物种更改这行即可
             pAdjustMethod = 'fdr',  #指定 p 值校正方法
             pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
             qvalueCutoff = 0.2
           )

dotplot(
  kegg_diff_Flor,
  showCategory=10
  ) 

res <- gseGO(
  gene_list,    # 根据logFC排序的基因集
  ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Mm.eg.db,    # 使用人的OrgDb
  keyType = "ENSEMBL",    # 基因id类型
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",    # p值校正方法
)


res <- gene_exp_with_symbol_sort  %>%
  dplyr::filter(sig != 'non') %>%
  dplyr::pull(Symbol, name=G2_vs_G1_log2FoldChange) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Mm.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

kegg_diff_Flor <- gene_exp_with_symbol_sort %>% 
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(Symbol, name=G2_vs_G1_log2FoldChange) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Mm.eg.db) %>% 
  .$ENTREZID %>%
  enrichKEGG(gene = .,
             organism = 'mmu',  #例如，oas 代表绵羊，其它物种更改这行即可
             pAdjustMethod = 'fdr',  #指定 p 值校正方法
             pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
             qvalueCutoff = 0.2
  )



################  

dotplot(res,title="EnrichmentGO_MF_dot")
dotplot(res, showCategory=20,title="EnrichmentGO_MF")
barplot(kegg_diff_Flor, showCategory=20,title="EnrichmentGO_MF")











########### USP2 差异基因

USP2_diff_gene <- USP2_diff_gene %>% 
  mutate(sig=case_when(avg_log2FC <= -0.1 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.1 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.1 | avg_log2FC < 0.1 ~ 'non'))



#####  Arac, SPG 铁死亡，自噬, 
##############



cnetplot(kegg_diff_Flor)



###########   GSEA
gene_exp_with_symbol <- gene_exp_with_symbol %>% 
  mutate(sig=case_when(G2_vs_G1_log2FoldChange <=  1 & G2_vs_G1_FDR < 0.05 ~ 'down',
                       G2_vs_G1_log2FoldChange >=  1 & G2_vs_G1_FDR < 0.05 ~ 'up',
                       G2_vs_G1_FDR > 0.05 ~ 'non',
                       G2_vs_G1_log2FoldChange > -1 | G2_vs_G1_log2FoldChange < 1 ~ 'non'))



res_gsea <- gseGO(
  gene_list_ensemble,    # 根据logFC排序的基因集
  ont = "BP",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Mm.eg.db,    # 使用人的OrgDb
  keyType = "ENSEMBL",    # 基因id类型
  pvalueCutoff = 1,
  pAdjustMethod = "BH",    # p值校正方法
)

res_gsea_kegg <- gseKEGG(
  gene_list_ensemble,    # 根据logFC排序的基因集
  organism = "mmu",    # 人的拉丁名缩写
  pvalueCutoff = 1,
  pAdjustMethod = "BH"
)

gene_list_ensemble <- gene_exp_with_symbol_sort$G2_vs_G1_log2FoldChange 
names(gene_list_ensemble) <- gene_exp_with_symbol_sort$Gene_ID
head(gene_list_ensemble)


gseaplot2(res_gsea, title = 'autophagic cell death', geneSetID = 1)














