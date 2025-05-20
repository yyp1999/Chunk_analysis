library(GSVA)
library(GSEABase) 
library(survival)
library(survminer)

#OSSC
expr_data <- read.csv("F:\\pythonproject\\Tensor\\Chunk_v2\\data\\OSSC\\GSVA/lr_expr_df.csv", row.names = 1, header = TRUE, check.names = FALSE)
lr_gene <- read.csv("F:\\pythonproject\\Tensor\\Chunk_v2\\data\\OSSC\\GSVA/lr_genes_res2.csv", row.names = 1, header = TRUE, check.names = FALSE)
expr_data <- as.matrix(expr_data) 

expr_genes <- rownames(expr_data)
#gene_list <- rownames(lr_gene)  
gene_list <- lr_gene$gene
#gene_set <- list("OSSC_relevant_LR_gene" = gene_list)  

filtered_gene_list <- gene_list %>% intersect(expr_genes)

gene_set <- list("OSSC_relevant_LR_gene" = filtered_gene_list)

param <- gsvaParam(
  exprData = expr_data,         
  geneSets = gene_set,         
  minSize = 5,                 
  maxSize = 1000,               
  kcdf = "Gaussian",           
)

set.seed(123)
gsva_results <- gsva(param)
gsva_df <- as.data.frame(gsva_results)
write.csv(gsva_df, "F:\\pythonproject\\Tensor\\Chunk_v2\\result\\benchmark\\GSVA\\GSVA_res\\bench_GSVA_Chunk_logtpm.csv")

