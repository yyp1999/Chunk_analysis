library(Seurat)
library(scriabin)
library(tidyverse)
library(pbapply)
library(future)
library(dplyr)
library(tidyr)
library(readr)

plan(sequential)
future::plan(sequential)

seurat_obj <- readRDS("seurat_COVID_processed.rds")

dim(LayerData(seurat_obj, assay = "RNA", layer = "data")) 
dim(seurat_obj[["pca"]]@cell.embeddings)  
head(seurat_obj@meta.data)  
table(seurat_obj$Status)  
table(seurat_obj$cell_type_coarse)  
table(seurat_obj$cell_type_fine)  

#interaction-programs
DimPlot(seurat_obj , label = T, repel = T) + NoLegend()
DimPlot(seurat_obj, label = T, repel = T, group.by = "cell_type_coarse") + NoLegend()
DimPlot(seurat_obj, group.by = "Status")
seurat_obj <- SeuratWrappers::RunALRA(seurat_obj)

seurat_obj_ip <- FindAllInteractionPrograms(seurat_obj, group.by = "Status", 
                                      cell_types = "cell_type_coarse", assay = "alra")
seurat_obj_ip_sig <- InteractionProgramSignificance(seurat_obj_ip)
seurat_obj <- ScoreInteractionPrograms(seurat_obj, mods = seurat_obj_ip_sig)

covid_ip_ligands <- FindMarkers(seurat_obj, group.by = "Status", ident.1 = "COVID", assay = "IPligands", test.use = "wilcox_limma")
covid_ip_receptors <- FindMarkers(seurat_obj, group.by = "Status", ident.1 = "COVID", assay = "IPreceptors", test.use = "wilcox_limma")

poi <- covid_ip_ligands %>% rownames_to_column("IP") %>% 
  top_n(n = 50, wt = avg_log2FC) %>% pull(IP)
poi2 <- covid_ip_receptors %>% rownames_to_column("IP") %>% 
  top_n(n = 50, wt = avg_log2FC) %>% pull(IP)

features <- seurat_obj_ip_sig %>% dplyr::filter(name==poi) %>%
  separate(lr_pair, into = c("ligand","receptor"), sep = "=") %>% pull(ligand) %>% unique()
features2 <- seurat_obj_ip_sig %>% dplyr::filter(name==poi2) %>%
  separate(lr_pair, into = c("ligand","receptor"), sep = "=") %>% pull(receptor) %>% unique()

LRI <- seurat_obj_ip_sig %>% dplyr::filter(name==poi) %>% unique()
LRI2 <- seurat_obj_ip_sig %>% dplyr::filter(name==poi2) %>% unique()

DotPlot(seurat_obj, features = features, group.by = "Status")
DotPlot(seurat_obj, features = features, group.by = "cell_type_coarse")

dotplot <- DotPlot(seurat_obj, features = features, group.by = "cell_type_coarse")
dotplot <- dotplot + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
print(dotplot)

dotplot2 <- DotPlot(seurat_obj, features = features2, group.by = "cell_type_coarse")
dotplot2 <- dotplot2 + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
print(dotplot2)

LRI_combined <- bind_rows(LRI, LRI2) %>%
  distinct(lr_pair, .keep_all = TRUE) %>%
  separate(lr_pair, into = c("ligand", "receptor"), sep = "=",remove = FALSE)


write_csv(LRI_combined %>% select(lr_pair, ligand, receptor),
          "combined_LR_pairs.csv")


