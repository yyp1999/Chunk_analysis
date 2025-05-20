ptm = Sys.time()
library(CellChat)
library(patchwork)
library(SingleCellExperiment)
setwd('//data//yypdata//NMF//Benchmark//cellchat//')
sce = readRDS('/data/yypdata/NMF/Benchmark/multinichenet/COVID.rds')
data <- as.matrix(counts(sce))  
meta <- as.data.frame(colData(sce))

meta_covid <- meta[meta$Status == "COVID", ]
meta_healthy <- meta[meta$Status == "Healthy", ]
data_covid <- data[, rownames(meta_covid)]
data_healthy <- data[, rownames(meta_healthy)]

cellchat_covid <- createCellChat(object = data_covid, 
                                 meta = meta_covid, 
                                 group.by = "cell_type_coarse")
cellchat_healthy <- createCellChat(object = data_healthy, 
                                   meta = meta_healthy, 
                                   group.by = "cell_type_coarse")
cellchat_covid@DB <- CellChatDB.human  
cellchat_covid <- subsetData(cellchat_covid)  
cellchat_covid <- identifyOverExpressedGenes(cellchat_covid)
cellchat_covid <- identifyOverExpressedInteractions(cellchat_covid)
cellchat_covid <- computeCommunProb(cellchat_covid, type = "triMean")  
cellchat_covid <- filterCommunication(cellchat_covid, min.cells = 10)  
cellchat_covid <- computeCommunProbPathway(cellchat_covid)
cellchat_covid <- aggregateNet(cellchat_covid)

cellchat_healthy@DB <- CellChatDB.human
cellchat_healthy <- subsetData(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedGenes(cellchat_healthy)
cellchat_healthy <- identifyOverExpressedInteractions(cellchat_healthy)
cellchat_healthy <- computeCommunProb(cellchat_healthy, type = "triMean")
cellchat_healthy <- filterCommunication(cellchat_healthy, min.cells = 10)
cellchat_healthy <- computeCommunProbPathway(cellchat_healthy)
cellchat_healthy <- aggregateNet(cellchat_healthy)

object.list <- list(Healthy = cellchat_healthy, COVID = cellchat_covid)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

save(object.list, file = "cellchat_object_list.RData")
save(cellchat, file = "cellchat_merged.RData")
load("cellchat_object_list.RData")
load("cellchat_merged.RData")

gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
gg1 + gg2  

par(mfrow = c(1, 2))
netVisual_diffInteraction(cellchat, weight.scale = TRUE) 
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")  

gg1 <- netVisual_heatmap(cellchat, measure = "count")  
gg2 <- netVisual_heatmap(cellchat, measure = "weight")  
gg1 + gg2

netVisual_bubble(cellchat, sources.use = c(0:14), targets.use =c(0:14),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, sources.use = c(0:14), targets.use = c(0:14),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in COVID", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(0:14), targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in COVID", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1

unique(gg1$data$interaction_name)
write.csv(unique(gg1$data$interaction_name), "cellchat_lri.csv", row.names = FALSE)
