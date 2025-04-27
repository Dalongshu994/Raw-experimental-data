library(CellChat)
cellchat = createCellChat(object=mydata, group.by="cell_type", meta=mydata@meta.data)
groupSize = as.numeric(table(cellchat@idents))
CellChatDB = CellChatDB.human
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search="Cell-Cell Contact")
cellchat@DB <- CellChatDB.use
#####
cellchat = subsetData(cellchat)
future::plan("multisession", workers=4)

cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat) 

cellchat = projectData(cellchat, PPI.human)

#####################################################################################

cellchat = computeCommunProb(cellchat, raw.use=FALSE, population.size=TRUE)

cellchat = filterCommunication(cellchat, min.cells=10)

cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
#####################################################################################
####### 
b=netVisual_bubble(cellchat, sources.use=c(1,3,4,5,6,7,8), targets.use=c(2), color.heatmap="Spectral")+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=15, hjust=0.5, face="bold", size=10), axis.text.y=element_text(face="bold", size=10))
ggsave("cell_chat.pdf", b, width=12, height=15)