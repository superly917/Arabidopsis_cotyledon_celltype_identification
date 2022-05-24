library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
seurat_ob = readRDS("example.rds")

# the interested marker
topn_markers2vis=c("RBCS-1A","LHCB1.1","HDG2","POLAR","SPCH","MUTE",
                    "EPF2","BASL","EPF1","HIC","FAMA","SCRM",
                    "TMM","IQD5")
# Palette define
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))
CustomCol2 <- function(n){
  my_palette=c(
    "#7fc97f","#beaed4","#fdc086","#386cb0","#f0027f","#a34e3b",
    "#666666","#1b9e77","#d95f02","#7570b3","#d01b2a","#43acde",
    "#efbd25","#492d73","#cd4275","#2f8400","#d9d73b","#aed4ff",
    "#ecb9e5","#813139","#743fd2","#434b7e","#e6908e","#214a00",
    "#ef6100","#7d9974","#e63c66","#cf48c7","#ffe40a","#a76e93",
    "#d9874a","#adc64a","#5466df","#d544a1","#54d665","#5e99c7",
    "#006874","#d2ad2c","#b5d7a5","#9e8442","#4e1737","#e482a7",
    "#6f451d","#2ccfe4","#ae6174","#a666be","#a32b2b","#ffff99",
    "#3fdacb","#bf5b17")
  return(my_palette[n])
}

# featureplot
ggfeature = lapply(topn_markers2vis, function(x) FeaturePlot(seurat_ob,
                    features = x,
                    reduction= "tsne",
                    ncol = 2, pt.size = 0.2) +
                    theme( plot.title = element_text(hjust = 0.5)) + sc)

nrow = ceiling(length(topn_markers2vis)/2)
pdf("topmarker_gene_featureplot.pdf",width = 10, height = length(topn_markers2vis)*2)
grid.arrange(grobs = ggfeature, nrow = nrow, ncol=2)
dev.off()
png("topmarker_gene_featureplot.png",width = 12, height = nrow*6, res = 96, units = "in")
grid.arrange(grobs = ggfeature, nrow = nrow, ncol=2)
dev.off()

# vlnplot
groupby = "clusters"
colors2use = CustomCol2(1:length(unique(seurat_ob@meta.data[,groupby])))
ggvln = lapply(topn_markers2vis, function(x) VlnPlot(seurat_ob, features = x, cols = colors2use,
            pt.size = 0,
            group.by = groupby) +
            labs(title = "",y = x) + 
            theme(legend.position = "none", 
            axis.title.x = element_text(size = 0),
            axis.title.y = element_text(size = 8),
            axis.text.y=element_text(size = 8)))
pdf("topmarker_gene_violin_plot.pdf",width = 12, height = length(topn_markers2vis)*2)
grid.arrange(grobs = ggvln, ncol=1)
dev.off()
png("topmarker_gene_violin_plot.png",width = 12, height = nrow*6, res = 96, units = "in")
grid.arrange(grobs = ggvln, ncol=1)
dev.off()


# Identify celltype based on the SingleR reference dataset
library(SingleR)
library(scater)
library(dplyr)

seurat_ob = readRDS("example.rds")
query.m = GetAssayData(seurat_ob,assay = "RNA",slot = "counts")
query.sce = SingleCellExperiment(assays = list(counts = query.m))
query.sce = logNormCounts(query.sce)
# 
ref.sce = readRDS("reference.rds")

pred = SingleR(query.sce,ref.sce,labels=factor(ref.sce$celltype),BPPARAM = MulticoreParam(workers = 10))
saveRDS(pred,"singleR.rds")

#add the celltyping result as the metadata for each cell
seurat_ob <- AddMetaData(object = seurat_ob,
                         metadata = pred$labels, col.name = "raw_celltype")
seurat_ob@meta.data$clusters = factor(seurat_ob@meta.data$clusters,levels=sort(as.numeric(unique(seurat_ob@meta.data$clusters))))
main_celltyping_df = seurat_ob@meta.data %>% mutate(cell_barcode = row.names(seurat_ob@meta.data)) %>%
                     select(cell_barcode, raw_celltype, clusters, sampleid)
write.table(main_celltyping_df, quote = F,
            "celltyping_results.xls", sep ="\t",row.names =F)
#count the cell number for each cell type
main_celltyping_stat = main_celltyping_df %>%
                       select( clusters, raw_celltype)%>%
                       group_by(clusters)%>%
                       count(raw_celltype) %>%
                       rename(cell_num=n)
write.table(main_celltyping_stat, quote = F, "ref_celltyping_statistics.xls", sep ="\t",row.names =F)
ggheat = plotScoreHeatmap(pred, clusters=seurat_ob@meta.data$clusters, max.labels = 10 ) 
ggsave("celltyping_heatmap.pdf",plot=ggheat)
ggsave("celltyping_heatmap.png",plot=ggheat)

seurat_ob = SetIdent(seurat_ob, value = "clusters")
top_celltype = main_celltyping_stat %>% group_by(clusters) %>% top_n(1,cell_num)
write.table(top_celltype, quote = F, "top_celltyping_statistics.xls", sep ="\t",row.names =F)

from.id = as.vector(top_celltype$clusters)
to.id = as.vector(top_celltype$raw_celltype)
seurat_ob = SetIdent(seurat_ob, value = plyr::mapvalues(x=Idents(seurat_ob),from=from.id,to=to.id))
seurat_ob = StashIdent( seurat_ob, save.name = "celltype")
color_level = length(levels(seurat_ob@meta.data$celltype))
ggtsne2 = DimPlot(object = seurat_ob, reduction = "tsne",pt.size = 1)+ theme( plot.title = element_text(hjust = 0.5)) + scale_colour_manual( values = CustomCol2(1:color_level))
ggsave("celltyping.pdf",plot= ggtsne2)
ggsave("celltyping.png",plot= ggtsne2)

