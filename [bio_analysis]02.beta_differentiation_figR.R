#P. Augsornworawat
#Sample codes used for integration and analysis of multiome GEX and ATAC sequencing datasets
#this is a utilization of Signac package. Please refer to https://satijalab.org/signac/ for more details.
#Figure 1
#SC ISLETS
library(FigR)
library(SummarizedExperiment)
library(FNN)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(S4Vectors)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(ComplexHeatmap)
library(networkD3)
library(BiocParallel)
library(chromVAR)
library(JASPAR2020)
library(motifmatchr)
library(SeuratWrappers)
library(monocle3)
library(Matrix)

library(ggseqlogo)
library(chromVAR)
library(TFBSTools)
library(patchwork)

#scislet1, scislet2, week2 are download from GEO with https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199636

set.seed(1234)


# Set working directory
setwd("../pygot_data/beta_diff/multi_omics/")

#Laod gene annotations for human hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

#set processes for  Quantifying peaks in each dataset using FeatureMatrix
library(future)
plan("multicore", workers = 8)

options(future.globals.maxSize = 200000 * 1024^2) # for 80 Gb RAM

#READ RNA DATA
RNA.scislet1 <- Read10X_h5("./scislet1/filtered_feature_bc_matrix.h5")
RNA.scislet2 <- Read10X_h5("./scislet2/filtered_feature_bc_matrix.h5")
RNA.scislet3 <- Read10X_h5("./week2/filtered_feature_bc_matrix.h5")

#READ in peak sets
peaks.scislet1 <- read.table(file = "./scislet1/atac_peaks.bed",col.names = c("chr", "start", "end"))
peaks.scislet2 <- read.table(file = "./scislet2/atac_peaks.bed",col.names = c("chr", "start", "end"))
peaks.scislet3 <- read.table(file = "./week2/atac_peaks.bed",col.names = c("chr", "start", "end")) 

# convert to genomic ranges
gr.scislet1 <- makeGRangesFromDataFrame(peaks.scislet1)
gr.scislet2 <- makeGRangesFromDataFrame(peaks.scislet2)
gr.scislet3 <- makeGRangesFromDataFrame(peaks.scislet3)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.scislet1, gr.scislet2, gr.scislet3))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

#CREATE FRAGMENT FILES
# load metadata
md.scislet1 <- read.table(file = "./scislet1/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 
md.scislet2<- read.table(file = "./scislet2/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 
md.scislet3 <- read.table(file = "./week2/per_barcode_metrics.csv",stringsAsFactors = FALSE,sep = ",", header = TRUE,row.names = 1)[-1, ] 

# create fragment objects
frags.scislet1 <- CreateFragmentObject(path = "./scislet1/atac_fragments.tsv.gz",cells = rownames(md.scislet1))
frags.scislet2 <- CreateFragmentObject(path = "./scislet2/atac_fragments.tsv.gz",cells = rownames(md.scislet2))
frags.scislet3 <- CreateFragmentObject(path = "./week2/atac_fragments.tsv.gz",cells = rownames(md.scislet3))

#Quantify peaks in each dataset
scislet1.counts <- FeatureMatrix(fragments = frags.scislet1, features = combined.peaks,cells = rownames(md.scislet1))
scislet2.counts <- FeatureMatrix(fragments = frags.scislet2, features = combined.peaks,cells = rownames(md.scislet2))
scislet3.counts <- FeatureMatrix(fragments = frags.scislet3, features = combined.peaks,cells = rownames(md.scislet3))

#Create Seurat object using RNA Data
scislet1 <- CreateSeuratObject(counts = RNA.scislet1$`Gene Expression`, assay = "RNA")
scislet2 <- CreateSeuratObject(counts = RNA.scislet2$`Gene Expression`, assay = "RNA")
scislet3 <- CreateSeuratObject(counts = RNA.scislet3$`Gene Expression`, assay = "RNA")

#remove non-common cells in multiome data
include_list1 <- colnames(scislet1)
include_list2 <- colnames(scislet2)
include_list3 <- colnames(scislet3)

scislet1.counts<-scislet1.counts[, include_list1]
scislet2.counts<-scislet2.counts[, include_list2]
scislet3.counts<-scislet3.counts[, include_list3]

#Add RNA assay to Seurat Object
scislet1[["ATAC"]] <- CreateChromatinAssay(counts = scislet1.counts, sep = c(":", "-"), fragments = frags.scislet1, annotation = annotation)
scislet2[["ATAC"]] <- CreateChromatinAssay(counts = scislet2.counts, sep = c(":", "-"), fragments = frags.scislet2, annotation = annotation)
scislet3[["ATAC"]] <- CreateChromatinAssay(counts = scislet3.counts, sep = c(":", "-"), fragments = frags.scislet3, annotation = annotation)

#QUALITY CONTROL
DefaultAssay(scislet1) <- "ATAC"
DefaultAssay(scislet2) <- "ATAC"
DefaultAssay(scislet3) <- "ATAC"

scislet1 <- NucleosomeSignal(scislet1)
scislet2 <- NucleosomeSignal(scislet2)
scislet3 <- NucleosomeSignal(scislet3)

scislet1 <- TSSEnrichment(scislet1)
scislet2 <- TSSEnrichment(scislet2)
scislet3 <- TSSEnrichment(scislet3)

#visualize QC features
VlnPlot(object = scislet1,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)
VlnPlot(object = scislet2,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)
VlnPlot(object = scislet3,features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 4,pt.size = 0)

# filter out low quality cells
scislet1 <- subset( x = scislet1,subset = nCount_ATAC < 50000 & nCount_RNA < 50000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.25 & TSS.enrichment > 2)
scislet2 <- subset( x = scislet2,subset = nCount_ATAC < 40000 & nCount_RNA < 40000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.5 & TSS.enrichment > 2)
scislet3 <- subset( x = scislet3,subset = nCount_ATAC < 40000 & nCount_RNA < 40000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nucleosome_signal < 1.5 & TSS.enrichment > 2)



#Peak calling for each data
# call peaks using MACS2
peakscislet1 <- CallPeaks(scislet1 , macs2.path = "/disk/xuruihong/miniconda3/bin/macs3")
peakscislet2 <- CallPeaks(scislet2, macs2.path = "/disk/xuruihong/miniconda3/bin/macs3")
peakscislet3 <- CallPeaks(scislet3, macs2.path = "/disk/xuruihong/miniconda3/bin/macs3")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peakscislet1 <- keepStandardChromosomes(peakscislet1, pruning.mode = "coarse")
peakscislet1 <- subsetByOverlaps(x = peakscislet1, ranges = blacklist_hg38_unified, invert = TRUE)
peakscislet2 <- keepStandardChromosomes(peakscislet2, pruning.mode = "coarse")
peakscislet2 <- subsetByOverlaps(x = peakscislet2, ranges = blacklist_hg38_unified, invert = TRUE)
peakscislet3 <- keepStandardChromosomes(peakscislet3, pruning.mode = "coarse")
peakscislet3 <- subsetByOverlaps(x = peakscislet3, ranges = blacklist_hg38_unified, invert = TRUE)

# convert to genomic ranges
gr.scislet1 <- makeGRangesFromDataFrame(peakscislet1)
gr.scislet2 <- makeGRangesFromDataFrame(peakscislet2)
gr.scislet3 <- makeGRangesFromDataFrame(peakscislet3)

# Create a unified set of peaks to quantify in each dataset
recalled.combined.peaks <- reduce(x = c(gr.scislet1, gr.scislet2, gr.scislet3))

# Filter out bad peaks based on length
peakwidths <- width(recalled.combined.peaks)
recalled.combined.peaks <- recalled.combined.peaks[peakwidths  < 10000 & peakwidths > 20]



# quantify counts in each peak
macs2_countsscislet1 <- FeatureMatrix(fragments = scislet1$ATAC@fragments, features = recalled.combined.peaks, cells = colnames(scislet1))
macs2_countsscislet2 <- FeatureMatrix(fragments = scislet2$ATAC@fragments, features = recalled.combined.peaks, cells = colnames(scislet2))
macs2_countsscislet3 <- FeatureMatrix(fragments = scislet3$ATAC@fragments, features = recalled.combined.peaks, cells = colnames(scislet3))

# create a new assay using the MACS2 peak set and add it to the Seurat object
scislet1[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet1, fragments = frags.scislet1, annotation = annotation)
scislet2[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet2, fragments = frags.scislet2, annotation = annotation)
scislet3[["peaks"]] <- CreateChromatinAssay(counts = macs2_countsscislet3, fragments = frags.scislet3, annotation = annotation)


#MERGING of DATASETS
#save condition info
scislet1$dataset <- 'scislet1'
scislet2$dataset <- 'scislet2'
scislet3$dataset <- 'scislet3'

#rename cells with ID
scislet1 = RenameCells(scislet1, add.cell.id = "scislet1")
scislet2 = RenameCells(scislet2, add.cell.id = "scislet2")
scislet3 = RenameCells(scislet3, add.cell.id = "scislet3")

#RUN SCTransform
scislet1 <- SCTransform(scislet1,assay = "RNA")
scislet2 <- SCTransform(scislet2,assay = "RNA")
scislet3 <- SCTransform(scislet3,assay = "RNA")

#Fomd variable feature and Run PCA
scislet1 <- FindVariableFeatures(scislet1)
scislet2 <- FindVariableFeatures(scislet2)
scislet3 <- FindVariableFeatures(scislet3)

scislet1 <- RunPCA(scislet1)
scislet2 <- RunPCA(scislet2)
scislet3 <- RunPCA(scislet3)

#ANCHORING datasets and merge objects using integrate
# find integration anchors 
integratefeatures<- SelectIntegrationFeatures(object.list = list(scislet1, scislet2, scislet3),nfeatures = 500)
anchors <- PrepSCTIntegration(object.list = list(scislet1, scislet2, scislet3))

anchors <- FindIntegrationAnchors(object.list = anchors,reduction = "rpca",dims = 1:20,normalization.method = "SCT", anchor.features = integratefeatures,  assay = c('SCT', 'SCT','SCT'))

#INTEGRATE data and create a new merged object using SCT data
integratedRNA <- IntegrateData(anchorset = anchors, normalization.method = "SCT",dims = 1:20,new.assay.name = "integratedRNA")

saveRDS(integratedRNA, './integratedRNA.rds')
saveRDS(anchors, './anchors.rds')

rm(macs2_countsscislet1)
rm(macs2_countsscislet2)
rm(macs2_countsscislet3)

rm(scislet1)
rm(scislet2)
rm(scislet3)

rm(scislet1.counts)
rm(scislet2.counts)
rm(scislet3.counts)

rm(RNA.scislet1)
rm(RNA.scislet2)
rm(RNA.scislet3)




#INTEGRATE data and create a new merged object using LSI 
integratedlsi<-integratedRNA
DefaultAssay(integratedlsi) <- "peaks"
integratedlsi <- FindTopFeatures(integratedlsi, min.cutoff = 5)
integratedlsi <- RunTFIDF(integratedlsi)
integratedlsi <- RunSVD(integratedlsi)
integratedlsi <- IntegrateEmbeddings(anchorset = anchors, new.reduction.name = "integratedLSI",reductions = integratedlsi@reductions[["lsi"]])

#Combined assays and reduction into one seurat object
integrated <- integratedlsi
integrated@assays[["integratedRNA"]] <- integratedRNA@assays[["integratedRNA"]]

DefaultAssay(integrated) <- "RNA"
integrated[["RNA"]] <- as(integrated[["RNA"]], "Assay")
sceasy::convertFormat(integrated, from="seurat", to="anndata", assay='RNA', main_layer ='counts',
                      outFile='./integratedRNA.h5ad', )


#JOINT UMAP VISUALIZATION
DefaultAssay(integrated) <- "integratedRNA"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = FALSE)

DefaultAssay(integrated) <- "peaks"
integrated <- FindTopFeatures(integrated, min.cutoff = 5)
integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)
# build a joint neighbor graph using both assays
integrated <- FindMultiModalNeighbors(object = integrated,reduction.list = list("pca", "integratedLSI"), dims.list = list(1:24, 2:24),modality.weight.name = "RNA.weight", verbose = TRUE)
# build a joint UMAP visualization
integrated <- RunUMAP(object = integrated, nn.name = "weighted.nn",assay = "SCT",verbose = TRUE)
integrated <- FindClusters(integrated,resolution = 0.08, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#PLOT UMAP figures
DefaultAssay(integrated) <- "SCT"
DimPlot(integrated, label = TRUE, repel = TRUE, reduction = "umap") 
DimPlot(integrated, label = TRUE, repel = TRUE, reduction = "umap", group.by = 'dataset') 


#Find gene expression markers for each cluster
DefaultAssay(integrated) <- "SCT"
integrated <- PrepSCTFindMarkers(integrated)
cluster0 <- FindMarkers(integrated, ident.1 = 0,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster1 <- FindMarkers(integrated, ident.1 = 1,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster2 <- FindMarkers(integrated, ident.1 = 2,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster3 <- FindMarkers(integrated, ident.1 = 3,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster4 <- FindMarkers(integrated, ident.1 = 4,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster5 <- FindMarkers(integrated, ident.1 = 5,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster6 <- FindMarkers(integrated, ident.1 = 6,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 
cluster7 <- FindMarkers(integrated, ident.1 = 7,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE)
cluster8 <- FindMarkers(integrated, ident.1 = 8,ident.2 = NULL, logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE) 


DotPlot(integrated, features = c('CHGA', 'PDX1', 'NEUROD1', 'TPH1', 'DDC', 'INS', 'ERO1B', 'GCG', 'LOXL4',
                                 'SST','HHEX', 'PURPL','PTCHD4', 'VIM', 'SPARC', 'MKI67', 'PCNA', 'KRT19', 'KRT7')) + RotatedAxis()

DotPlot(integrated, features = c('CHGA', 'PDX1', 'NEUROD1', 'TPH1', 'DDC', 'INS', 'ERO1B', 'GCG', 'LOXL4',
                                 'SST','HHEX', 'PURPL','PTCHD4', 'VIM', 'SPARC', 'MKI67', 'PCNA', 'KRT19', 'KRT7')) + RotatedAxis()

#RENAME clusters with identity based on expression markers
integrated <- RenameIdents(integrated, '0'='EC','1'='Beta','2'='Alpha','3'='EC', '4'='Delta','5'='EC2','6'='Mesenchyme', '7'='Exocrine', '8'='Proliferating alpha cell')
integrated$celltype <- integrated@active.ident

#ANNOTATED UMAP COMBINED PLOT WITH COLORS
DefaultAssay(integrated) <- "SCT"
DimPlot(integrated, label = FALSE, repel = TRUE,group.by="celltype", reduction = "umap", cols= c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327")) 
DimPlot(integrated, label = FALSE, repel = TRUE, reduction = "umap", group.by = 'dataset', cols = c("#2f4b7c", "#d45087", "#ffa600"), shuffle=TRUE) 

#ATAC PEAKS plots of key hormone genes
DefaultAssay(integrated) <- "peaks"
Idents(integrated) <- "celltype"
CoveragePlot(object = integrated,  region = "INS",  features = "INS", expression.assay = "SCT", extend.upstream = 5000,  extend.downstream = 2000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))
CoveragePlot(object = integrated,  region = "GCG",  features = "GCG", expression.assay = "SCT", extend.upstream = 1000,  extend.downstream = 7000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))
CoveragePlot(object = integrated,  region = "SST",  features = "SST", expression.assay = "SCT", extend.upstream = 5000,  extend.downstream = 2000) & scale_fill_manual(values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))


#Add the gene activity (promoter accessibility) matrix to the Seurat object as a new assay and normalize it
integrated[['ActivityRNA']] <- CreateAssayObject(counts = gene.activities)
gene.activities <- GeneActivity(integrated) 
integrated <- NormalizeData(object = integrated, assay = 'ActivityRNA', normalization.method = 'LogNormalize', scale.factor = median(integrated$nCount_ActivityRNA))


#MOTIF activity analysis
#Adding motif information to the Seurat object using chromVAR package
DefaultAssay(integrated) <- 'peaks'
pfm <- getMatrixSet(
  x = JASPAR2020,  # 
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

integrated <- AddMotifs(
  object = integrated,  
  genome = BSgenome.Hsapiens.UCSC.hg38,  
  pfm = pfm  
)

integrated <- RunChromVAR(object = integrated, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)


#Use this to find motif id from gene name
DefaultAssay(integrated) <- 'peaks'
motif.name <- ConvertMotifID(integrated, name = c('TPI1', 'PDX1', 'PAX6', 'FEV', 'MNX1'))
motif.name

#integrated <- RenameIdents(integrated, 'EC2'='EC')
DotPlot(integrated, motif.name, assay = 'chromvar', )

p1 <- DimPlot(integrated)
p2 <- FeaturePlot(integrated, features = motif.name, min.cutoff = 'q50',
            max.cutoff = 'q90',
            pt.size = 0.1)
p1 + p2

DotPlot(integrated, c('ERO1B', 'SLC30A8', 'ACVR1C', 'INS', 'ISL1'), assay = 'SCT')
CoveragePlot(object = integrated,  
             region = "PTF1A",  
             features = "PTF1A",
             extend.upstream = 15000,)

CoveragePlot(object = integrated,  
             region = "PTF1A",  
             features = "PTF1A", 
             expression.assay = "SCT", 
             extend.upstream = 5000,  
             extend.downstream = 2000) & scale_fill_manual(
               values = c("#FFA500", "#8E2C85", "#134D9C","#E0082D","#D877CF","#7D000E","#AAAAAA","#679327"))


integrated <- LinkPeaks(
  integrated,
  peak.assay = 'peaks',
  expression.assay = 'SCT',
  #genes.use = c("TENM2")
)
saveRDS(integrated, './integrated.rds')

integrated <- readRDS('./integrated.rds')

integrated_beta_ec <-  integrated[,Idents(integrated) %in% c('EC', 'Beta')]
integrated_beta_ec <- LinkPeaks(
  integrated_beta_ec,
  peak.assay = 'peaks',
  expression.assay = 'SCT',
  #genes.use = c("TENM2")
)
#Use this to find motif id from gene name
DefaultAssay(integrated_beta_ec) <- 'peaks'

integrated_beta_ec$RNA

#FigR



atac.se <- SummarizedExperiment(
  assays = c(counts=integrated_beta_ec$peaks$counts, data=integrated_beta_ec$peaks@data),
  rowRanges = integrated_beta_ec$peaks@ranges,
  colData = integrated_beta_ec@meta.data
)
atac.se

DefaultAssay(integrated_beta_ec) <- 'RNA'
integrated_beta_ec <- NormalizeData(integrated_beta_ec, normalization.method = "LogNormalize", scale.factor = 10000)

# Run using multiple cores if parallel support
cisCorr <- runGenePeakcorr(atac.se, integrated_beta_ec$RNA$data, 
                          genome = "hg38",
                          nCores = 10, 
                          p.cut=NULL)



cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)

dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 5, # No. sig peaks needed to be called a DORC
                       labelTop = 20,
                       returnGeneList = TRUE, # Set this to FALSE for just the plot
                       force=2)
dorcGenes
# Unfiltered
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

dorcMat <- getDORCScores(ATAC.se = atac.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 10)



cellkNN <- get.knn(integrated_beta_ec@reductions$integratedLSI@cell.embeddings,k = 30)$nn.index
rownames(cellkNN) <- colnames(integrated_beta_ec)
# Smooth dorc scores using cell KNNs (k=30)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 4)

RNAmat <- integrated_beta_ec$RNA$data
RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]
# Smooth RNA using cell KNNs
# This takes longer since it's all genes
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat, nCores = 10)


# Visualize on pre-computed UMAP
umap.d <- as.data.frame(integrated_beta_ec@reductions$umap@cell.embeddings)

gene2plot = 'INS'
# DORC score for Dlx3
p <- DimPlot(integrated_beta_ec)
dorcg <- plotMarker2D(umap.d,dorcMat.s,markers = gene2plot,maxCutoff = "q0.99",colorPalette = "brewer_heat") + ggtitle("DORC")
rnag <- plotMarker2D(umap.d,RNAmat.s,markers = gene2plot,maxCutoff = "q0.99",colorPalette = "brewer_purple") + ggtitle("RNA")
p + dorcg + rnag



figR.d <- runFigRGRN(ATAC.se = atac.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = 4)


figR.d %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

rankDrivers(figR.d,rankBy = "meanScore",interactive = FALSE)

plotDrivers(figR.d,score.cut = 2,marker = "Lef1")

write.table(figR.d, "./FigR_GRN.txt", sep='\t')

figR.d <- read.table("./FigR_GRN.txt", sep='\t')

plotfigRHeatmap(figR.d = figR.d,
                score.cut = 0.5,
                TFs = c("ISL1", "PDX1", "PAX6","MNX1", "FEV", "NEUROG3", "XBP1"),
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
)
DefaultAssay(integrated) <- 'SCT'
FeaturePlot(integrated, c('ASPH', 'ERO1B', 'GLS', 'IRX2', 'ISL1', 'NEFM', 'NPTX2', 'TMOD1'))

plotfigRNetwork(figR.d,
                score.cut = 1,
                TFs = c("ISL1", "PDX1", "PAX6","MNX1", "FEV", "NEUROG3", "XBP1", ),
                weight.edges = TRUE)

