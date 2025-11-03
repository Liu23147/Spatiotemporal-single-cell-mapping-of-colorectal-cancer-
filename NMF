#######NMF#######
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
scRNA1 <- runNMF(scRNA1, k = ndim, assay="SCT")
scRNA2.list <- SplitObject(scRNA1,
                           split.by ="tissue")

geneNMF.programs <- multiNMF(scRNA2.list, assay="SCT", slot="data", k=4:9, nfeatures = 1000)
Idents(scRNA1) <- 'tissue'
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=9,
                                        weight.explained = 0.7,
                                        max.genes=100)

ph <- plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)
library(msigdbr)
library(fgsea)library(clusterProfiler)
library(org.Hs.eg.db)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(scRNA1), category = "H")
})
head(top_p$MP6)

MP3fuji=top_p$MP3
dotplotGsea(data =top_p$MP3,topn = 20,
            order.by = 'NES',
            add.seg = T)

ggplot(top_p$MP9, aes(x = -log10(pval), y = pathway)) +
  geom_segment(aes(xend = 0, yend = pathway), color = "gray") +
  geom_point(aes(size = overlap), color = "#2c7bb6") +
  scale_size_continuous(range = c(4, 10)) +
  labs(
    x = expression(-log[10](pval)),
    y = NULL,
    size = "Overlap Genes",
    title = "Pathway Enrichment-MP9"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank())


mp.genes <- geneNMF.metaprograms$metaprograms.genes
scRNA1 <-AddModuleScore(scRNA1, features = mp.genes, assay="SCT", ncores=20, name = "MP")
DotPlot(scRNA1, features=names(mp.genes))


jjDotPlot(object = scRNA1,
          gene =names(mp.genes),
          xtree = F,
          rescale = T,
          point.geom = F,
          tile.geom = T,id='tissue',assay='SCT',ytree =F)
