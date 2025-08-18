library(mclust)

filepath <- "results/OrthoMAP_Seurat_Objects/mmseqs2/OOG/result2.Robj"
orthomap_obj <- get(load(filepath))

mclust::adjustedRandIndex(orthomap_obj[, "IDs"], orthomap_obj$orhtomap_clusters)
