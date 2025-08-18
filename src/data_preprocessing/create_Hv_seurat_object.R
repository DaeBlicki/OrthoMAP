#' ----------------------------------------------------------------------------
#' @title:      Create Hv.Robj
#'
#' @description:
#' It loads `aepAtlasDubInclude.rds` [https://pubmed.ncbi.nlm.nih.gov/36639202/]
#' and identifys the doublets using `curatedIdent` column in meta.data.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 14/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
library(readr)
library(Seurat)

# load .Robj included with doublets and nondubs list
data <- readRDS(file.path("download", "HVAEP.Robj"))
nondubs <- read_tsv(file.path("download", "HVAEP_nondubs.tsv"),
                    col_names = FALSE)

# remove the doublets as described from source code of the project
nondubs <- unlist(nondubs)
data <- subset(data, cells = nondubs)
nondubs <- rownames(
  data@meta.data[!(data@meta.data$seurat_clusters %in% c(37, 41)), ]
)
data <- subset(data, cells = nondubs)
data <- UpdateSeuratObject(data)
# -----------------------------------------------------------------------------
# create IDs column: Neurons, cnidocytes, glands,
#                    GC/SC, Endoderm, Ectoderm, Dubs
# -----------------------------------------------------------------------------
neuro <- c("I_Neuro", "I_Ec1N", "I_Ec2N", "I_Ec3N", "I_Ec4N", "I_Ec5N",
           "I_En1N", "I_En2N", "I_En3N")
nem <- c("I_EarlyNem",
         "I_StenoNB", "I_DesmoNB", "I_IsoNB",
         "I_StenoNC", "I_DesmoNC", "I_IsoNC", "I_maturingNC")
gl <- c("I_ZymogenGl", "I_GlProgen", "I_GranularGl", "I_SpumousGl")
sc_gc <- c("I_FemGc", "I_MaleGC", "Ec_SomaticGerm",
           "En_BodyCol/SC", "Ec_BodyCol/SC", "I_ISC")
dubs <- c("Ec_NB_Dubs", "En_Gl_Dubs", "Ec_NC_Dubs", "En_NC_Dubs")
endo <- c("En_Head", "En_Foot", "En_Tent", "En_Hypo")
ecto <- c("Ec_Head", "Ec_Peduncle", "Ec_Battery", "Ec_BasalDisk")
celltype_list <- list(neuro, nem, gl, sc_gc, endo, ecto, dubs)
tissue_list <- list("neurons", "cnidocytes", "glands", "GC/SC",
                    "endoderm", "ectoderm", "dubs")
# Create Hashmap (Named Vectors)
# -------------------------------------------------------------------------
#' @title create_celltype_to_tissue_map
#' @brief maps celltype to tissue in Hydra v.
#' @param celltype_vector Character Vector, cell type in CuratedIdent
#' @param tissue String, cell type mapped to tissue in IDs
#' @return Named Vector, mapping cell type to tissue
# -------------------------------------------------------------------------
create_celltype_to_tissue_map <- function(celltype_vector, tissue) {
  tissue_vector <- rep(tissue, length(celltype_vector))
  setNames(tissue_vector, celltype_vector)
}
celltype_to_tissue_map <- unlist(mapply(create_celltype_to_tissue_map,
                                        celltype_list, tissue_list))
celltype <- as.character(sapply(data$curatedIdent, function(x) {
  celltype_to_tissue_map[as.character(x)]
}))

data$IDs <- celltype

# update and save seurat object
save(data, file = file.path("data", "single_cell_atlas", "Hv.Robj"))

# Clean workspace
rm(list = ls())