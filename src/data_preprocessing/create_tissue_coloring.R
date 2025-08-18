#' ----------------------------------------------------------------------------
#' @title:      Create tissue palette for Seurat
#'
#' @description:
#' It generates a color palette in Seurat with similar coloring for each
#' tissue group.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 28/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
tissue_order <- c(
  # Neurons
  "Ac.neural", "Nv.neuronal", "Nv.neuroglandular", "Hv.neurons",
  # Cnidocytes
  "Ac.cnidocyte", "Nv.cnidocyte", "Nv.mat.cnido", "Hv.cnidocytes",
  # Glands
  "Ac.gland.dig", "Ac.gland.muc", "Nv.gland.mucous", "Hv.glands",
  # GC/SC
  "Hv.GC/SC", "Nv.pSC", "Nv.PGCs",
  # Ectoderm and Epidermis
  "Ac.epidermis", "Nv.epithelia.ect", "Hv.ectoderm",
  "Nv.pharyngeal.ect", "Nv.ectoderm.embryonic",
  # Endoderm and Gastrodermis
  "Ac.gastrodermis", "Nv.gastrodermis", "Nv.mesendoderm.embryonic",
  "Hv.endoderm",
  # Muscle
  "Ac.muscle.st", "Nv.retractor muscle",
  # Others
  "Ac.unchar.1", "Nv.unchar.immune", "Hv.dubs"
)

tissue_colors <- c(
  # 💙 Neurons – cool, electric tones
  "#00429d", "#2a7fff", "#6baed6", "#9ecae1",
  # 💛 Cnidocytes – sharp, energetic yellows & golds
  "#f7dc6f", "#f1c40f", "#ffd700", "#e6b800",
  # 💜 Glands – vibrant purples & violets
  "#6a0dad", "#8e44ad", "#a569bd", "#c39bd3",
  # 💖 GC/SC – bold magentas and corals
  "#d33682", "#ff00ff", "#ff77ff",
  # 🧡 Epidermis – warm orange & copper tones
  "#e67e22", "#ff8c00", "#ffa500", "#f5b041", "#f8c471",
  # 💚 Endoderm – nature-inspired greens
  "#1a9850", "#66c2a5", "#2ecc71", "#a1d99b",
  # ❤️ Muscle – striking reds
  "#c0392b", "#e74c3c",
  # 🩶 Other / Uncharacterized – neutral greys
  "#4d4d4d", "#999999", "#d9d9d9"
)

palette_df <- data.frame(
  Tissue = tissue_order,
  Color = tissue_colors,
  stringsAsFactors = FALSE
)

# Save as CSV
filename <- file.path("data", "tissue_palette.csv")
write.csv(palette_df, filename, row.names = FALSE)

# Clean workspace
rm(list = ls())
