from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samap.utils import (save_samap, load_samap)
from samalg import SAM
import scanpy as sc
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Step 1: Loading SAMap Object
    print(f"[Step 1] Loading SAMap Object, started at {datetime.now().strftime('%H:%M:%S')}")
    samap = load_samap("results/samap_objects/samap_result.pkl")
    # Step 2: Store as .h5ad file
    adata = samap.samap.adata
    adata.write("results/samap.h5ad")

# Standard Python main
if __name__ == "__main__":
    main()