from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samap.utils import (save_samap, load_samap)
from samalg import SAM
from datetime import datetime
import scanpy as sc
import pandas as pd
import numpy as np

def main():
    # Step 1: Loading SAMap Object and get core parameters
    samap = load_samap("results/samap_objects/samap_result.pkl")
    adata = samap.samap.adata
    cell_ids = adata.obs.index

    # Step 2: Loop over resolution and random seed
    # Get the best clustering based on ARI for each random seed
    resolutions = np.arange(0.5, 5.0 + 0.5, 0.5)
    random_states = np.arange(1, 20 + 1, 1)
    for res in resolutions:
        print(f"Current Resolution: {res}, started at {datetime.now().strftime('%H:%M:%S')}")
        # create data frame for resolution
        df = pd.DataFrame({"Cell_ID": cell_ids})
        for random_state in random_states:
            print(f"Current Random State: {random_state}")
            sc.tl.leiden(adata, resolution=res, random_state=random_state, key_added='integrated_leiden')
            df[f"Iteration_{random_state}"] = adata.obs.loc[df["Cell_ID"], 'integrated_leiden'].values
            del adata.obs['integrated_leiden']
        # Save csv for each resolution
        filename = f"results/leiden_clusters/resolution_{res}.csv"
        df.to_csv(filename, index=False)

# Standard Python main
if __name__ == "__main__":
    main()