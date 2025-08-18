from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samap.utils import save_samap
from samalg import SAM
import scanpy as sc
from datetime import datetime
import pandas as pd

def main():
    # Step 1: Loading in raw data
    print(f"[Step 1] Loading in raw data, started at {datetime.now().strftime('%H:%M:%S')}")
    fn_ac = "data/single_cell_atlas/Ac.h5ad"
    fn_hv = "data/single_cell_atlas/Hv.h5ad"
    fn_nv = "data/single_cell_atlas/Nv.h5ad"

    # Step 2: Loading in gene map
    print(f"[Step 2] Loading gene maps, started at {datetime.now().strftime('%H:%M:%S')}")
    annot_ac = pd.read_csv("data/annotation_table/Ac.csv")
    gene_ac = list(zip(annot_ac["geneID"], annot_ac["gene.name"]))
    print(gene_ac[:5])
    annot_hv = pd.read_csv("data/annotation_table/Hv.csv")
    gene_hv = list(zip(annot_hv["geneID"], annot_hv["gene.name"]))
    print(gene_hv[:5])
    annot_nv = pd.read_csv("data/annotation_table/Nv.csv")
    gene_nv = list(zip(annot_nv["geneID"], annot_nv["gene.name"]))
    print(gene_nv[:5])
    gene_names = {'Ac':gene_ac, 'Hv':gene_hv, 'Nv':gene_nv}

    # Step 3: Loading in SAM objects directly
    print(f"[Step 3.1] SAM on Aurelia coerulea : started at {datetime.now().strftime('%H:%M:%S')}")
    adata_ac = sc.read_h5ad(fn_ac)
    adata_ac.obs["orig.ident"] = adata_ac.obs["orig.ident"].astype(str)
    sam_ac=SAM(counts = adata_ac)
    sam_ac.run(batch_key = "orig.ident")
    sam_ac.save_anndata("results/sam_objects/Ac.h5ad")

    print(f"[Step 3.2] SAM on Hydra vulgaris: started at {datetime.now().strftime('%H:%M:%S')}")
    adata_hv = sc.read_h5ad(fn_hv)
    adata_hv.obs["orig.ident"] = adata_hv.obs["orig.ident"].astype(str)
    sam_hv=SAM(counts = adata_hv)
    sam_hv.run(batch_key = "orig.ident")
    sam_hv.save_anndata("results/sam_objects/Hv.h5ad")
    
    print(f"[Step 3.3] SAM on Nematostella vectensis: started at {datetime.now().strftime('%H:%M:%S')}")
    adata_nv = sc.read_h5ad(fn_nv)
    adata_nv.obs["orig.ident"] = adata_nv.obs["orig.ident"].astype(str)
    sam_nv=SAM(counts = adata_nv)
    sam_nv.run(batch_key = "orig.ident")
    sam_nv.save_anndata("results/sam_objects/Nv.h5ad")

    # Step 4: Run SAMap
    print(f"[Step 4] Run SAMap, started at {datetime.now().strftime('%H:%M:%S')}")
    sams = {'Ac':sam_ac, 'Hv':sam_hv, 'Nv':sam_nv}
    sm = SAMAP(sams, f_maps = 'samap_directory/maps/', names = gene_names)
    sm.run(pairwise=True)

    # Save SAMap result
    print(f"[Step 5] Save SAMap result, started at {datetime.now().strftime('%H:%M:%S')}")
    save_samap(sm, "results/samap_objects/samap_result.pkl")

# Standard Python main
if __name__ == "__main__":
    main()