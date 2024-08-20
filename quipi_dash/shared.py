from pathlib import Path

import pandas as pd


df = pd.read_pickle("./data/clean/cleaned_df_ALL.pi")
non_genes = ["Patient", "sample_name", "Indication", "Compartment", "Archetype", "X_umap1", "X_umap2"]

categoricals = ["Patient", "sample_name", "Indication", "Compartment", "Archetype"]


genes = list(set(df.columns) - set(non_genes))
#genes = ["NCAM1", "DLK1"]