import numpy as np
import pandas as pd
from scipy.stats import zscore
import shared as sh

def calculate_gene_factor_score(gene_set,compartment):
    return None

    input_arr = pd.read_csv("./quipi_log2_tpm.csv", usecols=sh.non_genes + gene_set)
    #log2_subset = sh.quipi_log2[sh.quipi_log2["archetype"] != "Unclassified"][sh.quipi_log2["compartment"] == compartment][sh.non_genes+gene_set]
    #raw_subset = sh.quipi_raw[sh.quipi_raw["archetype"] != "Unclassified"][sh.quipi_raw["compartment"] == compartment][sh.non_genes+gene_set]

    z_subset = log2_subset[gene_set].apply(zscore)
    z_subset["factor_score"] = z_subset.mean(axis=1)
    log2_subset_full = log2_subset[sh.non_genes].merge(z_subset,left_index=True,right_index=True)

    return log2_subset_full
