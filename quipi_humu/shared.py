import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import csv
import scanpy as sc


flow_cats = ["Species", "Group"]
flow_scores = ['CD45_Live', 'Tcells_Live','Myeloid_Live', 'Stroma_Live', 
               'Treg_Live', 'CD4_Tcells', 'CD8_Tcells', 'Macrophages_Myeloid', 
               'Monocytes_Myeloid', 'cDC1_Myeloid', 'cDC2_Myeloid', 'Ki67_Tumor', 
               'NKcells_Live', 'Bcells_Live','Neutrophiles_Live', 'Exhaustion_CD8']

adata_vars = list(pd.read_feather("./quipi_humu_data/quipi_humu_vars.feather", columns=["gene"])["gene"])
categoricals = ["Tumor_Line", "Experiment", "Compartment","PanCan_Compartment","Coarse_Annot","Fine_Annot"]
#adata = sc.read_h5ad("./quipi_humu_data/250305_Tumor_Combined_clean_HuMu.h5ad", backed="r")