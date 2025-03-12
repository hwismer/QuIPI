import pandas as pd
import pyarrow.feather as feather
import csv

flow_cats = ["Species", "Group"]
flow_scores = ['CD45_Live', 'Tcells_Live','Myeloid_Live', 'Stroma_Live', 
               'Treg_Live', 'CD4_Tcells', 'CD8_Tcells', 'Macrophages_Myeloid', 
               'Monocytes_Myeloid', 'cDC1_Myeloid', 'cDC2_Myeloid', 'Ki67_Tumor', 
               'NKcells_Live', 'Bcells_Live','Neutrophiles_Live', 'Exhaustion_CD8']

categoricals = ['orig.ident','nCount_RNA','nFeature_RNA','Tumor_Line','Experiment','Age','Drink',
                      'Diet','Drink_Diet','Microbiota','Housing','Tumor_site','PanCan_Compartment',
                      'Coarse_Annot', 'Fine_Annot']

categoricals_opts = ['Tumor_Line','Experiment','PanCan_Compartment','Coarse_Annot', 'Fine_Annot']

# Get all the gene names
with open("./quipi_humu_data/quipi_humu_genes.csv", "r") as f:
    genes = next(csv.reader(f))
