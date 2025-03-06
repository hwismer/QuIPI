import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import csv


flow_cats = ["Species", "Group"]
flow_scores = ['CD45_Live', 'Tcells_Live','Myeloid_Live', 'Stroma_Live', 
               'Treg_Live', 'CD4_Tcells', 'CD8_Tcells', 'Macrophages_Myeloid', 
               'Monocytes_Myeloid', 'cDC1_Myeloid', 'cDC2_Myeloid', 'Ki67_Tumor', 
               'NKcells_Live', 'Bcells_Live','Neutrophiles_Live', 'Exhaustion_CD8']