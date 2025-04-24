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

categorial_opts_dict = {"Tumor_Line": ['B16F10', 'LLC', 'MC38', 'RENCA', '4T1', 'CT26', 'KPCY-F', 'KPCY-C2', 'PYMT-F'],
                        "Experiment" : ['B16 ABx/HFD','LLC ABx/HFD','MC38 ABx/HFD','RENCA ABx/HFD','4T1 ABx/HFD',
                                        'CT26 ABx/HFD','PDAC 6','LLC Young/Aged','MC38 Young/Aged','B16 Young/Aged 2',
                                        '4T1 Young/Aged','RENCA Young/Aged','CT26 Young/Aged','SCPYMT1','SCPYMT2',],
                        "PanCan_Compartment" : ['T_cells','Myeloids','Stroma','Tumor','Treg',],
                        "Coarse_Annot" : ['CD8_T_cells','TAMs','Fibroblasts','CD4_Tconv','Tumor',
                                          'CD4_Treg','cDC','Mono_Mac','Monocytes','pDC','gdT_other_T_cells',],
                        "Fine_Annot" : ['CD8_eff','Doublets_TAM_prolif–Tex','Fibro_3','CD4_Teff_aging','Tumor_1','CD4_Treg','TAM_C1q_Apoe','CD4_Teff','CD4_Teff_naive','CD8_naive','cDC1','MonoMac_Cd14_Ccr2_C1q','CD8_eff_ex','Mono_ISG','TAM_Mrc1_Cd163','Melano_1','Myofibro_1','Fibro_2','TAM_Ms4a7','cDC2','Mono_Cd14_Ccr2','cDC2_Mgl2','TAM_prolif','TAM_Atp6v0d2_Trm_like','DC_prolif','TAM_Arg1','Melano_2','cDC3','Mo_DC','Tumor_Prolif_2','CD4_Klf2','Tumor_2','DC_ISG','TAM_prolif_Mgl2_Cst3','MonoMac_ISG_C1q','pDC','Tumor_Prolif_1','TAM_prolif_matrix','Melano_Prolif','gdT_other','DC_Cst3']}

# Get all the gene names
with open("./quipi_humu_data/quipi_humu_genes.csv", "r") as f:
    reader = csv.reader(f)
    genes = next(reader)
