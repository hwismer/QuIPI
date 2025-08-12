import pandas as pd
import pyarrow.feather as feather
import csv
import scanpy as sc


## HUMU
flow_cats = ["Species", "Human Archetype / Mouse Tumor Line", "Human Indication / Mouse Tumor Line" ]
flow_scores = ['CD45_Live', 'Tcells_Live','Myeloid_Live', 'Stroma_Live', 
               'Treg_Live', 'CD4_Tcells', 'CD8_Tcells', 'Macrophages_Myeloid', 
               'Monocytes_Myeloid', 'cDC1_Myeloid', 'cDC2_Myeloid', 'Ki67_Tumor', 
               'NKcells_Live', 'Bcells_Live','Neutrophiles_Live', 'Exhaustion_CD8']

categoricals = ['orig.ident','nCount_RNA','nFeature_RNA','Tumor Line','Experiment','Age','Drink',
                      'Diet','Drink_Diet','Microbiota','Tumor_site','Compartment',
                      'Coarse Annotation', 'Fine Annotation']

categoricals_opts = ['Tumor Line','Compartment','Coarse Annotation', 'Fine Annotation']

categorial_opts_dict = {"Tumor Line": ['B16-F10', 'LLC', 'MC38', 'RENCA', '4T1', 'CT26', 'KPCY-F', 'KPCY-C2', 'PYMT-F'],
                        "Experiment" : ['B16 ABx/HFD','LLC ABx/HFD','MC38 ABx/HFD','RENCA ABx/HFD','4T1 ABx/HFD',
                                        'CT26 ABx/HFD','PDAC 6','LLC Young/Aged','MC38 Young/Aged','B16 Young/Aged 2',
                                        '4T1 Young/Aged','RENCA Young/Aged','CT26 Young/Aged','SCPYMT1','SCPYMT2',],
                        "Compartment" : ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma"],
                        "Coarse Annotation" : ['CD8_T_cells','TAMs','Fibroblasts','CD4_Tconv','Tumor',
                                          'CD4_Treg','cDC','Mono_Mac','Monocytes','pDC','gdT_other_T_cells',],
                        "Fine Annotation" : ['CD8_eff','Doublets_TAM_prolif–Tex','Fibro_3','CD4_Teff_aging','Tumor_1','CD4_Treg','TAM_C1q_Apoe','CD4_Teff','CD4_Teff_naive','CD8_naive','cDC1','MonoMac_Cd14_Ccr2_C1q','CD8_eff_ex','Mono_ISG','TAM_Mrc1_Cd163','Melano_1','Myofibro_1','Fibro_2','TAM_Ms4a7','cDC2','Mono_Cd14_Ccr2','cDC2_Mgl2','TAM_prolif','TAM_Atp6v0d2_Trm_like','DC_prolif','TAM_Arg1','Melano_2','cDC3','Mo_DC','Tumor_Prolif_2','CD4_Klf2','Tumor_2','DC_ISG','TAM_prolif_Mgl2_Cst3','MonoMac_ISG_C1q','pDC','Tumor_Prolif_1','TAM_prolif_matrix','Melano_Prolif','gdT_other','DC_Cst3']}

#adata = sc.read_h5ad("./quipi_humu_data/quipi_humu_adata.h5ad", backed="r")

## QUIPI

humu_compartments = ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma"]
quipi_compartments = ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma", "Live"]

quipi_cats = ["indication","compartment", "archetype"]
quipi_cats_opts = ["Indication", "Compartment", "Archetype"]
quipi_cats_dict = {"Indication" : "indication", "Compartment":"compartment", "Archetype" : "archetype", "---" : "compartment"}

with open("./quipi_humu_data/quipi_cols.csv", 'r') as file:
    reader = csv.reader(file)
    # Read the first row only (header)
    quipi_all_columns= next(reader)

quipi_cats = ['patient','sample_name','indication', 'archetype', 'compartment', 'sample_type', 'sample_type_cat', 'x_umap1', 'x_umap2']
quipi_genes = list(set(quipi_all_columns) - set(quipi_cats))



## COLORS

compartment_colors = {"T Cell": "#009292",
                      "Myeloid" : "#FF6DB6",
                      "Stroma" : "#006DDB",
                      "Tumor" : "#DB6D00",
                      "T-Reg" : "#B66DFF",
                      "Live" : "#004949"}

compartment_order = ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma"]

indic_to_color = {'LUNG':'rgb(15, 8, 58)',
                  'HEP':'rgb(55, 77, 161)', 
                  'ADR':'rgb(248, 156, 116)',
                  'GBM':'rgb(69, 70, 75)', 
                  'CRC':'rgb(74, 135, 145)',
                  'BRC':'rgb(158, 185, 243)',
                  'KID':'rgb(124, 20, 20)',
                  'MEL':'rgb(15, 2, 2)',
                  'PNET':'rgb(82, 26, 113)',
                  'GYN':'rgb(153, 99, 87)',
                  'HNSC':'rgb(39, 111, 36)',
                  'SI':'rgb(220, 176, 242)',
                  'SRC':'rgb(281, 23, 110)',
                  'GALL':'rgb(246, 207, 113)',
                  'PDAC':'rgb(165, 171, 44)',
                  'BLAD':'rgb(212, 160, 64)'}

indication_order = ['LUNG', 'HEP', 'ADR', 'GBM', 'CRC', 'BRC', 'KID', 'MEL',
                    'PNET', 'GYN', 'HNSC', 'SI', 'SRC', 'GALL', 'PDAC', 'BLAD']

colors_pancan = {'Unclassified' : '#d6d893',
                 'IR CD8 Mac' : '#ed1e21',
                 'IR CD8 Mono' : '#f06ba8',
                 'IR CD4 Mac' : '#7e1515',
                 'IS CD8' : '#128042',
                 'IS CD4' : '#98ca3a',
                 'TC Mac' : '#2c276b',
                 'TC DC' : '#4a87c7',
                 'MC DC2' : '#7f7f7f',
                 'MC DC1' : 'black',
                 'ID CD4 Mac' : '#fdd80d',
                 'ID Mono' : '#b8882c',
                 'ID CD8 Mac' : '#f68c20'}

mouse_colors =  {'4T1' : '#1E90FF',
                 'B16-F10' : '#8B7500',
                 'CT26' : '#A020F0',
                 'ID8' : '#7F7F7F',
                 'KPCY-C2' : '#006400',
                 'KPCY-C5' : '#00CD00',
                 'KPCY-F' : '#EEC591',
                 'LLC' : '#008B00',
                 'MC38' : '#CD3700',
                 'PyMT-F' : '#DA70D6',
                 'PyMT-M' : '#8B4789',
                 'RENCA' : '#104E8B',
                 'YUMM1.G1' : '#FFC125',
                 'YUMM3.3' : '#FF69B4',
                 'YUMM5.2' : '#6B8E23'}

groupby_colors = {"Tumor Line" : mouse_colors,
                  "Compartment" : compartment_colors,
                  "Coarse Annotation": None,
                  "Fine Annotation": None}


quipi_colors = {"compartment" : compartment_colors,
                "indication" : indic_to_color,
                "archetype" : colors_pancan}



humu_orders = {"Compartment" : compartment_order,
               "Tumor Line" : [],
               "Coarse Annotation" : [],
               "Fine Annotation" : []}

quipi_compartment_order = ["Tumor", "T Cell", "T-Reg", "Myeloid", "Stroma", "Live"]
archetype_order = ['Unclassified', 'IR CD8 Mac', 'IR CD8 Mono', 'IR CD4 Mac', 'IS CD8','IS CD4', 'TC Mac','TC DC','MC DC2','MC DC1','ID CD4 Mac','ID Mono','ID CD8 Mac']
sample_type_order = ["T1", "T2", "N1", "N2"]
indication_order = ['LUNG', 'HEP', 'CRC', 'KID', 'MEL', "MELB", 'PNET', 'GYN', 'HNSC', 'SRC', 'PDAC', 'BLAD', "GBM"]
mu_tumor_order =  ['4T1','B16-F10' ,'CT26' , 'ID8' ,'KPCY-C2' , 'KPCY-C5','KPCY-F','LLC' ,'MC38','PyMT-F','PyMT-M','RENCA','YUMM1.G1' ,'YUMM3.3','YUMM5.2']


humu_flow_orders = {"Species" : [],
                    "Human Archetype / Mouse Tumor Line" : archetype_order + mu_tumor_order,
                    "Human Indication / Mouse Tumor Line" : indication_order + mu_tumor_order
                    }

quipi_orders = {"indication" : indication_order,
               "sample_type" : sample_type_order,
               "archetype" : archetype_order,
               "compartment" : quipi_compartment_order}

# Get all the gene names
with open("./quipi_humu_data/quipi_humu_genes.csv", "r") as f:
    reader = csv.reader(f)
    genes = next(reader)
