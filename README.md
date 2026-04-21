# Creating an atlas of CD8 T Cell sc data and using said data for mapping their inhibitory receptor distribution

#### David Theodor Nimrichtr within the [Oxford WIMM Sharma Group](https://www.imm.ox.ac.uk/research/research-groups/sharma-group-mechanistic-t-cell-genomics)

## Documents
- [draft project description](https://unioxfordnexus-my.sharepoint.com/:w:/r/personal/rdom0281_ox_ac_uk/Documents/project%20plans/student%20project/David%20project%20draft.docx?d=wd93b3192ff1d45a6a30ebb6606caca1d&csf=1&web=1&e=NgBtd4)
- [todo list](https://unioxfordnexus-my.sharepoint.com/:w:/r/personal/rdom0281_ox_ac_uk/Documents/project%20plans/student%20project/David-%20to%20do%20list.docx?d=wad99e5995b984b7da27bf1fbebd36035&csf=1&web=1&e=a8Qa2K)

# Checklist

1. [ ] **Generation of the T cell atlas across diseases and tissues (excel sheet, in progress)**

    1. [ ] create checklist
    2. [ ] Check the platform for sequencing, the reference genome used (if they mapped differently, highlight them), gene names (if they are from different nomenclature then you cannot map them- HUGO is preferred/ HGNC or the Ensembl gene IDs needs to be mapped to HUGO). 

    3. [ ] Check the meta data, especially column name across studies, such as batch, donor, disease names and collected tissues etc. They need to be concordant with each other. In the excel sheet I am preparing I am making them concordant so, you need to annotate them accordingly. (There are examples within the GitHub page I sent you above) 
    


2. [ ] **Inhibitory receptor analyses (list of receptors: [data/protein_list.csv])**
	1. [x] Start with smaller scale CD8+ T cell data. (You can have CD4+ T cells as well but in this case, we are not really looking for them. You can get rid of them once you have log1p transformed data: You should have cells with CD8 expression and no CD4 expression)  
	2. [x] Check QC and filter your dataset with good quality cells  
	3. [x] Optimise your parameter and collect your statistics about the data (important to justify why you have used specific number of clusters etc., check what are the highly expressed genes does that make biological sense) 
	4. [ ] Once you have T cell state clusters, use markers to understand what biology each cluster is representing. (You can use the markers in this [paper](https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full) Figure 1D). But please also check literature for how you can better name them. You can use the paper of the data.  
	5. [ ] Then, you can start to check where the inhibitory receptors are expressed (as a sanity check, you need to see them in effector, cytotoxic and exhausted cells and not in naive cells). You can check what meta data you have (disease, tissue collected, age etc) and then try to compare if they have an effect on any receptor expression: Navigate the data with statistics (consider using pseudobulking as well): 
		1. [ ] Boxplots: x-axis diseases and y-axis scaled expression of receptors among cells. You can make statistical test between them.  
		2. [ ] Odds ratio: You can binarise expression as exist or not exist and test with Fisher test across tissues/diseases. 
		3. [ ] More sophisticated a linear model: You can use expression value as y, and x will be any meta data or combination of them, and you can apply OLS and p-value correction to collect significant associations.  
	6. [ ] Later step would be transcriptional regulation. I will explain later.

Filter out ACTB, parameters

0. Use FindAllMarkers - most expressed for each cluster -> dot plot
1.	Inhibitory receptors - they will be not highly expressed - join them into families?

2. Celltypist and compare against the dataset
Add detailed gene annotation
Only T cells
⦁	only 8 and not 4
⦁	Only 3 (both 4 and 8)
Inhibitory receptors - they will be not highly expressed - join them into families?


2. Cellchat
⦁	Liana_path julio group
⦁	Same tissue, same disease, different chats
⦁	Cancer vs autoimmune
⦁	https://github.com/sqjin/CellChat
⦁	https://github.com/saezlab/liana-py
https://doi.org/10.1038/s41467-022-30755-0
https://doi.org/10.1038/s41556-024-01469-w
PBMC cellchat
Finad a dataset with TIL + tumour tissue -> cellchat

3. Calculate cell cycle stage
⦁	With seurat (in R?)
sc.tl.rank_genes_groups(adata, groupby=clustered_name)
https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pl.rank_genes_groups_dotplot.html





# Data Sources
- https://www.nature.com/articles/s41592-024-02530-0
	- https://zenodo.org/records/13382785


# LabBook

### 20. 3. 2026
- [x] created the labbook
- [x] created the bioreport library
- [ ] adding more granular gene filtering
	- [x] copied exclusion patterns from sumana, used cc genes from here https://www.nature.com/articles/s41592-025-02793-1#Sec30 - top 50 genes from the gene expression programmes CellCycle-G2M, CellCycle-S, CellCycle-Late-S
- [x] creating parameter iteration script
### 24. 3. 2026
- [x] improving the test_pipeline algorithm
	- [x] dont selecting cc genes in a smarter way
	- [x] changing order: the full gene exclusion is final step of QC
- [x] meeting with Cansu and Sumana, summary bellow
### 26. 3. 2026
- [x] added ACTB to filtering based on SUmanas recommendation

### *easter break*

### 17. 4. 2026
- [x] working findallmarkers implementation (launching the function in R)

### 21. 4. 2026
- [x] determined best parameters based on my own judgement evaluating the [parameter iteration report](cd8atlas/output/parameter_iteration/parameter_iteration_20260421_035352_report.pdf)
- [x] implemented state determination based on https://www.biorxiv.org/content/10.1101/2024.07.17.603780v1.full
- [ ] 





# Meetings
### Tasks after meeting on 24. 3. 2026

inhibitory club - end of april

Filter out ACTB, parameters

0. Use FindAllMarkers - most expressed for each cluster -> dot plot
1. Filter Inhibitory receptors - they will be not highly expressed - join them into families?

2. Celltypist and compare against the dataset
Add detailed gene annotation
Only T cells
⦁	only 8 and not 4
⦁	Only 3 (both 4 and 8)
Inhibitory receptors - they will be not highly expressed - join them into families?


2. Cellchat
⦁	Liana_path julio group
⦁	Same tissue, same disease, different chats
⦁	Cancer vs autoimmune
⦁	https://github.com/sqjin/CellChat
⦁	https://github.com/saezlab/liana-py
https://doi.org/10.1038/s41467-022-30755-0
https://doi.org/10.1038/s41556-024-01469-w
PBMC cellchat
Finad a dataset with TIL + tumour tissue -> cellchat

3. Calculate cell cycle stage
⦁	With seurat (in R?)
sc.tl.rank_genes_groups(adata, groupby=clustered_name)
https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pl.rank_genes_groups_dotplot.html