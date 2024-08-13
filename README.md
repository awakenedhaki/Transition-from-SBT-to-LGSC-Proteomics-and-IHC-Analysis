# Changes in the tumour microenvironment mark the transition from serous borderline tumour to low-grade serous carcinoma

Publication: https://doi.org/10.1002/path.6338

Authors: **Rodrigo Vallejos**, Almira Zhantuyakova, Gian Luca Negri, Spencer D Martin, Sandra E Spencer, Shelby Thornton, Samuel Leung, Branden Lynch, Yimei Qin, Christine Chow, Brooke Liang, Sabrina Zdravko, J Maxwell Douglas, Katy Milne, Bridget Mateyko, Brad H Nelson, Brooke E Howitt, Felix KF Kommoss, Lars-Christian Horn, Lien Hoang, Naveena Singh, Gregg B Morin, David G Huntsman, Dawn Cochrane

Corresponding author: Dawn Cochrane

The work performed in this project was primarily conducted at the BC Cancer Research Centre (BCCRC), with affiliations to the Vancouver General Hospital and the University of British Columbia. The project was overseen by Dr. Dawn Cochrane, a staff scientist in the Huntsman Lab within the Department of Molecular Oncology at the BCCRC. This project was a collaborative effort between research groups led by Dr. David Huntsman, Dr. Gregg Morin, Dr. Brad Nelson, and Dr. Brooke Howitt.

## Project Overview

Low-grade serous ovarian carcinoma (LGSC) is a rare and aggressive subtype of ovarian cancer, distinct from the more common high-grade serous ovarian carcinoma (HGSC) in terms of pathology, biology, and clinical outcomes. LGSC often arises from serous borderline ovarian tumours (SBTs), but the mechanisms underlying this transformation remain poorly understood. 

To elucidate the biology of LGSC, this project conducted a proteomic analysis of formalin-fixed, paraffin-embedded tissue samples, including LGSC (n = 11), HGSC (n = 19), and SBTs (n = 26). The proteomic data revealed that the protein expression profiles could differentiate between the histotypes of ovarian epithelial tumours. Among the differentially expressed proteins, Fibroblast Activation Protein (FAP), commonly found in cancer-associated fibroblasts, was identified as significantly more abundant in LGSC compared to SBTs.

In addition to proteomics, multiplex immunohistochemistry (mIHC) was performed to assess the immune landscape within these tumours. The LGSC samples exhibited a tumour microenvironment enriched with Tregs and M2 macrophages, in contrast to SBTs. These findings suggest that changes in the tumour microenvironment may facilitate LGSC tumourigenesis and progression, highlighting the tumour microenvironment as a potential therapeutic target in LGSC.

This repository contains the R code used for the proteomic and immunohistochemical analyses, as well as the generation of the figures presented in the publication.

## Data Availability

The proteomic data is available in ProteomeXchange at [PXD046360](https://proteomecentral.proteomexchange.org/). The mass spectrometry proteomics data have been deposited in the ProteomeXchange Consortium via the PRIDE partner repository with the dataset identifier PXD046360.

No data has been made available in this repository.

## Figures

All figures were generated using `ggplot2`, exported into a vector format using `ggsave`, and assembled into the final figures using Affinity Designer 2 (v2.4.2).

## Code

The code is organized into the following directories:

- `analysis`: Contains the R code used to analyze the data (e.g., proteome, IHC, and multiplex IHC).
- `render`: Contains the R code used to generate the figures.
- `revisions`: Contains the R code used for the revisions requested by the reviewers.

The scripts within each directory are grouped according to the figure and panel they are associated with.

## Citation

        Vallejos R, Zhantuyakova A, Negri GL, et al. Changes in the tumour microenvironment mark the transition from serous borderline tumour to low-grade serous carcinoma.  J Pathol 2024

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Session Info

```
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.0

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Vancouver
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.2.1                ggrepel_0.9.3               RColorBrewer_1.1-3         
 [4] viridisLite_0.4.2           pheatmap_1.0.12             ggforce_0.4.1              
 [7] org.Hs.eg.db_3.17.0         AnnotationDbi_1.62.1        IRanges_2.34.1
[10] S4Vectors_0.38.1            Biobase_2.60.0              BiocGenerics_0.46.0 
[13] msigdbr_7.5.1               fgsea_1.26.0                DEqMS_1.18.0 
[16] matrixStats_1.0.0           limma_3.56.2                glue_1.6.2
[19] magrittr_2.0.3              ConsensusClusterPlus_1.64.0 tibble_3.2.1
[22] purrr_1.0.1                 broom_1.0.5                 stringr_1.5.0
[25] tidyr_1.3.0                 readxl_1.4.3                forcats_1.0.0
[28] ggplot2_3.4.2               janitor_2.2.0               dplyr_1.1.2
[31] readr_2.1.4                 here_1.0.1                 

loaded via a namespace (and not attached):
 [1] DBI_1.1.3               bitops_1.0-7            rlang_1.1.1             snakecase_0.11.0       
 [5] compiler_4.3.1          RSQLite_2.3.1           mgcv_1.8-42             png_0.1-8              
 [9] vctrs_0.6.3             pkgconfig_2.0.3         crayon_1.5.2            fastmap_1.1.1          
[13] backports_1.4.1         XVector_0.40.0          labeling_0.4.2          utf8_1.2.3             
[17] tzdb_0.4.0              bit_4.0.5               zlibbioc_1.46.0         cachem_1.0.8           
[21] GenomeInfoDb_1.36.1     blob_1.2.4              BiocParallel_1.34.2     tweenr_2.0.2           
[25] parallel_4.3.1          cluster_2.1.4           R6_2.5.1                stringi_1.7.12         
[29] lubridate_1.9.2         cellranger_1.1.0        Rcpp_1.0.10             Matrix_1.5-4.1         
[33] splines_4.3.1           timechange_0.2.0        tidyselect_1.2.0        rstudioapi_0.14        
[37] codetools_0.2-19        lattice_0.21-8          withr_2.5.0             KEGGREST_1.40.0        
[41] polyclip_1.10-4         Biostrings_2.68.1       pillar_1.9.0            generics_0.1.3         
[45] vroom_1.6.3             rprojroot_2.0.3         RCurl_1.98-1.12         hms_1.1.3              
[49] munsell_0.5.0           tools_4.3.1             data.table_1.14.8       babelgene_22.9         
[53] fastmatch_1.1-3         cowplot_1.1.1           grid_4.3.1              colorspace_2.1-0       
[57] nlme_3.1-162            GenomeInfoDbData_1.2.10 cli_3.6.1               fansi_1.0.4            
[61] rematch_1.0.1           gtable_0.3.3            farver_2.1.1            memoise_2.0.1          
[65] lifecycle_1.0.3         httr_1.4.6              bit64_4.0.5             MASS_7.3-60  
```
