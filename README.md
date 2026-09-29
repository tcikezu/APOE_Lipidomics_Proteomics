# Linda_Lipidomics_Proteomics

## System Requirements
- R (>=4.2.2)
- The packages in Analysis/libraries.r:
  - SummarizedExperiment
  - ggplot2
  - patchwork
  - cowplot
  - ComplexHeatmap
  - dplyr
  - magrittr 
  - WGCNA (1.72.5)
  - org.Hs.eg.db (3.16.0)
  - limma
  - rrcov
  - msigdbr (7.5.1)
- For complete version information see sessionInfo.txt
- No non-standard hardware required.

## Installation guide
- To install R, see `https://cloud.r-project.org/`. For installing individual packages, we recommend using the `librarian` package (see `https://www.rdocumentation.org/packages/librarian/versions/1.8.1`). Run `librarian::shelf` on all packages listed in `Analysis/libraries.r`.

## Demo
- Run `Analysis/wgcna.r`
- Expected output are results in the folders `Results-WGCNA/lipids`, `Results-WGCNA/proteins` and `Results-WGCNA/lipid-proteins`

## Instructions
Codes needed to reproduce WGCNA analyses:
- Analysis/data.r
- Analysis/wgcna.r
- Analysis/colors.r

Data needed:
- Data/lipid\_demo.xlsx
- Data/proteins\_demo.xlsx
