# Oxytocin_FND

Here we provide the explanatory code for the manuscript "Salivary Oxytocin and Amygdalar Alterations in Functional Neurological Disorders - Weber & Stoffel et al., 2024"

This code works with example data but won't produce the same results as reported in the paper, as the data can only be shared on request. This code should give a braod overview on how all the analyses in the manuscript were conducted. 

## 1. Oxytocin, OXTR Methylation and Genotyping Data
The paper contains various analyses all based on regression models conducted in R. The code here provides example data with the corresponding code in order to facilitate reproducibility of our analyses. For further questions, feel free to contact: samantha.weber@bli.uzh.ch or natascha.stoffel@unifr.ch


## 2. Functional Connectivity
The code for functional connectivity and to calculate the correlation between FC and Oxytocin levels is found in the folder "Functional Connectivity". Figures have been created using BrainNet Viewer (https://www.nitrc.org/projects/bnv/). Edges and Node Files to use as input are created within the script "FunctionalConnectivity_v2.m".

## 3. Structural Analysis
Amygdalar volume was extracted using CAT12 toolbox (https://neuro-jena.github.io/cat/). Subject-wise estimates of mean amygdala volume were extracted according to the standard procedure as described in the Manual. These estimates were further used in the R code to examine the association between amygdalar volume and oxytocin, OXTR methylation and OXTR genotype. 

