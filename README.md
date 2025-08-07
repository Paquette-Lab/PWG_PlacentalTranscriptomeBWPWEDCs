This repository contains code from the Paquette Lab related to the manuscript, “The Placental Transcriptome Serves as a Mechanistic Link between Prenatal Phthalate Exposure and Placental Efficiency”. This documents presents a summary of code document names and purpose.

Phthalate and Covariate Data Cleaning: 

•	Phthalate_Covar_RNAseq_DataCleaning_N222.Rmd
  	This code clean covariate data, RNAseq data, and phthalate data resulting in a final N=222 for phthalate RNAseq analyses.
•	DEHP_DataCleaning.Rmd
    This code cleans and performs the geometric mean for DEHP data which was processed separately and then added back with the other phthalates data. It includes cleaning for the full and sex-stratified subsets. 
•	Female_Phthalate_Covar_RNAseq_DataCleaning_N109.Rmd
    This code clean covariate data, RNAseq data, and phthalate data for a female subset of the sample resulting in a final N=109 for the female-stratified phthalate RNAseq analyses.
•	Male_Phthalate_Covar_RNAseq_DataCleaning_N113.Rmd
o	This code clean covariate data, RNAseq data, and phthalate data for a male subset of the sample resulting in a final N=113 for the male-stratified phthalate RNAseq analyses.
Differential Expression Analysis Phthalates: 
•	Phthalate_RNAseqAnalysis_N222.Rmd
o	This code performs the RNAseq analysis for associations between phthalates and the placental transcriptome.
•	Female_Phthalate_RNAseqAnalysis_N109.Rmd
o	This code performs the RNAseq analysis for associations between phthalates and the placental transcriptome in a female-stratified subset (N=109) of the sample.
•	Male_Phthalate_RNAseqAnalysis_N113.Rmd
o	This code performs the RNAseq analysis for associations between phthalates and the placental transcriptome in a female-stratified subset (N=113) of the sample.
Phthalates Mixture Analysis: 
•	Phthalates_QGComp_N222.Rmd
o	This code runs the quantile g-computation mixtures analysis for 15 metabolites to determine the association with placental gene expression in the full sample (N=222).
•	Phthalate_WGCNA_QGComp_N222.Rmd
o	This code runs the quantile g-computation mixtures analysis for 15 metabolites to determine the association with placental gene expression summarized as WGCNA gene modules in the full sample (N=222).
•	QGComp_Phthalate_Interpretation_N222.Rmd
o	This code evaluates QGComp mixture results for significance in the full sample (N=222) and saves/plots results. 
•	Male_Phthalate_QGComp_N113.Rmd
o	This code runs the quantile g-computation mixtures analysis for 15 metabolites to determine the association with placental gene expression in the male-stratified sample (N=113).
•	Female_Phthalate_QGComp_N109.Rmd
o	This code runs the quantile g-computation mixtures analysis for 15 metabolites to determine the association with placental gene expression in the female-stratified sample (N=109).
•	Sex-stratified_Phthalate_QGComp_Interpretation.Rmd
o	This code evaluates sex-stratified mixture results for significance in the male (N=113) and female (N=109) subsets. 
Differential Expression Analysis of placental weight, birthweight adjusted for placental weight, and birthweight:placental weight ratio
•	PW_BWadj_ratio_analyses.Rmd
o	This code a) cleans covariate and RNA-sequencing data, b) conducts descriptive statistics of the population characteristics for the N=253 participants included in the placental efficiency analyses, c) visualizes and evaluates the relationship between placental weight (PW) and birthweight (BW), d) evaluates and visualizes differential gene expression associated with PW, BW adjusted for a nonlinear association with PW (BWadj) and the BW:PW ratio, and e) evaluates and visualizes associations between WGCNA modules and PW, BWadj, and BW:PW ratio. 
Weighted Gene Co-expression Network Analysis (WGCNA)
•	WGCNA/cqn_prep_code.R
o	This code calculates the length and GC content information needed to conduct Conditional Quantile Normalization (CQN) prior to conducting WGCNA.
•	WGCNA/cqn_for_wgcna.Rmd
o	This code generates the CQN-normalized gene expression data used in the WGCNA analysis.
•	WGCNA/conductWGCNA.Rmd
o	This code conducts WGCNA. 
•	WGCNA/wgcna_moduleInterpretation.Rmd
o	This code is used to interpret the WGCNA modules by identifying hubgenes by correlation analysis and KEGG pathways enriched in each module by over-representation analysis.  
WGCNA-Phthalates: 
•	Phthalates_WGCNA.Rmd
o	This code evaluates the association between GAPPS WGCNA modules and phthalates.
Linear Models- Phthalates and Placental Efficiency
•	Phthalates_BW_BWPW_LinearModels.Rmd
o	Here we are using linear models to identify associations between phthalates, BW and BW:PW ratio
Meet in the Middle Analysis
•	MeetintheMiddle_Phthalates_PlacentalEfficiency.Rmd
o	This code performs the meet in the middle (overlap) analysis for the phthalate and placental efficiency analyses (BWadj and BW:PW).
•	Module_DEG_Overlap.Rmd
o	In this code, we take the genes associated with Phthalates, Birthweight(adj), or BW:PW (DEGs=FDR<0.05) and overlap them with the genes inside modules that are either 1) Meet in the middle modules or 2) Modules that are significant in the high-dimensional mediation analysis. The resulting overlap is used to create a supplemental table.
•	mitm_circos_plot.R
o	This code generates the Circos plot of differentially expressed genes (FDR<0.05) associated with phthalate metabolites, BWadj, or BW:PW ratio and highlights meet-in-the-middle genes with arcs/links that are color-coded according to the phthalate metabolite in question. 
•	mitm_wgcna_and_pathways.Rmd
o	This code generates two panels of a figure: one depicting the association between WGCNA modules and phthalate metabolites, DEHP, and the mixture, as well as measures of placental efficiency (BW:PW and BWadj) and the other depicting the KEGG pathways enriched among each WGCNA modules’ member genes. 
Formal Mediation Methods
•	mediation_mitm_genes.Rmd
o	Using the meet in the middle (MITM) genes that were associated with at least 1 phthalate (FDR<0.05) and BWadj (FDR<0.05), this code evaluates the role of the gene as a mediator individually using linear models adjusted for covariates. 
•	mediation_mitm_wgcna.Rmd
o	Using the MITM WGCNA modules that were associated with at least 1 phthalate metabolite (p<0.05) and BWadj or BW:PW ratio (p<0.05), this code evaluates the role of each module as a mediator individually using linear models adjusted for covariates. 
•	mediation_hima_genes.Rmd
o	Using the placental transcriptome as a high-dimensional mediator, this code employs the hima package to conduct high dimensional mediation (HIMA) for the relationship between each phthalate that was marginally associated with BWadj or BW:PW ratio. 
•	mediation_hima_wgcna.Rmd
o	Using the WGCNA coexpression modules as a high-dimensional mediator, this code employs the hima package to conduct high dimensional mediation (HIMA) for the relationship between each phthalate that was marginally associated with BWadj or BW:PW ratio. 
<img width="468" height="624" alt="image" src="https://github.com/user-attachments/assets/0c2bab72-ae13-4083-9fe9-93120fe4d3f3" />
