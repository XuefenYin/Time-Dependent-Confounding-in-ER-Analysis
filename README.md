# Time-Dependent-Confounding-in-ER-Analysis
This repository contains the R model files used to investigate how conventional exposure–response (ER) analyses—such as logistic regression, Kaplan-Meier plots, and Cox proportional hazards models—may be affected by time-dependent confounding factors, including exposure accumulation, dose modification patterns, and event onset timing. Due to the large size of the simulation outputs, we have only included results for the true exposure–response curves to allow readers to quickly run a representative simulation. This repository also serves as supplementary material for our forthcoming publication.
# A.Identifying Bias Sources in ER Analyses
#### 1 ER1_DH1_1Dose.Rmd
Defines an ER scenario (ER1) where the response is independent of drug exposure, under a **fixed dosing regimen** without any dose modifications (DH1).
#### 2 ER1_DH3_1Dose.Rmd
Defines the same exposure-independent scenario (ER1), but under an **empirical dosing regimen** with significant dosing modifications(DH3).
#### 3 ER2_ORR_DH1_1D_truncate by PFS_1000ID
Defines a **positive ER scenario** (ER2) where the response is driven by drug exposure under **fixed dosing** (DH1). The endpoint is ****overall response rate** (ORR)**, and observations are truncated by progression-free survival (PFS).
#### 4 ER2_ORR_DH2_AE_LDR_truncate by PFS_1D_1000ID.Rmd
Defines a **positive ER scenario** (ER2) under a **dynamic dosing regimen** (DH2) where plasma exposure triggers adverse events (AEs) leading to dose reductions (LDR). Endpoint: **ORR**, truncated by PFS.
#### 5 ER2_PFS_DH1_1D_1000ID.Rmd
Defines an ER scenario (ER2, higher exposure will have a longer survival time) under **fixed dosing** (DH1). Endpoint: **progression-free survival (PFS)**.
#### 6 ER2_PFS_DH2_AE_LDR_1D_1000ID.Rmd
Defines an ER scenario (ER2, higher exposure will have a longer survival time) under **dynamic dosing (DH2)** with AE-driven dose reductions. Endpoint: **PFS**.
# B.Evaluating Bias-Mitigation Strategies
## Comparing Conventional vs Modified method to derive time-dependent exposures
These files test the performance of conventional vs modified methods for deriving time-dependent exposures under dynamic dosing (DH2):
#### 1. ER1_DH2_1Dose_ORR_Conventional.Rmd
#### 2. ER2_ORR_DH2_AE_LDR_truncate by PFS_1D_1000ID.Rmd
#### 3. ER1_DH2_1Dose_Modified.Rmd
## Investigating the Influence of Dose Level and Range
These models explore how dose levels and dose ranges influence ER relationships:
#### 4. ER2_ORR_DH2_AE_LDR_truncate by PFS_2D_500ID for each arm.Rmd
#### 5. ER2_PFS_DH2_AE_LDR_2D_500ID for each arm.Rmd
#### 6. ER1_DH1_2Dose.Rmd
#### 7. ER1_DH3_2Dose.Rmd
#### 8. ER1_DH2_2Dose_ORR_TTE.Rmd
# C.Impact of Sample Size and Dropout (Censoring) Rate
These files assess how sample size and censoring rate affect bias:
#### 9.ER1_DH1_1Dose_sample_size_censor_rate.Rmd
#### 10.ER1_DH3_1Dose_sample_size_censor_rate.Rmd
# D.Establishing the True Exposure–Response Curves
These scripts were used to generate the true ER relationships for ORR and PFS:
#### 1.Explore true ER relationship_DH1
#### 2.Explore true ER relationship_DH2_LDR.Rmd
#### 3.ER2_DH2_AE_HL24_ii24.Rmd – A sample file used to extract the true ER curve under DH2
