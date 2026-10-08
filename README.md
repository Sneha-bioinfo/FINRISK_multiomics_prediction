# Latent variable software for integrating multi-omics data for incident disease risk prediction in microbiome-based cohort studies

**Beta version**

This repository contains the FINRISK evaluation code for the Hologen deliverable *"Latent variable software for integrating multi-omics data for incident disease risk prediction in microbiome-based cohort studies."*

## About the software

The software is the extended **IntegratedLearner** framework [IntegratedLearner](https://github.com/himelmallick/IntegratedLearner), developed by Himel Mallick and Nalin Arora. The original IntegratedLearner supported cross-sectional analysis of binary and continuous outcomes. The extended version adds entirely new functionalities:

- **Survival outcomes**: multi-omics integration for prospective cohorts with time-to-event data, enabling incident disease risk prediction.
- **Multiclass outcomes**: evaluation of more than two outcome categories.

Currently, I'm working  on the  manusscript of the paper with Himel Mallick and NAKON aora, to oisnryicude the enw functionaltiies and sciemtic onbservations form the test on the firnisk cohort obstained from the Finrisk and my supervisor Leo Lahti, finalizng, data ansalsyss.

## This repository

The software was tested on the FINRISK cohort using taxonomic (metagenomic) and metabolomic profiles. Every version of the software, up to the current final version, was evaluated on FINRISK data, and the results guided its progression. This repository contains the analysis code for that testing.

## Scripts
- [https://github.com/Sneha-bioinfo/FINRISK_multiomics_prediction/blob/main/scripts/FINRISK_IL_survival.R](FINRISK_IL_survival.R)- Script for survival outcome prediction on the FINRISK data using Integratedlearner framework
- [https://github.com/Sneha-bioinfo/FINRISK_multiomics_prediction/blob/main/scripts/FINRISK_IL_survival_null.R](FINRISK_IL_survival_null.R)- Script for survival outcome prediction on negative control to test if the model accurately predicts biolgical association and not jist by chane
- [https://github.com/Sneha-bioinfo/FINRISK_multiomics_prediction/blob/main/scripts/FINRISK_IL_muticlass.R](FINRISK_IL_multiclass.R)- Script for multiclass outcome prediction on the FINRISK


## Data

FINRISK data is sensitive data and needs to be assessed using the CSC SD desktop services with permission obtained to asses the data. 


## Funding
This work is part of my (DC2) PhD project and has been funded by the *European Union under Grant Agreement 101169005 (Hologen Consortium 2025)
<img width="1205" height="122" alt="image" src="https://github.com/user-attachments/assets/7eea7d1a-0adc-4004-93c2-6a89e606984b" />











