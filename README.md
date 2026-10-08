# Latent variable software for integrating multi-omics data for incident disease risk prediction in microbiome-based cohort studies

**Beta version**

This repository contains the FINRISK analysis scripts used to test the beta version of the Hologen deliverable *"Latent variable software for integrating multi-omics data for incident disease risk prediction in microbiome-based cohort studies."*

## Background

The latent variable software is the extended [IntegratedLearner](https://github.com/himelmallick/IntegratedLearner) framework, developed by Himel Mallick and Nalin Arora and evaluated on the FINRISK cohort by me. The original IntegratedLearner supported cross-sectional analysis of binary and continuous outcomes. The extended version adds entirely new functionalities:

- **Survival outcomes**: multi-omics integration for prospective cohorts with time-to-event data, enabling incident disease risk prediction.
- **Multiclass outcomes**: evaluation of more than two outcome categories.

A manuscript presenting these new functionalities and the FINRISK findings is being prepared as a collaborative effort with my supervisor, Leo Lahti, and with Himel Mallick and Nalin Arora from the original IntegratedLearner team.

## This repository

The software was tested on the FINRISK cohort using taxonomic (microbiome data )and metabolomic profiles. Every version of the software, up to the current final version, was evaluated on FINRISK data, and the results guided its progression. This repository contains the analysis scripts for the same.

## Scripts

- [`FINRISK_IL_survival.R`](scripts/FINRISK_IL_survival.R): survival outcome prediction on FINRISK data using the IntegratedLearner framework.
- [`FINRISK_IL_survival_null.R`](scripts/FINRISK_IL_survival_null.R): negative control analysis for the survival model, testing that predictions reflect true biological association and not chance.
- [`FINRISK_IL_multiclass.R`](scripts/FINRISK_IL_muticlass.R): multiclass outcome prediction on FINRISK data.

## Data

FINRISK data is sensitive and thus not included in this repository. The data needs to be accessed through the CSC SD Desktop service, with permission granted to access the data.

## Funding

This work is part of my (DC2) PhD project and has been funded by the European Union under Grant Agreement 101169005 (Hologen Consortium 2025).

<img width="1205" height="122" alt="Hologen funding banner" src="https://github.com/user-attachments/assets/7eea7d1a-0adc-4004-93c2-6a89e606984b" />
