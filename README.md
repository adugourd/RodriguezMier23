# Inference of signaling networks from multi-condition and multi-timepoint data



This repository contains examples of how CORNETO can be used to infer signaling networks from multiple independent interventions (e.g multiple perturbations with drugs), or from longitudinal data (e.g RNA-seq measured at different time-points after a drug perturbation).



## PANACEA dataset

CORNETO subsumes the functionality of CARNIVAL and extends it by modelling the problem of inferring the network given a Prior Knowledge and a set of observations (e.g transcription factor activities per condition) by using Mixed Integer Optimization on top of Network Flows. By leveraging information about multiple independent interventions at the same time, the intracellular signal of cells can be much better characterized.

We illustrate the use of CORNETO with a prostate cancer cell line (DU-145) from PANACEA [1]. The drug has been perturbed with 32 different kinase inhibitors, and RNA-seq data is available for each perturbation. Using this data, we estimated transcription factor (TFs) activities for each perturbation with DecoupleR [2], and we inferred the sparsest sub-graph from Omnipath [3] that is able to recapitulate the signaling events from perturbations to TFs.

Code and data are available in the [panacea](https://github.com/saezlab/RodriguezMier23/tree/main/panacea) folder. The main scripts are:

- [`preprocess.R`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/processing.R): R script to process the GSE186341 DU-145 transcriptomics data, estimate TF activities and prepare input data for the method.
- [`panacea_multicondition.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/panacea-multicondition.ipynb): Python notebook using CORNETO to infer the intracellular signaling of DU-145
- [`process_results.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/process_results.ipynb): Python notebook that reads output results from multiple runs and summarizes the results.



## References

1. Douglass Jr, Eugene F., et al. "A community challenge for a pancancer drug mechanism of action inference from perturbational profile data." *Cell Reports Medicine* 3.1 (2022): 100492.
2. Badia-i-Mompel, Pau, et al. "decoupleR: ensemble of computational methods to infer biological activities from omics data." *Bioinformatics Advances* 2.1 (2022): vbac016.
3. Türei, Dénes, Tamás Korcsmáros, and Julio Saez-Rodriguez. "OmniPath: guidelines and gateway for literature-curated signaling pathway resources." *Nature methods* 13.12 (2016): 966-967.



