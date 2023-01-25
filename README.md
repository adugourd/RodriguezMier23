# Inference of signaling networks from multi-condition and multi-timepoint data



This repository contains examples of how CORNETO can be used to infer signaling networks from multiple independent interventions (e.g multiple perturbations with drugs), or from longitudinal data (e.g RNA-seq measured at different time-points after a drug perturbation).



## PANACEA dataset

CORNETO subsumes the functionality of CARNIVAL and extends it by modelling the problem of inferring the network given a Prior Knowledge and a set of observations (e.g transcription factor activities per condition) by using Mixed Integer Optimization on top of Network Flows. By leveraging information about multiple independent interventions at the same time, the intracellular signal of cells can be much better characterized.

We illustrate the use of CORNETO with a prostate cancer cell line (DU-145) from PANACEA [1]. The drug has been perturbed with 32 different kinase inhibitors, and RNA-seq data is available for each perturbation. Using this data, we estimated transcription factor (TFs) activities for each perturbation with DecoupleR [2], and we inferred the sparsest sub-graph from Omnipath [3] that is able to recapitulate the signaling events from perturbations to TFs.

Code and data are available in the [panacea](https://github.com/saezlab/RodriguezMier23/tree/main/panacea) folder. The main scripts are:

- [`preprocess.R`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/processing.R): R script to process the GSE186341 DU-145 transcriptomics data, estimate TF activities and prepare input data for the method.
- [`panacea_multicondition.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/panacea-multicondition.ipynb): Python notebook using CORNETO to infer the intracellular signaling of DU-145
- [`process_results.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/panacea/process_results.ipynb): Python notebook that reads output results from multiple runs and summarizes the results.



## A375 BRAF V600E longitudinal RNA-seq data

We used data from Gerosa et al. [4] to illustrate how CORNETO can be used in a different complete setting, without changing the underlying method. In this dataset, A375 Braf V600E mutant cells were treated with a BRAF inhibitor for 24h, and then stimulated with EGF. Transcriptomics data was obtained after EGF stimulation at different timepoints (0h, 0.5h, 1h, 2h, 3h, 4h, 8h). 

We estimated TF activities from changes in gene expression between consecutive time points (0.5h vs 0h, 1h vs 0.5h, ...), and then we use CORNETO to sequentially infer signaling networks that explain the TF activities from EGFR.

The following scripts are available in the [A375-BRAFi](https://github.com/saezlab/RodriguezMier23/tree/main/A375-BRAFi) folder:

- [`preprocess.R`](https://github.com/saezlab/RodriguezMier23/blob/main/A375-BRAFi/preprocess.R): R script that processes the transcriptomics data, estimates differential expression between timepoints and then TF activities.
- [`a375-longitudinal.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/A375-BRAFi/a375-longitudinal.ipynb): Implements the sequential loop with CORNETO to infer signaling networks that explains the observed changes in TFs. Network for the diff. time point _t+1_ is fitted to the data minimizing the structural distance of the network with respect the network for _t_.
- [`process_results.ipynb`](https://github.com/saezlab/RodriguezMier23/blob/main/A375-BRAFi/process_results.ipynb): Process and analyzes the results and generates the plots.



## References

1. Douglass Jr, Eugene F., et al. "A community challenge for a pancancer drug mechanism of action inference from perturbational profile data." *Cell Reports Medicine* 3.1 (2022): 100492.
2. Badia-i-Mompel, Pau, et al. "decoupleR: ensemble of computational methods to infer biological activities from omics data." *Bioinformatics Advances* 2.1 (2022): vbac016.
3. Türei, Dénes, Tamás Korcsmáros, and Julio Saez-Rodriguez. "OmniPath: guidelines and gateway for literature-curated signaling pathway resources." *Nature methods* 13.12 (2016): 966-967.
4. Gerosa, Luca, et al. "Receptor-driven ERK pulses reconfigure MAPK signaling and enable persistence of drug-adapted BRAF-mutant melanoma cells." Cell systems 11.5 (2020): 478-494.



