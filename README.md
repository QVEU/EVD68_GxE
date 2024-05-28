# EVD68_GxE

## This repository serves as a home for the pipelines for both cell-state and viral adaptation analysis of Enterovirus D68 (EV-D68). 

### RNAseq
Cell/temperature/infection-dependent gene expression analysis is in the RNAseq folder, where RNAseq data was analyzed using Kallisto and Sleuth packages, followed by gene p-value and beta-value thresholding, and connections formed between genes via Cytoscape. Top gene hits were used for downstream molecular biology assays.

### Surfseq and Liquid Handling
High-thoroughput passaging of two EV-D68 strains in three cell types at two temperatures, in octuplet, was handled by an open-source, automated liquid handler built by OpenTrons. The protocol and robot commands are both in the liquid_handling folder. Modify locations of input fastq files and reference files in Surfseq.sh then run from repo directory.

Analysis of whole-genome amplicons generated from RNA extracted from each passage replicate was performed by Surfseq, which uses the high-replicate, consensus-sequence produced from each sample. This was based on the amplicon-based protocol ARTIC (https://artic.network), which was used here for tracking the accumulation of mutations that potential allow for condition-specific adaptations. Resulting mutations are combined into a tidy dataframe for downstream analysis. 