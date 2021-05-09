# Data and Code from the "Biology Meets Subduction" microbial diversity paper across the subduction zone of Costa Rica

[![DOI](https://zenodo.org/badge/199023313.svg)](https://zenodo.org/badge/latestdoi/199023313)

The paper can be accessed at [https://doi.org/10.1038/s41561-021-00725-0](https://doi.org/10.1038/s41561-021-00725-0)

Cite as:
>Fullerton, K.M., Schrenk, M.O., Yücel, M. et al. Effect of tectonic processes on biosphere–geosphere feedbacks across a convergent margin. Nat. Geosci. (2021). https://doi.org/10.1038/s41561-021-00725-0

This repository contains the raw data and the code used to analyze the interaction between the 16S rRNA microbial diversity and geochemistry across the hot spring of Costa Rica collected in the framework of the "Biology Meets Subductions" project, financed by the Deep Carbon Observatory and the Alfred P. Sloan Foundation.

Data analyses have been performed using the freely available open source software [Mothur](https://www.mothur.org/) and the R Statistical Software [https://www.R-project.org/](R development core team, 2010).

The file bms_mothur.txt contains the mothur commands, while the file Fullerton_et_al_BMS_16S_final_analysis.r contains the R code used for the analysis of the 16S and Fullerton_et_al_BMS_carbon_fiaxation_metagenome.r the code used for the analysis of the metagenome following [mi-faser annotation](https://services.bromberglab.org/mifaser/). The raw sequences for the 16S rRNA tag amplicon data can be obtained from the NCBI SRA archive under Bioproject number PRJNA579365, while the raw sequences for the metagenome are under Bioproject PRJNA627197. The folder 16S_rRMA_data contains the ASV count table, the taxonomy, tree and sample_data files for the 16S rRNA analysis. The folder mi-faser_metagenomes contains the _mi-faser_ annotations and the list of EC involved in carbon metabolism. The file SubductCR_bac_sample_table.csv contains all the other environmental data.

The preprint of the paper can be found on EarthXiv as preprint at [https://doi.org/10.31223/osf.io/gyr7n](https://doi.org/10.31223/osf.io/gyr7n).

For questions please contact me.

The code and data are released under the Creative Commons License Attribution 4.0 International (CC BY 4.0). Ream more about his license at https://creativecommons.org/licenses/by/4.0/.

![CC-BY](https://www.fosteropenscience.eu/learning/open-licensing/course/en/assets/b4467d3769dbd7a80e8d641361ff364b505d118d.png)
