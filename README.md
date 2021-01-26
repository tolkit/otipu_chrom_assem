# Scripts and files associated with the chromosomal assembly of _O. tipulae_ manuscript
---

Prerequisites can be installed via:
```
conda create -n otipu -c r -c conda-forge r-base r-tidyverse r-reshape2 r-scales r-gtools r-ggpubr r-cluster r-tsne
conda activate otipu
```

After cloning the repository, R scripts can executed from the master directory.
Figures can be found under `report/figure` while analysis results will be generated under `analysis`.



```
git clone https://github.com/tolkit/otipu_chrom_assem.git
cd otipu_chrom_assem
Rscript scripts/nigon_definition.R
```

To regenerate the circos plot, circos must have been installed.

```
cd otipu_chrom_assem/analyses/circos
circos -conf etc/Oscheius_tipulae_3.3_all_chromosomes_Pif.circos.conf
circos -conf etc/Oscheius_tipulae_3.3_telomeres.conf
```
