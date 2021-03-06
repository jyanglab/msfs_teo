# Welcome to the `Evolution of DNA methylation ` project

This is a research repo for our project, entitled "**Evolutionary and functional genomics of DNA methylation in maize domestication and improvement**".

## Introduction
In this study, we leveraged whole-genome sequencing (WGS) and whole-genome bisulfite sequencing (WGBS) on populations of modern maize, landrace, and teosinte Zea mays ssp. parviglumis to investigate the adaptive and phenotypic consequences of methylation variation in maize.

## Architecture about this Repo
This project contains ~540 commits. A `largedata` directory was intentionally ignored by adding to `gitignore` because of the large size of the files within the folder. To guide the visitors having a better idea about the repo, here we briefly introduce the functions or sepecific purposes of the directory system. The layout of directories is based on the idea from [ProjectTemplate](http://projecttemplate.net/architecture.html). 

1. **cache**: Here we store intermediate data sets that are generated during a preprocessing step.
2. **data**: Here we store our raw data of small size. Data of large size, i.e. > 100M, store in a `largedata` folder that has been ignored using `gitignore`.
3. **reports**: Documentation codes (i.e. Rmd files) for generating the figures.
4. **graphs**: Graphs produced during the study.
5. **lib**: Some functions for our work.
6. **profilling**: Analysis scripts for the project. It contains some sub-directories.
7. **table**: Table produced during the study.

## Figures:

Rmd code to generate Figures in the paper.

1. **Figure 1**: Panels of this figure can be produced with code [`reports/Fig1_sfs.Rmd`](https://https://github.com/jyanglab/msfs_teo/tree/master/reports/Fig1_sfs.Rmd). And the [Figure1](https://github.com/jyanglab/msfs_teo/tree/master/graphs/fig1_mSFS.pdf) in a pdf version.
2. **Figure 2**: Panels of this figure can be produced with code [`reports/Fig2_Circos_Plot.Rmd`](https://github.com/jyanglab/msfs_teo/tree/master/reports/Fig2_Circos_Plot.Rmd). And the [Figure2](https://github.com/jyanglab/msfs_teo/tree/master/graphs/Fig2_circos.pdf) in a pdf version.
3. **Figure 3**: Panels of this figure can be produced with code [`reports/Fig3_CG_DMR_Function.Rmd`](https://github.com/jyanglab/msfs_teo/tree/master/reports/Fig3_CG_DMR_Function.Rmd). And the [Figure3](https://github.com/jyanglab/msfs_teo/tree/master/graphs/fig3_DMR_feature.pdf) in a pdf version.
2. **Figure 4**: Panels of this figure can be produced with code [`reports/Fig4_VCAP_Vgt1.Rmd`](https://github.com/jyanglab/msfs_teo/tree/master/reports/Fig4_VCAP_Vgt1.Rmd). And the [Figure4](https://github.com/jyanglab/msfs_teo/tree/master/graphs/fig4_vgt1.pdf) in a pdf version.
3. **Figure 5**: Panels of this figure can be produced with code [`reports/Fig5_tb1.Rmd`](https://github.com/jyanglab/msfs_teo/tree/master/reports/Fig5_tb1.Rmd). And the [Figure5](https://github.com/jyanglab/msfs_teo/tree/master/graphs/fig5_tb1.pdf) in a pdf version.

## Supplementary Tables:


## License
This repo is free and open source for research usage, licensed under [GPLv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).

