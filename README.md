
# treeclimbR_article

This repository provides toy examples to understand treeclimbR and codes to reproduce figures in treeclimbR article.

### Installation

* Data container: [TreeSummarizedExperiment](https://github.com/fionarhuang/TreeSummarizedExperiment)
* Algorithm: [treeclimbR](https://github.com/fionarhuang/treeclimbR)
* Data visualization: [TreeHeatmap](https://github.com/fionarhuang/TreeHeatmap). 

```
BiocManager::install("fionarhuang/TreeSummarizedExperiment")
BiocManager::install("fionarhuang/treeclimbR")
BiocManager::install("fionarhuang/TreeHeatmap")
```

### Get started with toy examples ([click here](https://fionarhuang.github.io/treeclimbR_toy_example/))

Below are results of one toy dataset. 

* To capture signal patterns on the tree, `treeclimbR` firstly propose multiple candidates that are generated by tuning a parameter `t`:

<p align="center"> 
<img src="https://github.com/fionarhuang/treeclimbR_toy_example/blob/master/output/signal_cands.gif">
</p>

Heatmap shows counts of entities (rows) in samples (columns) split by groups. Branches that include entities with signals are colored in orange.

* Nodes detected by from `treeclimbR` are compared to those detected by `BH` under FDR 0.05.
<p align="center"> 
<img src="https://github.com/fionarhuang/treeclimbR_toy_example/blob/master/output/signal_result.png">
</p>


### Reproduce figures in the manuscript of treeclimbR
1. Parametric synthetical microbial data

2. Non-parametric synthetical microbial data

3. AML-sim and BCR-XL-sim ([see here](https://github.com/fionarhuang/treeclimbR_article/tree/master/cytof))
 The semi-simulated data are download in the `cytof/data/` folder from `HDCytoData` using [Download.R](https://github.com/fionarhuang/treeclimbR_article/tree/master/cytof/data) 
   - AML-sim (DA folder)
      - install snakemake
      - Set directory to `DA/`
      - Update paths to input and output files specified in the configuration file ([config.yaml](https://github.com/fionarhuang/treeclimbR_article/blob/master/cytof/DA/config.yaml))
      - Specify paths to your R libraries in the [.Renviron](https://github.com/fionarhuang/treeclimbR_article/blob/master/cytof/DA/.Renviron). If there is less or more than 3 library paths, then files in `analysis/` folder that have the code below also need to be updated correspondingly.
      ```
      .libPaths(c(
            Sys.getenv('R_LIBS_1'), 
            Sys.getenv('R_LIBS_2'),
            Sys.getenv('R_LIBS_3')))
      ```            
      - dry run the pipeline using `snakemake -npr` 
      - run the pipeline using `snakemake --cores n` (n is the number of cores to be used)
   - BCR-XL-sim (DS folder). Similar to run AML-sim pipeline.
 [Figues 3]() is generated using [all_figure.R](https://github.com/fionarhuang/treeclimbR_article/blob/master/cytof/summary/all_figure.R)
   
4. Infant gut microbial data ([see here](https://htmlpreview.github.io/?https://github.com/fionarhuang/treeclimbR_article/blob/master/microbe/docs/index.html))

5. Mouse miRNA data ([see here](https://htmlpreview.github.io/?https://github.com/fionarhuang/treeclimbR_article/blob/master/miRNA/docs/index.html))
6. Mouse cortex scRNAseq data([see here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/fionarhuang/treeclimbR_article/master/LPS/docs/index.html))

