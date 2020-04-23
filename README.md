# treeclimbR

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
3. AML-sim and BCR-XL-sim
4. Infant gut microbial data
5. Mouse miRNA data
6. Mouse cortex scRNAseq data

