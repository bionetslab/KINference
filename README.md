# KINference

This repository contains the source code of **KINference: Data-Driven Inference of Kinase Interaction Networks** alongside two example datasets.

<img src="README_files/imgs/workflow.png" style="display: block" />

## Installation instruction

1. Install with:
```R
devtools::install_github('bionetslab/KINference')
```

## How to run KINference

### 1. Input provided as intensity matrices
```R
library(Kinference)

run_KINference(
    x0.path = 'data-raw/example_data/example_x0.tsv',
    x1.path = 'data-raw/example_data/example_x1.tsv',
    apply.CORR = TRUE,
    m = 9,
    output.id = 'example_matrix_run', 
    output.path = 'example_results/'
)
``` 

### 2. Input provided in vector format
```R
library(Kinference)

run_KINference(
    f.path = 'data-raw/example_data/example_f.tsv',
    output.id = 'example_vector_run', 
    output.path = 'example_results/'
)
```
### 3. Only computing kinase enrichments
If you only want to compute kinase enrichments, set all filters to `FALSE`: `apply.DIFF = FALSE`, `apply.FS = FALSE`, `apply.PCST = FALSE`, and `apply.CORR = FALSE`.

## Input specification

KINference can be run with either 2 input matrices of intensity measurements of phosphorylation sites or 1 input vector of log2FC-transformed intensities. The inputs have to be provided as tab-separated files. The matrices and the vector must contain one column named `Protein`. This column contains the annotations of the UniprotIDs and the phosphorylated amino acid (AA) and its sequence position (Pos) in UniprotID_AAPos format. See the example datasets in `data-raw/example-data/` for reference. 

## Parameters of KINference

- `x0.path`: Path to the first intensity matrix file. Default is NA.
- `x1.path`: Path to the second intensity matrix file. Default is NA.
- `f.path`: Path to the intensity vector file. Default is NA.
- `output.path`: Directory where the results will be saved. Default is 'results'.
- `output.id`: Identifier for the output files. Default is 'key'.
- `species`: Species name, either 'Homo sapiens' or 'Mus musculus'. Default is 'Homo sapiens'.
- `translate_uniprots`: If set to false, the Uniprot IDs will not be translated to human Uniprot IDs. If set to FALSE, the PCST filter will be automatically disabled because the PCST filter needs a connected network and, therefore, needs all Uniprot IDs to be human Uniprot IDs because the provided kinase data is only for human kinases. Default is TRUE.
- `paired.samples`: Logical indicating if the samples are paired (see paper for difference in paired vs unpaired computation of log2fc). Default is TRUE.
- `apply.log2`: Logical indicating if log2 transformation should be applied to the data. Default is FALSE.
- `n`: Number of top kinases to infer. Default is 15.
- `alpha`: Parameter for baseline kinase inference. Default is 0.9.
- `apply.DIFF`: Logical indicating if the DIFF filter should be applied. Default is TRUE.
- `apply.FS`: Logical indicating if the FS filter should be applied. Default is TRUE.
- `apply.CORR`: Logical indicating if the CORR filter should be applied. Default is FALSE.
- `apply.PCST`: Logical indicating if the PCST filter should be applied. Default is TRUE.
- `beta`: Parameter for node filter computation. Default is 0.4.
- `gamma`: Parameter for kinase enrichment computation. Default is 1.0.
- `delta`: Parameter for edge filter computation. Default is 0.8.
- `epsilon`: CORR significance threshold. Default is 0.05.
- `m`: CORR minimum sample threshold. Default is 10.
- `multiple_testing_correction`: Method for multiple testing correction for CORR filter. Default is 'BH'.
- `custom_serine_threonine_kinase_data.path`: Path to a custom serine/threonine kinase data file. Default is NULL.
- `custom_tyrosine_kinase_data.path`: Path to a custom tyrosine kinase data file. Default is NULL.