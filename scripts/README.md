# README

## Container
A singularity container was used for the `meffil` environment due to complicated installation / dependencies. `plink` and `king` are also available in the singularity container. 

The singularitry definition files [`meffil.def`](meffil.def) and [`tensorqtl.def`](tensorqtl.def) include container specifications. A pre-built images of our `meffil` environment is available at `quay.io` and can be  pulled using the following commands if you have singularity available on your system:

```bash
singularity pull oras://quay.io/datatecnica/meffil
```

We use a `tensorqtl` image pulled from Broad Institute (`francois4/tensorqtl:1.0.9`).

## Genotype preparation

[`prep-plink.sh`](prep-plink.sh) conducts LD-pruning of variants using `plink` and generates PCs with `king` to use as covariates for `meffil` and `tensorQTL`.

## EWAS

[`run-meffil.sh`](run-meffil.sh) is a wrapper for executing  [`meffil.R`](meffil.R) in the prebuilt singularity container for:
- importing `idat` files
- performing `mefffil` QC and sample validation
- calculating methylation beta estimates
- calculating methylation PCs
- normalizing methylation beta estimates
- subsetting to samples with genotype data
- running EWAS for iPSC and PBMC datasets

[`plot-EWAS.R`](plot-EWAS.R) generates QQ plots, manhattan plots, and methylation plots as a function of age (Figure 3).



## MethQTL analysis


[`build-covariates.R`](build-covariates.R) imports methylation values and generates properly-formatted phenotype and covariates files for tensorQTL.

[`plot-methQTL.R`](plot-methQTL.R) generates scatter plots in Figure 4B,C,D.

[`cell-specific-methQTL.R`](cell-specific-methQTL.R.R) conducts enrichment (GO, KEGG) analysis and plotting for PBMC or iPSC-specific methQTL (Figure 5A).

[`plot-genotype-specific-methylation.R`](plot-genotype-specific-methylation.R) generates methylation allele dosage plots for sites with shared or cell-type-specific effects of genotype on methylation (i.e. significant slope across genotypes) (Figure 5B,C,D).

[`plot-locus-tracks.R`](plot-locus-tracks.R) generates methylation + locus plots for genes of interest. `ADARB2`, `B3GNTL1`, `HLADPB2` and `SNTG2` (Figure 6).


## Clock analysis

[`get-clock-cpgs.R`](get-clock-cpgs.R) is a convenience script to extract lists of CpGs from R's `methylclockData` datasets in `ExperimentHub`. It generates a text file of CpGs for each clock analyzed. Used to generate plots by [`plot-methQTL-clocks.R`](scripts/check-EWAS-clocks.R) (Figure 2).


[`plot-methQTL-clocks.R`](plot-methQTL-clocks.R) generates stacked showing proportions of non-methQTL, iPSC- or PBMC-specific methQTL, or shared methQTL for lists of CpGs in various methylation clocks. (figures not included in final manuscript).

## Imprinting analysis

A list of imprinting genes, [`imprinting-genes.txt`](imprinting-genes.txt) was generated from `https://www.geneimprint.com/site/genes-by-species` by making the following modifications:
- Manually found and replaced unknown `?` characters
- Replaced comma with semicolon, e.g. `sed -i 's/, /;/g' imprinting-genes-before-correction.txt > imprinting-genes.txt`

[`imprinting-check.R`](imprinting-check.R) generates density plots for the subset of probes within imprinting genes compared to all non-imprinting genes (Figure 4F).


