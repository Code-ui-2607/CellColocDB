

CellColocDB connects disease-associated GWAS regions with molecular QTL evidence across tissues and cell contexts. It builds on two resources maintained by our group:

- **CAUSALdb2** provides curated GWAS summary statistics and GWAS fine-mapping results. Its data processing and expanded fine-mapping framework are described in the [CAUSALdb2 publication](https://doi.org/10.1093/nar/gkae1096). The [original CAUSALdb fine-mapping pipeline](https://github.com/mulinlab/CAUSALdb-finemapping-pip) is also available.
- **QTLbase2** provides curated molecular QTL data across bulk tissues and other molecular contexts. Its sources and processing procedures are described in the [QTLbase2 publication](https://doi.org/10.1093/nar/gkac1020).

The upstream GWAS and QTL processing is not repeated in this repository. CellColocDB defines disease-associated GWAS windows, integrates the corresponding GWAS and QTL data, and evaluates colocalisation using **coloc.abf** as the primary method, with **coloc-SuSiE** and **SMR/HEIDI** as complementary analyses where suitable inputs are available.

## Scope of the example

This repository provides a reproducible **one-window test**. It contains saved data, scripts and settings for three GWAS–QTL tasks from one previously defined window. It demonstrates:

1. Reconstruction and validation of matched GWAS–QTL inputs, including variant identity and allele alignment.
2. Computation of coloc.abf, coloc-SuSiE and SMR/HEIDI results with saved method settings and the supplied LD data where required.
3. Comparison with expected numerical results and output of quality-control records.

The supplied example uses GWAS dataset **F900126** and a window on **chromosome 20:44,240,196–45,240,196**, with three molecular-feature/QTL-context combinations recorded in `config.json`. 

## Software

Tested with Python 3.8.10, R 4.4.3, coloc 5.2.3, susieR 0.14.2, jsonlite 2.0.0 and SMR 1.3.1. The method versions and SuSiE `maxit=2000` match the website documentation.

## Install

Install Python 3.8 or newer and R 4.4.3. In R, install the pinned packages:

```sh
Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_version("susieR", version="0.14.2", repos="https://cloud.r-project.org"); remotes::install_version("coloc", version="5.2.3", repos="https://cloud.r-project.org"); remotes::install_version("jsonlite", version="2.0.0", repos="https://cloud.r-project.org")'
```

Install the SMR 1.3.1 binary separately. Set `RSCRIPT` and `SMR_BIN` to the executable paths if they are not on `PATH`. The Python runner uses only the standard library.

## Run

From the repository directory:

```sh
RSCRIPT=/path/to/Rscript SMR_BIN=/path/to/smr python3 run.py
```

Read `results/results.tsv` and `results/QA.json`. The script will not overwrite an existing output directory. For another run, use `python3 run.py --output another_results` with the same environment variables.

The runner checks the stored input identities and data properties, verifies the LD-file checksum, runs the methods and compares the resulting statistics with values in `config.json`. 

## References

### Data resources

1. Wang J, et al. **CAUSALdb2: an updated database for causal variants of complex traits.** *Nucleic Acids Research.* 2025;53(D1):D1295–D1301. https://doi.org/10.1093/nar/gkae1096
2. Huang D, et al. **QTLbase2: an enhanced catalog of human quantitative trait loci on extensive molecular phenotypes.** *Nucleic Acids Research.* 2023;51(D1):D1122–D1128. https://doi.org/10.1093/nar/gkac1020
3. Wang J, et al. **CAUSALdb: a database for disease/trait causal variants identified using summary statistics of genome-wide association studies.** *Nucleic Acids Research.* 2020;48(D1):D807–D816. https://doi.org/10.1093/nar/gkz1026. Original fine-mapping pipeline: https://github.com/mulinlab/CAUSALdb-finemapping-pip

### Statistical methods

4. Giambartolomei C, et al. **Bayesian test for colocalisation between pairs of genetic association studies using summary statistics.** *PLOS Genetics.* 2014;10:e1004383. https://doi.org/10.1371/journal.pgen.1004383 — coloc.abf.
5. Wallace C. **A more accurate method for colocalisation analysis allowing for multiple causal variants.** *PLOS Genetics.* 2021;17:e1009440. https://doi.org/10.1371/journal.pgen.1009440 — coloc-SuSiE.
6. Zhu Z, et al. **Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets.** *Nature Genetics.* 2016;48:481–487. https://doi.org/10.1038/ng.3538 — SMR/HEIDI.
7. Wang G, Sarkar A, Carbonetto P, Stephens M. **A simple new approach to variable selection in regression, with application to genetic fine mapping.** *Journal of the Royal Statistical Society: Series B.* 2020;82:1273–1300. https://doi.org/10.1111/rssb.12388 — SuSiE model.
