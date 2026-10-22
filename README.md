# One-window test

## Software

Tested with Python 3.8.10, R 4.4.3, coloc 5.2.3, susieR 0.14.2,
jsonlite 2.0.0 and SMR 1.3.1. The method versions and SuSiE
`maxit=2000` match the website documentation.

## Install

Install Python 3.8 or newer and R 4.4.3. In R, install the three pinned
packages:

```sh
Rscript -e 'install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_version("susieR", version="0.14.2", repos="https://cloud.r-project.org"); remotes::install_version("coloc", version="5.2.3", repos="https://cloud.r-project.org"); remotes::install_version("jsonlite", version="2.0.0", repos="https://cloud.r-project.org")'
```

Install the SMR 1.3.1 binary separately. Set `RSCRIPT` and `SMR_BIN` to
the executable paths if they are not on `PATH`. The Python runner uses only
the standard library.

## Run

From this directory:

```sh
RSCRIPT=/path/to/Rscript SMR_BIN=/path/to/smr python3 run.py
```

Read `results/results.tsv` and `results/QA.json`. The script will not
overwrite an existing output directory. For a separate run, use
`python3 run.py --output another_results` with the same environment variables.
