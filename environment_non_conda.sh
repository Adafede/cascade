#!/usr/bin/env bash

echo "Installing packages missing from conda"

R --slave -e "install.packages('readr', repo = \"https://cran.rstudio.com/\")"

R --slave -e "install.packages('MSnbase', repo = \"https://bioconductor.org/packages/release/bioc/\")"

R --slave -e "install.packages('nucleR', repo = \"https://bioconductor.org/packages/release/bioc/\")"
