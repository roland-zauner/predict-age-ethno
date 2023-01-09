## generate docker image for age and ethno prediction
## source: https://github.com/roland-zauner/predict-age-ethno.git
## zaro 2023-01-06


## get R base from posit in Ubuntu linux image
## https://github.com/rstudio/r-docker
FROM rstudio/r-base:4.2-jammy

## provide required folders - these should be mapped to volumes in the container
RUN mkdir /home/scripts /home/input /home/output

## install linux libs required by R packages
RUN apt-get update && apt-get install -y libxml2 \
libxml2-dev \
libssl-dev \
libpng-dev \
&& rm -r /var/lib/apt/lists/*

## the following R libraries are required for R script predict-ethno.R
RUN R -q -e 'install.packages(c("magrittr","dplyr","readr","BiocManager"),repos = "http://cran.us.r-project.org")'

## the following R libraries are required in addition for R script predict-age.R
RUN R -q -e 'BiocManager::install(c("minfi","wateRmelon","IlluminaHumanMethylationEPICmanifest"))'

## R scripts and essential data
COPY ./scripts /home/scripts