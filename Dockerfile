## generate docker image for age and ethno prediction
## source: https://github.com/roland-zauner/predict-age-ethno.git
## zaro 2023-01-06


## get R base from posit in Ubuntu linux image
## https://github.com/rstudio/r-docker
FROM rstudio/r-base:4.2-jammy

## provide required folders - these should be mapped to volumes in the container
RUN mkdir /home/scripts /home/input /home/output

## the following R libraries are required
RUN R -q -e 'install.packages(c("magrittr","dplyr","readr"),repos = "http://cran.us.r-project.org")'

## R scripts and essential data
COPY ./scripts /home/scripts