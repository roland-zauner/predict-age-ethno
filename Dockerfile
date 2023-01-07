## generate docker image for age and ethno prediction
## source: https://github.com/roland-zauner/predict-age-ethno.git
## zaro 2023-01-06


## get R base from posit in Ubuntu linux image
## https://github.com/rstudio/r-docker
FROM rstudio/r-base:4.2-jammy

RUN mkdir /home/scripts /home/data

COPY ./predict-age-ethno /home/scripts