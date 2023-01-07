# predict-age-ethno

Purpose of this example code is to illustrate how ethnicity can be predicted based on AISNPs (Kidd Seldin Panel) and biological and phenological age can be predicted based on DNA methylation data (Horvath & colleague clocks).


Run script from command line:
sudo docker run --rm -v input:/home/input -v output:/home/output pred-ethno:latest Rscript /home/scripts/predict-ethno.R /home/input/SNP-test-dataset.txt /home/scripts/aisnp.1kg.RData
