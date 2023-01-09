################################################################################
##                                                                            ##
## DNAm age/lifespan prediction                                               ##
## -------------------------------------------------------------------------- ##
## input data: Illumina DNA methylation (EPIC array ) data, .idat files       ##
## algorithm:  Horvath/Levine DNAm SkinBlood age and PhenoAge algorithm       ##
##             wateRmelon implementation Gorrie-Stone, Schalkwyk, El Khoury   ##
## output data: pred-age-result.txt                                           ##
## -------------------------------------------------------------------------- ##
## version: 2023-01-08-2045                                                   ##
## Copyright (C) Roland Zauner 2023                                           ##
##                                                                            ##
################################################################################

library(magrittr)

predict.age <- function(snp.dat.file="~/EthnoPredict/predict-age-ethno/input") {
  t.start <- Sys.time()
  cat("\ninitiating R script for DNA methylation (DNAm) data based age prediction ...")
  cat("\n-------------------------------------------------------------------------------------------------------\n")
  cat("Copyright (C) 2023 Roland Zauner\n
    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software for academic demonstration purpose, and you are welcome to redistribute it
    under conditions of the GNU GENERAL PUBLIC LICENSE (see <https://www.gnu.org/licenses/>) 
    which also applies for underlying software packages used by this script, eg statistical software R.
    The GNU General Public License does not permit incorporating your program
    into proprietary programs.\n
    This program utilises algorithms refered to as Horvath age clocks implemented in wateRmelon R package.
    Non-academic use may result in royalties and may be subject to licensing conditions and or restrictions.
    It is the user's responsibility to verify any regulations in this regard and comply with respective agreements.
    For further details see README.md file on https://github.com/roland-zauner/predict-age-ethno.")
  cat("\n-------------------------------------------------------------------------------------------------------\n")
  cat("Loading DNA methylation raw data provided in green and red channel intensity data (.idat) files ...\n")

  ## load .idat files - there are always two files for red/green channel ending with _Red.idat and _Grn.idat
  rgSet <- minfi::read.metharray.exp(base=snp.dat.file)
  
  ## data background correction and dye bias correction following Triche et al 2013 (https://doi.org/10.1093/nar/gkt090)
  cat("Preprocessing DNA methylation data ...\n")
  mSet.noob <- minfi::preprocessNoob(rgSet,dyeMethod="single")
  
  ## perform age prediction using wateRmelon implementation: Tyler Gorrie-Stone, Leo Schalkwyk, Louis El Khoury
  cat("Calculating age predictions based on algorithms derived by Steve Horvath and colleagues,\n
      based on wateRmelon implementation of Gorrie-Stone and colleagues \n")
  library(wateRmelon)
  age.pred <- agep(minfi::getBeta(mSet.noob), method='all')
  
  age.pred <- age.pred %>% 
    dplyr::select(dplyr::contains(c("skinblood","phenoage"))) %>% 
    dplyr::select(!dplyr::contains("missing")) %>% round(digits=1) %>%
    tibble::rownames_to_column("SampleID") %>% 
    dplyr::rename(DNAm.SkinBloodAge=dplyr::contains("skinblood"),DNAm.PhenoAge=dplyr::contains("phenoage"))
  
  cat("\nPrediction results:")
  print(age.pred)
  
  cat("\n-------------------------------------------------------------------------------------------------------")
  cat("\nprediction results were stored in output/pred-age-result.txt")
  t.end <- Sys.time()
  cat("\nprocessing time:",t.end-t.start,"s")
  cat("\nfinished R script for DNA methylation (DNAm) based age predictions...")
  cat("\n========================================================================================================\n")
  return(age.pred)
}


## retrieve arguments from command line
args <- commandArgs(trailingOnly = TRUE)
arg <- args

## debug
#arg <- NULL
#arg[1] <- "~/EthnoPredict/predict-age-ethno/input/"


if(length(arg)==0){stop("Rscript aborted due to missing argument, pls provide path to .idat files (_Red.idat and _Grn.idat).")}

if(length(list.files(arg[1],pattern="*.idat"))<2){stop("Rscript aborted due to missing .idat files, pls provide path to EPIC raw data (_Red.idat and _Grn.idat) ...")}

res <- predict.age(snp.dat.file=arg[1])

write.table(res, file="/home/output/pred-age-result.txt",row.names=F)

