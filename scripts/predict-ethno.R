################################################################################
##                                                                            ##
## genotype-based Ethnicity prediction                                        ##
## -------------------------------------------------------------------------- ##
## input data: Illumina SNP array data                                        ##
## algorithm: 170 AISNP panel described by Kidd & Seldin et al 2019           ##
##            see https://www.nature.com/articles/s41598-019-55175-x          ##
##            is used to calculate a summary metric ratio based on            ##
##            the product of allele frequencies derived from                  ##
##            1000 genome data, see: ftp://ftp.1000genomes.ebi.ac.uk/         ##
## output data: pred-ethno-result.txt                                         ##
## -------------------------------------------------------------------------- ##
## version: 2023-01-08-1600                                                   ##
## Copyright (C) Roland Zauner 2023                                           ##
##                                                                            ##
################################################################################

library(magrittr)

predict.ancestry <- function(snp.dat.file, aisnp.1k.file) {
  t.start <- Sys.time()
  cat("\ninitiating R script for AISNP-based ethnicity prediction ...")
  cat("\n-------------------------------------------------------------------------------------------------------\n")
  cat("Copyright (C) 2023 Roland Zauner\n
    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software for academic demonstration purpose, and you are welcome to redistribute it
    under conditions of the GNU GENERAL PUBLIC LICENSE (see <https://www.gnu.org/licenses/>) 
    which also applies for underlying software packages used by this script, eg statistical software R.
    The GNU General Public License does not permit incorporating your program
    into proprietary programs.\n
    This program utilises allele frequencies derived from 1000 genome project data and applies
    the Kidd and Seldin AISNP panel for ethnicity prediction.
    Non-academic use may result in royalties and may be subject to licensing conditions and or restrictions.
    It is the user's responsibility to verify any regulations in this regard and comply with respective agreements.
    For further details see README.md file on https://github.com/roland-zauner/predict-age-ethno.")
  cat("\n-------------------------------------------------------------------------------------------------------\n")
  cat("Loading SNP data file & calculating probability scores...\n")
  library(magrittr)
  load(file=aisnp.1k.file)
  
  gen.dat <- readr::read_delim(snp.dat.file, 
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE,skip = 11,
                               col_names=c("Sample_ID","SNP_Name","Chr","Position","RsID",
                                           "Allele1_Plus","Allele2_Plus"),show_col_types = FALSE) %>% 
    dplyr::mutate(id=RsID) %>% dplyr::left_join(aisnp.1kg[,c(3,5,9:13)],by="id")
  
  cat("\n",sum(aisnp.1kg$id %in% gen.dat$RsID),"out of",nrow(aisnp.1kg),"AISNPs (Kidd-Seldin panel, 170 SNPs) were detected in genomic sample data.\n")
  
  gen.snp.dat <- gen.dat[gen.dat$RsID %in% aisnp.1kg$id,]
  gen.snp.dat <- gen.snp.dat[!duplicated(gen.snp.dat$RsID),]
  
  ## combine info to genotype based on 1kg info
  gen.snp.dat$genotype <- paste0(ifelse(gen.snp.dat$Allele1_Plus==gen.snp.dat$alt,"1","0"),"|",
                                 ifelse(gen.snp.dat$Allele2_Plus==gen.snp.dat$alt,"1","0"))
  
  ind.snps <- matrix(nrow=nrow(gen.snp.dat),ncol=5,dimnames= list(gen.snp.dat$RsID,
                                                                  c("s.EUR","s.AFR","s.AMR","s.EAS","s.SAS")))     
  
  for(snp.id in gen.snp.dat$RsID) {
    ##single snp analysis within one sample
    qEUR= gen.snp.dat$EUR_AF[gen.snp.dat$id==snp.id]; pEUR=1-qEUR
    qAFR= gen.snp.dat$AFR_AF[gen.snp.dat$id==snp.id]; pAFR=1-qAFR
    qAMR= gen.snp.dat$AMR_AF[gen.snp.dat$id==snp.id]; pAMR=1-qAMR
    qEAS= gen.snp.dat$EAS_AF[gen.snp.dat$id==snp.id]; pEAS=1-qEAS
    qSAS= gen.snp.dat$SAS_AF[gen.snp.dat$id==snp.id]; pSAS=1-qSAS
    
    switch (gen.snp.dat$genotype[gen.snp.dat$RsID %in% snp.id],
            "0|0"= {pat.pEUR=pEUR*pEUR; pat.pAFR=pAFR*pAFR; pat.pAMR=pAMR*pAMR; pat.pEAS=pEAS*pEAS; pat.pSAS=pSAS*pSAS},    
            "1|0"= {pat.pEUR=pEUR*qEUR; pat.pAFR=pAFR*qAFR; pat.pAMR=pAMR*qAMR; pat.pEAS=pEAS*qEAS; pat.pSAS=pSAS*qSAS},
            "0|1"= {pat.pEUR=qEUR*pEUR; pat.pAFR=qAFR*pAFR; pat.pAMR=qAMR*pAMR; pat.pEAS=qEAS*pEAS; pat.pSAS=qSAS*pSAS},
            "1|1"= {pat.pEUR=qEUR*qEUR; pat.pAFR=qAFR*qAFR; pat.pAMR=qAMR*qAMR; pat.pEAS=qEAS*qEAS; pat.pSAS=qSAS*qSAS}
    )
    ind.snps[snp.id,] <- c(pat.pEUR,pat.pAFR,pat.pAMR,pat.pEAS,pat.pSAS)
  }#end-for: single snp analysis within one patient
  
  ## avoid 0
  ind.snps[ind.snps<=0] <- 0.001
  ## combine across AISNPs and log transform for scores
  pat.ai.scores <- log10(apply(ind.snps,2,prod))
  pat.ai.ratios <- apply(ind.snps,2,prod)/sum(apply(ind.snps,2,prod))
  
  ## nominate ancestry
  nom.ancestry <- function(dat) {
    return(c("EUR","AFR","AMR","EAS","SAS")[which.max(dat)])
  }
  pred.ancest <- nom.ancestry(pat.ai.scores)
  
  cat("\nCalculated ancestry prediction scores:\n")
  print(pat.ai.scores)
  
  cat("\nBased on AISNPs the prediction with highest likelihood for this dataset is:",pred.ancest,"descendance.")
  
  cat("\n-------------------------------------------------------------------------------------------------------")
  cat("\nprediction results were stored in output/pred-ethno-result.txt")
  t.end <- Sys.time()
  cat("\nprocessing time:",t.end-t.start,"s")
  cat("\nfinished R script for AISNP-based ethnicity prediction...")
  cat("\n========================================================================================================\n")
  return(list(ratios=pat.ai.ratios,scores=pat.ai.scores,pred=pred.ancest))
}


## retrieve arguments from command line
args <- commandArgs(trailingOnly = TRUE)
arg <- args

## debug
#arg <- NULL
#arg[1] <- "~/EthnoPredict/predict-age-ethno/input/SNP-test-dataset.txt"
#arg[2] <- "~/EthnoPredict/predict-age-ethno/scripts/aisnp.1kg.RData"

if(length(arg)==0){stop("Rscript aborted due to missing arguments, pls provide genotype and 1000 genome allele frequency (aisnp.1kg.RData) file...")}

if(file.exists(arg[2])==FALSE){stop("Rscript aborted due to missing 1000 genome allele frequency file, pls provide aisnp.1kg.RData file ...")}

if(file.exists(arg[1])==FALSE){stop("Rscript aborted due to missing genotype file, pls provide Infinium SNP array txt file ...")}

res <- dplyr::as_tibble(as.data.frame(predict.ancestry(snp.dat.file=arg[1], aisnp.1k.file=arg[2])),rownames="ethno") %>% 
  dplyr::mutate(ethno=gsub("s.","",ethno))

write.table(res, file="/home/output/pred-ethno-result.txt",row.names=F)

