# Title     : TODO
# Objective : TODO
# Created by: Qi Zhao
# Created on: 2018/8/27

message = "#=========================================================================================
#        USAGE: Rscript Sequenza_analysis.R <stage2.seqz.gz> sampleid
#
#  DESCRIPTION: Do copynumber analaysis with Sequenza
#
#       INPUTS: <stage2.seqz.gz> : result from seqeunza-utils
#               <Sample_info.txt> :
#
# REQUIREMENTS: [R Packages] : sequenza
#       AUTHOR: Qi Zhao, zhaoqi@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.0
#      CREATED: 2018/8/27
#     REVISION: ---
#================================================================================================"

# arg[1]: stage2.seqz.gz  arg[2]: sampleID
library("sequenza")
args = commandArgs(T)
suppressMessages(library("sequenza"))
data.file <- args[1]
seqzdata <- sequenza.extract(data.file)
CP.example <- sequenza.fit(seqzdata)
sequenza.results(sequenza.extract = seqzdata, cp.table = CP.example, sample.id = args[2], out.dir="./")