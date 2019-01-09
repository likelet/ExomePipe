#!/bin/sh
tumorPileup=$1
normalPileup=$2
samplename=$3
chrLen=$4
chrFile=$5
SNPfile=$6
captureRegions=$7

#unset LD_LIBRARY_PATH
echo "[general]" > ${samplename}_freec_config.txt
echo "chrLenFile = ${chrLen}" >> ${samplename}_freec_config.txt
echo "window = 0" >> ${samplename}_freec_config.txt
echo "ploidy = 2" >> ${samplename}_freec_config.txt
echo "outputDir = . " >> ${samplename}_freec_config.txt
echo "sex=XX" >> ${samplename}_freec_config.txt
echo "breakPointType=4" >> ${samplename}_freec_config.txt
echo "chrFiles = ${chrFile}" >> ${samplename}_freec_config.txt
echo "bedtools = bedtools" >> ${samplename}_freec_config.txt
echo "sambamba = sambamba" >> ${samplename}_freec_config.txt
echo "SambambaThreads = 23" >> ${samplename}_freec_config.txt
echo "samtools = samtools" >> ${samplename}_freec_config.txt
echo "maxThreads=23" >> ${samplename}_freec_config.txt
echo "breakPointThreshold=0.8" >> ${samplename}_freec_config.txt
echo "noisyData=TRUE" >> ${samplename}_freec_config.txt
echo "printNA=FALSE" >> ${samplename}_freec_config.txt
echo "readCountThreshold=10" >> ${samplename}_freec_config.txt
echo "[sample]" >> ${samplename}_freec_config.txt
echo "mateFile = ${tumorPileup}" >> ${samplename}_freec_config.txt
echo "inputFormat = pileup" >> ${samplename}_freec_config.txt
echo "mateOrientation = FR" >> ${samplename}_freec_config.txt
echo "[control]" >> ${samplename}_freec_config.txt
echo "mateFile = ${normalPileup}" >> ${samplename}_freec_config.txt
echo "inputFormat = pileup" >> ${samplename}_freec_config.txt
echo "mateOrientation = FR" >> ${samplename}_freec_config.txt
echo "[BAF]" >> ${samplename}_freec_config.txt
echo "SNPfile = ${SNPfile}" >> ${samplename}_freec_config.txt
echo "minimalCoveragePerPosition = 5" >> ${samplename}_freec_config.txt
echo "[target]" >> ${samplename}_freec_config.txt
echo "captureRegions = ${captureRegions}" >> ${samplename}_freec_config.txt