#! /bin/sh
#default:freec result of each sample is stored in FREEC_Dir/sampleID_freec/
#snv maf file of each patient is stored in snv_imputed_maf/

patientID=$1
FREEC_Dir=$2
snv_maf=$3
OutDir=$4

if [ -d $OutDir ]
then
	mkdir -p $OutDir
fi

if [ -d $OutDir/$patientID ]
then
	mkdir -p $OutDir/$patientID
fi 

if [ -d $OutDir/SNV ]
then
	mkdir -p $OutDir/SNV
fi

if [ -d $OutDir/SampleMaf ]
then
	mkdir -p $OutDir/SampleMaf
fi


cat ${snv_maf}/*.maf  > merged_imputed.maf
perl get_sample_maf.pl merged_imputed.maf $OutDir/SampleMaf

>$OutDir/${patientID}.tsv;
for maf in $OutDir/SampleMaf/*.maf
do
	sampleID=$(basename $maf|cut -d "." -f1)
	###get snv.txt
	perl generate_pyclone_snv_maf.pl $maf $sampleID
	SNV=$(ls $OutDir/SNV/$sampleID.snv.txt)

	##
	CNV=$(ls ${FREEC_Dir}/${sampleID}_freec/$sampleID.pileup_CNVs)
	intersectBed -a $SNV -b $CNV -wo|awk -F "\t" -v pre=$pre 'OFS="\t"{if($12>0){print $8,$6,$7,2,0,$12,pre}}' >> $OutDir/${patientID}.tsv	
done

Rscript sub_pycloneInput_FREEC.R $OutDir/$patientID $patientID
rm $OutDir/${patientID}.tsv
rm merged_imputed.maf