#! /usr/bin/env nextflow

// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"

// requirement:
// - bwa
// - picard
// - samtools/sambamba
// - annovar
// - bcftools
// - freec 
// - delly
// - bedtools
// - pyclone 
// - R packages : data.table, reshape2

version="0.0.3"

//user options
if (params.help) {
    this.helpMessage()
    exit 1
}

//
outdir=params.outdir
inputdir=params.inputdir

// setting annotation file path 
referenceMap = defineReferenceMap()

// parsing tsvPath 
if (params.sample) tsvPath = params.sample

// set fastq or bam file for further analysis
fastqFiles = Channel.empty()
bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  fastqFiles = LikeletUtils.extractFastq(tsvFile)
} else{
    log.info LikeletUtils.print_red("no sample table specified, plz specified a sample table ") 
    exit 1
}
fastqFiles.into{FastqforMapping;FastqforCheckMate}

//start message 
startMessage()

if(fastqFiles.length()==0) {
    LikeletUtils.print_red("No sample loaded with yout sample tsv file, plz check the input")
    exit 1
}

/*
// Step 1 Alignment
*/
// alignment 
process bwa_aligment{
        tag { file_tag }

        input:
            tuple idPatient, gender, status, idSample, file(fastqFile1), file(fastqFile2) from FastqforMapping
            tuple file(genomeFile), file(bwaIndex) from Channel.value([
                                                            referenceMap.genomeFile,
                                                            referenceMap.bwaIndex
                                                            ])


        output:
            tuple idPatient, status, idSample, file("${idSample}.bam") into MappedBam, MappedBamForQC
        
        script:
            file_tag_new=idSample
            file_tag=file_tag_new

            readGroup = "@RG\\tID:niID\\tPU:noPU\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"

            """
            bwa mem -t ${task.cpus} ${genomeFile} -M -R \"${readGroup}\" ${fastqFile1} ${fastqFile2} | \
            samtools sort --threads ${task.cpus} -m 2G - > ${idSample}.bam
            """
        }


/*
// Step 2 BAM QC
*/
// bamQC
process BAMqcByQualiMap {
        tag {idPatient + "-" + idSample}
        publishDir  path: { params.outdir + "/Result/Qualimap" }, mode: 'move', overwrite: true

        input:
            tuple idPatient, status, idSample, file(bam) from MappedBamForQC

        output:
            file(idSample) into bamQCmappedReport

        when: params.runQualimap

        script:
        """
        qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${idSample} \
        -outformat HTML
        """
        }


/*
// Step 3 BAM markDuplicates
*/
// markDuplicates
process MarkDuplicates {
        tag {idPatient + "-" + idSample}

        input:
            tuple idPatient, status, idSample, file("${idSample}.bam") from MappedBam

        output:
            tuple idPatient, status, idSample, file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicateMarkedBams,mdBamToJoin
            tuple idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
            file ("${idSample}.bam.metrics") into markDuplicatesReport

        script:
        """
        gatk --java-options ${params.markdup_java_options} \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 50000 \
        --INPUT ${idSample}.bam \
        --METRICS_FILE ${idSample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${idSample}_${status}.md.bam
        """
}

/*
// Step 4 GATK base recalliberation
*/
process CreateRecalibrationTable {

        tag {idPatient + "-" + idSample}

        input:
            tuple idPatient, status, idSample, file(bam), file(bai) from duplicateMarkedBams // realignedBam
            tuple file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
            referenceMap.genomeFile,
            referenceMap.genomeIndex,
            referenceMap.genomeDict,
            referenceMap.dbsnp,
            referenceMap.dbsnpIndex,
            referenceMap.knownIndels,
            referenceMap.knownIndelsIndex,
            referenceMap.intervals,
            ])

        output:
            tuple idPatient, status, idSample, file("${idSample}.recal.table") into recalibrationTable
            tuple idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

        script:
        known = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
        """
        gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        --input ${bam} \
        --output ${idSample}.recal.table \
        -R ${genomeFile} \
        -L ${intervals} \
        --known-sites ${dbsnp} \
        ${known} \
        --verbosity INFO
        """
}

recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0,1,2])


process RecalibrateBam {
        tag {idPatient + "-" + idSample}

        publishDir path:{params.outdir + "/Result/recal_bam"}, mode: 'link'


        input:
            tuple idPatient, status, idSample, file(bam), file(bai), file(recalibrationReport) from recalibrationTable
            tuple file(genomeFile), file(genomeIndex), file(genomeDict), file(intervals) from Channel.value([
            referenceMap.genomeFile,
            referenceMap.genomeIndex,
            referenceMap.genomeDict,
            referenceMap.intervals,
            ])

        output:
            tuple idPatient, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam, recalibratedBamForStats
            tuple idPatient, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV
        script:
        """
        gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${genomeFile} \
        --input ${bam} \
        --output ${idSample}.recal.bam \
        -L ${intervals} \
        --create-output-bam-index true \
        --bqsr-recal-file ${recalibrationReport}
        """
}


//=======================================================================
// call tumor and normal analysis with mutect2 and other softwares
//=======================================================================


/*
Step 5  Mutect call somatic variations
*/
// separate recalibrateBams by status
bamsNormal = Channel.create()
bamsTumor = Channel.create()

recalibratedBam.choice(bamsNormal, bamsTumor) {it[1] == 0 ? 1 : 0}
// remove status

bamsNormal = bamsNormal.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }
bamsTumor = bamsTumor.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

bamsAll = bamsNormal.combine(bamsTumor, by : 0 )
    // Since idPatientNormal and idPatientTumor are the same
    // It's removed from bamsAll Channel (same for genderNormal)
    // /!\ It is assumed that every sample are from the same patient

(BamForMutect2,BamForMutect2Pon,BamforFreeC,BamforFACET,BamforMSIsensor,BamforGATK4CNV,BamforGATKallelicCounts,BamForDelly) = bamsAll.map {
    idPatient, idSampleTumor, bamTumor, baiTumor, idSampleNormal, bamNormal, baiNormal ->
    [idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}.into(8)

//=======================================================================
// start analysis for futher analysis 
//=======================================================================


/*
Mutect2
*/


if(params.runMutect2){
            if(params.ponfile){
                pon.path=params.ponfile
                PONfile=Channel.fromPath(pon.path,pon.path+".tbi")
            }else{

                BamForMutect2Pon.map{
                    idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor -> [idSampleNormal,bamNormal,baiNormal]
                }.unique().set{BamForMutect2PonUnique}
               process tumorModeForPON {

                input:
                tuple idSampleNormal, file(bamNormal), file(baiNormal) from BamForMutect2PonUnique
                            file(intervals) from Channel.value(referenceMap.intervals)
                            tuple file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
                                referenceMap.genomeFile,
                                referenceMap.genomeIndex,
                                referenceMap.genomeDict,
                                referenceMap.dbsnp,
                                referenceMap.dbsnpIndex
                            ])
                
                output:
                tuple idSampleNormal, file("${idSampleNormal}.vcf.gz") into NormalIndPonVCF

                script:
                """
                 gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                            Mutect2 \
                        -R ${genomeFile} \
                        -I ${bamNormal} \
                        -tumor ${idSampleNormal} \
                        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
                        -L ${intervals} \
                        -O ${idSampleNormal}.vcf.gz

                """
            }



             process createPONbyCollapesePON {
                input:
                    file vcffiles from NormalIndPonVCF.collect()
                output:
                    tuple file("PON.vcf.gz"),file("PON.vcf.gz.tbi") into PONfile
                script:
                """
                    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                    CreateSomaticPanelOfNormals \
                    -vcfs 3_HG00190.vcf.gz \
                    -vcfs 4_NA19771.vcf.gz \
                    -vcfs 5_HG02759.vcf.gz \
                    -O PON.vcf.gz
                
                """
            }
     
            }
            
            // Mutect 
            process RunMutect2 {
                    tag {idSampleTumor + "_vs_" + idSampleNormal}


                    errorStrategy { task.exitStatus == 143 ? 'retry' : 'terminate' }
                    maxRetries 3

                    input:
                        tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from BamForMutect2
                        file(intervals) from Channel.value(referenceMap.intervals)
                        tuple file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
                            referenceMap.genomeFile,
                            referenceMap.genomeIndex,
                            referenceMap.genomeDict,
                            referenceMap.dbsnp,
                            referenceMap.dbsnpIndex
                        ])
                        tuple file(ponfile),file(ponfiletbz) from PONfile

                    output:
                        tuple val("mutect2"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output
                        
                    script:

                    bamforReview = "${idSampleTumor}_${idSampleTumor}_for_review.bam"
                    """
                        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                            Mutect2 \
                            -R ${genomeFile}\
                            -I ${bamTumor}  -tumor ${idSampleTumor} \
                            -I ${bamNormal} -normal ${idSampleNormal} \
                            -L ${intervals} \
                            -pon PON.vcf.gz \
                            --af-of-alleles-not-in-resource 0.0000025 \
                            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
                            -O ${idSampleTumor}_vs_${idSampleNormal}.vcf
                    """
            }



            /*
            Step 6  Filter Mutect output 
            */
            // filtered VCF
            process filtered_Mutect_VCFs{
                    tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

                    publishDir path: {params.outdir +"/Result/Mutect2_filtered_VCF"}, mode: "copy"
                    input:
                        tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcf) from mutect2Output
                        tuple file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                            referenceMap.genomeFile,
                            referenceMap.genomeIndex,
                            referenceMap.genomeDict
                        ])

                    output:
                            // we have this funny *_* pattern to avoid copying the raw calls to publishdir
                        tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${outputFile}") into VcfHardFiltered_Annovar,VcfHardFiltered_VEP
                        tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.somatic.filter_mark.vcf.gz") into vcfMarkedFiltered

                    script:
                    outputFile = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.filtered.vcf"

                    """
                        #run filter Mutect
                        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                            FilterMutectCalls \
                            -V ${vcf} \
                            -O ${idSampleTumor}.somatic.filter_mark.vcf.gz

                        # hard filter 

                        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                        SelectVariants \
                        -R ${genomeFile} \
                        -V ${idSampleTumor}.somatic.filter_mark.vcf.gz \
                        -select "vc.isNotFiltered()" \
                        -O ${outputFile}

                    """
            }

            /*
            Step 7  Annotated filter variants by Annovar/Vep
            */

            if(param.annotater="annovar"){
                process Annotate_Mutect_VCF{
                        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

                        publishDir path: {params.outdir +"/Result/Annotaed_Filtered_VAFs"}, mode: "copy"

                        input:
                            tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcfFiltered) from VcfHardFiltered_Annovar
                            tuple file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                                referenceMap.genomeFile,
                                referenceMap.genomeIndex,
                                referenceMap.genomeDict
                            ])

                        output:
                                // we have this funny *_* pattern to avoid copying the raw calls to publishdir
                            tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${outName}.somatic.anno.*_multianno.txt") into annovarTXT,annovarTXTforSamplefile,AnnovarForSingleMaf
                            
                        script:
                        outName = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}"

                        """
                            table_annovar.pl ${vcfFiltered} ${params.annovarDBpath} -buildver ${params.annovarDB} \
                                    -out ${outName}.somatic.anno -remove \
                                    -otherinfo \
                                    -protocol ${params.annovarProtocol} \
                                    -operation ${params.annovarOperation} \
                                    -nastring . \
                                    -vcfinput

                        """
                }
                 /*
            Step 8  Collapse and convert annoted file into MAF file for futher analysis 
            */

            annovarTXTforSamplefile.map { variantCaller,idPatient, idSampleNormal, idSampleTumor, vcfFiltered ->
                "${vcfFiltered}\t${idSampleTumor}\n"
            }.collectFile(name: 'filelist.txt').set {annoFile_sample}

            annovarTXT2=annovarTXT.map{variantCaller, idPatient, idSampleNormal, idSampleTumor, annotatedTxt ->
                [annotatedTxt]
            }

            process Convert2MAF{
                    tag "multiple maf"

                    publishDir path: {params.outdir +"/Result/Mutect2_merged_MAF"}, mode: "copy"

                    input:
                        file annotatedTxt from annovarTXT2.collect()
                        file mapFile from annoFile_sample
                    output:
                        file "Mutect_Merged.MAF"  into Mutect_mergedMAF

                    script:
                    """
                        perl ${baseDir}/bin/Mutect_Annovar_to_MAF.pl filelist.txt > Mutect_Merged.MAF

                    """
            }

            process Convert2singleMAF{
                    tag "singleMaf"

                    input:
                        tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file(annovarOut) from AnnovarForSingleMaf
                    output:
                        tuple idPatient,idSampleTumor,file("${idSampleTumor}.maf") into SingleMafForPyclone

                    script:
                    """
                        echo -e "${annovarOut}\t${idSampleTumor}" > filelist.txt
                        perl ${baseDir}/bin/Mutect_Annovar_to_MAF.pl filelist.txt > ${idSampleTumor}.maf

                    """
            }

            }

            if(param.annotater=="vep"){

                process Annotate_Mutect_VCF_withVeP{
                        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

                        publishDir path: {params.outdir +"/Result/VEPannotation_result"}, mode: "copy"

                        input:
                            tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcfFiltered) from VcfHardFiltered_VEP
                            tuple file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                                referenceMap.genomeFile,
                                referenceMap.genomeIndex,
                                referenceMap.genomeDict
                            ])

                        output:
                                // we have this funny *_* pattern to avoid copying the raw calls to publishdir
                            tuple variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${idPatient}_vep.maf") into VepMAF
                            
                        script:
                        outName = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}"

                        """
                            
                            perl vcf2maf.pl --input-vcf  ${vcfFiltered}       \
                                            --output-maf ${idPatient}_vep.maf           \
                                            --filter-vcf 0              \
                                            --ref-fasta  ${genomeFile}      \
                                            --tumor-id   ${idSampleTumor}       \
                                            --normal-id  ${idSampleNormal}      \
                                            --vep-data   ${vepDB}    \
                                            --species    homo_sapiens   \
                                            --ncbi-build ${param.genomebuild}   
                                
                                    
                        """
                }


                 

                process MergeMaf {
                  input:
                    file maffiles from VepMAF.collect()
                  output:
                    file "merged.maf" into Mutect_mergedMAF
                  script:
                  """ 
                    cat *.maf | awk '$NR==1 || $0!~/header pattern/' > merged.maf
                  
                  """
                }
            }



           

            /*
            Step 9 MAF summary analysis with MAFtools 
                    sample size should be more than 5.
            */

            process MAFsummaryAnalysis{
                tag "Mutect"

                publishDir path: {params.outdir +"/Result/MAFtools_summary"}, mode: "copy"

                input: 
                    file MAFsfile  from Mutect_mergedMAF
                output: 
                    file "*" into MAFtools_Output

                when:
                params.runMAFsummary

                script:
                """
                    Rscript ${baseDir}/bin/MAFanalysis.R ${MAFsfile}
                """

            }
}




/*
Step 10 FreeC analysis pipe (optional)
*/

if(params.runFreeC){
            FreeCreferenceMap = LikeletUtils.defineFREECref()


            // Generate pileup files
            process Generate_Pileup_file_For_freeC{
                tag {idPatient}

                input:
                    tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from BamforFreeC
                    file bedfile from Channel.value([
                            referenceMap.intervals
                        ])

                output:
                    tuple idPatient, idSampleNormal, file("${idSampleNormal}.pileup.gz"), idSampleTumor, file("${idSampleTumor}.pileup.gz") into pileupFileForConfig

                script:
                
                """
                samtools mpileup -l ${bedfile} -q 20 ${bamNormal} | gzip >  ${idSampleNormal}.pileup.gz
                samtools mpileup -l ${bedfile} -q 20 ${bamTumor} | gzip >  ${idSampleTumor}.pileup.gz
                """
            }

            // Generate pileup files
            process run_FreeC{
                tag {idPatient}

                publishDir path: {params.outdir +"/Result/FREEC_OUTPUT/"+idPatient}, mode: "link"

                input:
                    tuple idPatient, idSampleNormal, file(pileNormal), idSampleTumor, file(pileTumor)from pileupFileForConfig
                    tuple file(chrLenFile), chrFile, file(snpfile),file(bedfile) from Channel.value([
                                                            FreeCreferenceMap.chrLenFile,
                                                            FreeCreferenceMap.chrFile,
                                                            FreeCreferenceMap.snpfile,
                                                            FreeCreferenceMap.bedfile
                                                        ])

                output:
                    file "*"
                    tuple idPatient, idSampleTumor, file("${pileTumor}_ratio.txt") into CnvFileForfreeC
                    tuple idPatient, idSampleTumor, file("*.pileup_CNVs") into CnvFileForPyclone
                script:
                
                """
                sh ${baseDir}/bin/generate_Freec_Config_WEX.sh ${pileTumor} ${pileNormal} ${idPatient} ${chrLenFile} ${chrFile} ${snpfile} ${bedfile}

                freec -conf ${idPatient}_freec_config.txt > ${idPatient}_freec.log 

                """
            }

            process FREEC_2_GISTIC_coversion{
                tag {idPatient}

                publishDir path: {params.outdir +"/Result/GISTIC_INPUT/"}, mode: "move"

                input: 
                    tuple idPatient, file(ratioFIle) from CnvFileForfreeC
                    file chrLenFile from Channel.value([FreeCreferenceMap.chrLenFile])
                output:
                    tuple idPatient,file("${idPatient}.gistic.seg") into GISTICsegFile
                
                script:
                
                """
                perl ${baseDir}/bin/freec_result_to_GISTIC_seg.pl ${ratioFIle} ${chrLenFile} ${idPatient} > ${idPatient}.gistic.seg
                """
            }

            

           
}

/*
Step 11 Run analysis by FACET(optional)
*/
if(params.runFACET){    
            FACETrefMAP = LikeletUtils.defineFACETref()

            process prepareFACETfile{
                tag ""

                input:
                        tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from BamforFACET
                        file snpfile from Channel.value([FACETrefMAP.facetVcf])
                output:
                        tuple idPatient, file("${idPatient}.rm.chrM.csv.gz") into preparedFACETfile
                script:
                        """
                        # get pileup files 
                        snp-pileup -g -q15 -Q20 -P100 -r25,0 ${snpfile} ${idPatient}.csv  ${bamNormal} ${bamTumor}  

                        # remove chrX fiel 
                        zcat ${idPatient}.csv | grep -v chrM | gzip - > ${idPatient}.rm.chrM.csv.gz


                        """
            }

            process runFACETanalysis{
                tag ""

                publishDir path: {params.outdir +"/Result/FACET/"+idPatient}, mode: "move"


                input:
                        tuple idPatient, file(csvFile) from preparedFACETfile

                output:
                        file "${idPatient}*"

                script:
                        """
                        # run FACET 
                        Rscript ${baseDir}/bin/FACETS.R ${csvFile} ${idPatient}
                        """
            }

}

/*
Step 12 MSI-sensor (optional)
*/
if(params.runMSIsensor){

            process scanMSIfromGenome{
                tag "run once"
                
                storeDir params.storedir


                input:
                    file genomeFile from Channel.value([referenceMap.genomeFile])
                output:
                    file "microsatellites.list" into MicrosatellitesListFile 

                script:

                    """
                    msisensor scan -d ${genomeFile} -o microsatellites.list
                    """


            }
            process runMSIsensor{
                tag {idPatient}

                publishDir path: {params.outdir +"/Result/MSIsensor/"}, mode: "move"

                input:
                    tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor)from BamforMSIsensor
                    file MSIlistfile from MicrosatellitesListFile 
                    file bedfile from Channel.value([referenceMap.intervals])

                output:
                    file "${idPatient}*"

                script:
                    """
                    msisensor msi -d ${MSIlistfile} -n ${bamNormal} -t ${bamTumor} -e ${bedfile} -o ${idPatient}
                    """ 
            }
}

/*
 Step 13 Run CNV analysis by GATK4 (optional)
*/
// make sure the relative package should be installed in your system if you are going to plot cnv result by GATK
//https://github.com/broadinstitute/gatk-protected/blob/master/scripts/install_R_packages.R
if(params.rungatk4CNV){

            process preprocessTargegetlist {
                tag ""
                input:
                    tuple file(genomeFile),file(genomeIndex),file(genomeDict),file(intervals) from Channel.value([
                        referenceMap.genomeFile,
                        referenceMap.genomeIndex,
                        referenceMap.genomeDict,
                        referenceMap.intervals
                        ])

                output:
                    file "targets_C.preprocessed.interval_list" into ProcessedInterval
                script:
                """
                    gatk --java-options "-Xmx${task.memory.toGiga()}g"  PreprocessIntervals \
                        -L ${intervals} \
                        -R ${genomeFile} \
                        --bin-length 0 \
                        --interval-merging-rule OVERLAPPING_ONLY \
                        -O targets_C.preprocessed.interval_list
                """
            }
            process GAKT4_CollectReadCounts{
                tag {idPatient}
                    input:
                        tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor)from BamforGATK4CNV
                        file Pintervals from ProcessedInterval

                        output:

                        tuple idPatient,file("${idPatient}.tumor.counts.hdf5") into tumorCovFiles
                        file "${idPatient}.normal.counts.hdf5"  into normalCovFiles

                    script:
                            """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectReadCounts  \
                                      -I ${bamTumor} \
                                      -L ${Pintervals} \
                                      --interval-merging-rule OVERLAPPING_ONLY \
                                      -O ${idPatient}.tumor.counts.hdf5
                            gatk --java-options "-Xmx${task.memory.toGiga()}g" CollectReadCounts  \
                                    -I ${bamNormal} \
                                    -L ${Pintervals} \
                                    --interval-merging-rule OVERLAPPING_ONLY \
                                    -O ${idPatient}.normal.counts.hdf5
                            """  
                
            }
            
            (NormalCovFilesForfile, NormalCovFilesForMerge) = normalCovFiles.into(2)
            NormalCovFilesForfile.collectFile { file -> ['argument.file', "-I "+file.name + ' '] }.set { normalPONargmentFile }


             process GAKT4_CreateReadCountPanelOfNormals{
                tag "creat PON"
                input:
                    file normalCNVfilelist from normalPONargmentFile
                    file normallist from NormalCovFilesForMerge.collect()

                    output:
                        
                    file  "cnvponC.pon.hdf5"  into NormalCNVpon
                    script:
                            """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g" CreateReadCountPanelOfNormals \
                                --minimum-interval-median-percentile 5 \
                                --arguments_file ${normalCNVfilelist} \
                                --output ./
                                -O cnvponC.pon.hdf5
                            """  
                
            }

            process GAKT4_DenoiseReadCounts{
                tag {idPatient}
                input:
                    file NormalCNVpon
                    tuple idPatient,file(tumorPcov) from tumorCovFiles

                output:
                    tuple idPatient,file("${idPatient}.denoisedCR.tsv") into TumorDenoiseReadCounts,TumorDenoiseReadCountsForCallSeg,TumorDenoiseReadCountsForPlot
                    tuple idPatient,file("${idPatient}.standardizedCR.tsv"),file("${idPatient}.denoisedCR.tsv") into DiagnosisCRforPlot
                script:
                            """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g"   DenoiseReadCounts \
                                -I ${tumorPcov}  \
                                --count-panel-of-normals ${NormalCNVpon} \
                                --standardized-copy-ratios ${idPatient}.standardizedCR.tsv \
                                --denoised-copy-ratios ${idPatient}.denoisedCR.tsv
                            """  
                
            }

            process GAKT4_PlotDenoisedCopyRatios{
                tag {idPatient}
                publishDir path: {params.outdir +"/Result/GATK4_CNV/plot"}, mode: "move"

                input:
                    tuple idPatient,file(standardizedCR),file(denoisedCR) from DiagnosisCRforPlot
                    file genomeDict from Channel.value([referenceMap.genomeDict])
                output:
    
                   file "*"

                script:
                            """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g"  PlotDenoisedCopyRatios \
                                    --standardized-copy-ratios ${standardizedCR}\
                                    --denoised-copy-ratios ${denoisedCR}\
                                    --sequence-dictionary ${genomeDict} \
                                    --minimum-contig-length 46709983 \
                                    --output sandbox/plots \
                                    --output-prefix ${idPatient}
                            """  
                
            }
            


            process GATK4_CollectAllelicCounts{
                tag {idPatient}
                input:
                    tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor)from BamforGATKallelicCounts
                    tuple file(genomeFile),file(genomeDict),file(genomeIndex),file(knownTargetSnp),file(knownTargetSnpIndex) from Channel.value([
                        referenceMap.genomeFile,
                        referenceMap.genomeDict,
                        referenceMap.genomeIndex,
                        referenceMap.knownTargetSnp,
                        referenceMap.knownTargetSnpIndex
                        ])
                output:
                    tuple idPatient,file("${idPatient}.tumor.allelicCounts.tsv"),file("${idPatient}.normal.allelicCounts.tsv") into AllelicCountTSV

                script:
                """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g"  CollectAllelicCounts \
                                -L  ${knownTargetSnp}\
                                -I  ${bamNormal} \
                                -R ${genomeFile} \
                                -O ${idPatient}.normal.allelicCounts.tsv
                            gatk --java-options "-Xmx${task.memory.toGiga()}g"  CollectAllelicCounts \
                                -L  ${knownTargetSnp}\
                                -I  ${bamTumor} \
                                -R ${genomeFile} \
                                -O ${idPatient}.tumor.allelicCounts.tsv
                """
            }

            AllelicCountTSV.join(TumorDenoiseReadCountsForCallSeg).set{AllelicCountTSV_TumorDenoiseReadCountsForCallSeg}

            process GAKT4_ModelSegments{
                tag {idPatient}

                input:
                    tuple idPatient, file(tumorAlliecCount),file(normalAlliecCount),file(tumorDenosiseCount) from AllelicCountTSV_TumorDenoiseReadCountsForCallSeg

                output:
                    tuple idPatient, file("${idPatient}.cr.seg") into TumorCRsegment
                    tuple idPatient, file("${idPatient}.hets.tsv"),file("${idPatient}.modelFinal.seg") into TumorSegmentFileForPlot
                script:
                """
                 gatk --java-options "-Xmx${task.memory.toGiga()}g"  ModelSegments \
                        --denoised-copy-ratios ${tumorDenosiseCount} \
                        --allelic-counts ${tumorAlliecCount} \
                        --normal-allelic-counts ${normalAlliecCount}\
                        --output ./ \
                        --output-prefix ${idPatient}
                """
            }
            process GAKT4_CallCopyRatioSegments{
                 tag {idPatient}
                 input:
                    tuple idPatient, file(tumorCRsegment) from TumorCRsegment

                output:
                    tuple idPatient, file("${idPatient}.called.seg") into TumorCallSegment

                script:
                """
                 gatk --java-options "-Xmx${task.memory.toGiga()}g"  CallCopyRatioSegments \
                      --input ${tumorCRsegment}\
                      --output ${idPatient}.called.seg
                """
            }


            TumorSegmentFileForPlot.join(TumorDenoiseReadCountsForPlot).set{ TumorSegment_TumorDenoiseReadCountsForPlot}

             process GAKT4_PlotModeledSegments{
                tag {idPatient}

                publishDir path: {params.outdir +"/Result/GATK4_CNV/plot"}, mode: "move"
                input:
                    tuple idPatient,file(tumorHetTsv),file(tumoFinalSeg),file(tumorDenosiseCount) from TumorSegment_TumorDenoiseReadCountsForPlot
                    file genomeDict from Channel.value([referenceMap.genomeDict])
                output:
    
                   file "*"

                script:
                            """
                            gatk --java-options "-Xmx${task.memory.toGiga()}g"  PlotModeledSegments \
                                --denoised-copy-ratios ${tumorDenosiseCount} \
                                --allelic-counts ${tumorHetTsv} \
                                --segments ${tumoFinalSeg}\
                                --sequence-dictionary ${genomeDict} \
                                --minimum-contig-length 46709983 \
                                --output ./ \
                                --output-prefix ${idPatient}
                            """  
                
            }
            
}

/*
 Step 14 Run SV analysis by Delly (optional)
*/

// plz make sure that the delly installed in your system 
// to conver bcf to vcf, we also need bcftools
if(params.runDelly){
    delly_exc = file(params.delly_exclute)
    process runDellyCallSV {

     
      input:
    
            tuple idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from BamForDelly
            file delly_exc  
            tuple file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                            referenceMap.genomeFile,
                            referenceMap.genomeIndex,
                            referenceMap.genomeDict
                        ])          
      output:
            tuple idPatient, file("${idPatient}.delly.bcf") into DellyOutBcf
      script:
            """
            delly call -x ${delly_exc} -o ${idPatient}.bcf -g ${genomeFile} ${bamTumor} ${bamNormal}
            """
    }
    process runDellyFilter {
        
      input:
        tuple idPatient, file(dellyBcffil) from DellyOutBcf
      output:
        tuple idPatient, file("${idPatient}.filter.vcf") into DellyfilterdBcf
      script:
            """
                delly filter -f somatic -o ${dellyBcffil} -s ${delly_exc} ${idPatient}.filter.vcf
            """
    }

    process ConvertBcf2vcf {
        publishDir path: {params.outdir +"/Result/Delly/"}, mode: "move"
        input:
            tuple idPatient, file(dellyfilterBcf) from DellyfilterdBcf
        output:
            file "${idPatient}.delly.sv.vcf"
        script:
        """
            bcftools view ${dellyfilterBcf} > ${idPatient}.delly.sv.vcf
        """
    }

}




// NGScheckmate
if(params.runNGScheckmate){

            NGScheckmateRefMap = defineNGScheckMateRef()

            (FastqforCheckMateFile,FastqforCheckMate)=FastqforCheckMate.into(2)
            FastqforCheckMateFile.map{
                idPatient, gender, status, idSample, fastqFile1, fastqFile2 ->
                "${fastqFile1.name}\t${fastqFile2.name}\t${idSample}\n"
            }.collectFile(name: "inputfile_forCheckmate.txt").set{InputFastqList}

            FastqforCheckMate.map{
                idPatient, gender, status, idSample, fastqFile1, fastqFile2 ->
                [fastqFile1,fastqFile2]
            }.set{FastqforCheckMate}


            process runNGScheckmateAnalysis{
                tag ""

                publishDir path: {params.outdir +"/Result/NGScheckMate/"}, mode: "move"

                input:
                    file fastqFile from FastqforCheckMate.collect()
                    file fastqFilelist from InputFastqList
                    file snpPTfile from Channel.value([NGScheckmateRefMap.NGScheckmatePTfile])
                output:

                    file "*"

                script:
               
                    """
                        ncm_fastq.py -l ${fastqFilelist} -pt ${snpPTfile} -p ${task.cpus} -O checkmate
                    """

            }
}



// Downstream analysis
// Copy number analysis
// Facet.R analysis for purity analysis
// Freec for WES
/*
process runADTex{
    tag {file_tag}

    input:
        tuple val(samplename),file(tumorBam),file(normalBam)  from tumorBam_for_ADTex
        file targetBed

    output:

        tuple val(file_tag_new),file("${file_tag_new}_ADTex_cnv") into ADTexOutFolder

    when: !params.runADTex
    shell:
        file_tag=samplename
        file_tag_new=file_tag

        '''
        mkdir tmp 
        python !{adtexpath} --normal !{normalBam} --tumor !{tumorBam} --bed !{targetBad} --out !{file_tag_new}_ADTex_cnv
        '''

        }



// strelka2 analysis
process runStrelka2{
    tag {file_tag}

    input:
        tuple val(samplename),file(tumorBam)  from tumorBam_for_strelka
        tuple val(normalfilename),file(normalBam) from normalBam_for_strelka
        file targetBed
        file genome_ref_exome

    output:

        tuple val(file_tag_new),file("${file_tag_new}_somatic.*") into strelka_vcf_gz

    when: !params.runStrelka2

    shell:
        file_tag=samplename
        file_tag_new=file_tag
        thread=${strelka_thread}


        """
    
            # run configuration 
            
            python !{strelkapath}/bin/configureStrelkaSomaticWorkflow.py \\
                --normalBam ${normalBam} \\
                --tumorBam ${tumorBam} \\
                --ref ${genome_ref_exome} \\
                --callRegions $targetBed \\
                --runDir ./ \\
                --exome
            python ./runWorkflow.py -m local -j !{thread}
            
            mv results/variants/somatic.snvs.vcf.gz ${file_tag_new}_somatic.snv.vcf.gz
            mv results/variants/somatic.indels.vcf.gz ${file_tag_new}_somatic.indels.vcf.gz 
            mv results/variants/somatic.snvs.vcf.gz.tbi ${file_tag_new}_somatic.snv.vcf.gz.tbi
            mv results/variants/somatic.indels.vcf.gz.tbi ${file_tag_new}_somatic.indels.vcf.gz.tbi
        """

        }
// gc base file for  run Sequenza
process prepareSequenza{
    tag "produce GC file"
    storeDir { params.outdir + "/RequiredFiles" }

    input:

        file genome_ref_exome

    output:

        file " hg19.gc50Base.txt.gz" into genome_gc_file
    when: !params.runSequenza

    shell:
        """
        sequenza−utils.py GC−windows −w 50 ${genome_ref_exome} | gzip > hg19.gc50Base.txt.gz
        """

        }

process runSequenza_python{
        tag {file_tag}

        input:

            tuple val(samplename),file(tumorbam) from tumorBam_for_sequenza
            tuple val(normalName),file(normalbam) from normalBam_for_sequenza
            file targetBed
            file genome_gc_file
            file genome_ref_exome

        output:
            tuple val(samplename),file("stage2.seqz.gz") into sequenza_bin_data

        when: !params.runSequenza
        shell:
            file_tag = samplename
            """
            #stage 1 
            sequenza-utils bam2seqz  -n ${normalBam} -t ${tumorbam}\
                    --fasta ${genome_ref_exome} -gc ${genome_gc_file} > stage1.seqz.gz
            #stage 2  binning 
            sequenza-utils seqz-binning -w 200 -s stage1.seqz.gz | gzip > stage2.seqz.gz
            """

        }
process runSequenza_R{
        tag {file_tag}

        input:

        tuple val(samplename),file(bin_data) from sequenza_bin_data

        output:

        when: !params.runSequenza
        shell:
        """
        
        """

        }

*/

// run EXCAVATOR2

if(params.runEXCAVATOR2){
    process EXCAVATOR2_Prepare_Target {
      storeDir params.storedir

      input:
        tuple file(genomeFile),file(genomeIndex),file(genomeDict),file(intervals) from Channel.value([
                        referenceMap.genomeFile,
                        referenceMap.genomeIndex,
                        referenceMap.genomeDict,
                        referenceMap.intervals
                        ])
      output:
        // not applied 
        file "MyTarget_w50000" into EXCAVATOR2_targetFile
      script:
      """
        perl ${params.excavatarHOME}/TargetPerla.pl ${params.excavatarHOME}/SourceTarget.txt ${intervals} MyTarget_w50000 50000 hg19
      """
    }
    
    // from the manual 
    // Before running EXCAVATORDataPrepare.pl you need to create a space delimited file with three
    // fields: The absolute path to the .bam file you want to analyse, the path to the main sample output
    // folder and the sample name. The sample name will be used as a prefix/suffix for output files.

    process EXCAVATOR2_Prepare_Data {
    

      input:
      
      output:
      
      script:
      """
        perl EXCAVATORDataPrepare.pl ExperimentalFilePrepare.w50000.txt --processors ${task.cpus}} --target MyTarget_w50000 --assembly hg19
      """
    }

    process EXCAVATORDataAnalysis {
      input:
      
      output:
      
      script:
      """
      perl EXCAVATORDataAnalysis.pl ExperimentalFileAnalysis.w50K.txt \
      --processors 6 --target MyTarget_w50K --assembly hg19 \
      --output /.../OutEXCAVATOR2/Results_MyProject_w50K --mode paired
      """
    }
}




//pipeline log
workflow.onComplete {

    log.info LikeletUtils.print_green("SYSUCC bioinformatics center Exome analysis Pipeline Complete !")

    //email information
    if(params.mail){
        recipient=params.mail
        def subject = 'My Exome execution'

        ['mail', '-s', subject, recipient].execute() <<
                """

    LncPipe execution summary
    ---------------------------
    Your command line: ${workflow.commandLine}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}

    """
    }


}
workflow.onError {
    println LikeletUtils.print_yellow("Oops... Pipeline execution stopped with the following message: ")+LikeletUtils.print_red(workflow.errorMessage)
}




/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def helpMessage() {
  // Display help message
  log.info "    Usage:"
  log.info "       nextflow run ExomeSeqPipe --sample <file.tsv> "
  log.info "    --sample <file.tsv>"
}

def minimalInformationMessage() {
  // Minimal information message
  println LikeletUtils.print_green("-------------------------------------------------------------")
  println LikeletUtils.print_green("                       Checking Parameters                   ")
  println LikeletUtils.print_green("-------------------------------------------------------------")
  print_parameter("\tCommand Line:   ", workflow.commandLine)
  print_parameter("\tProfile:        ", workflow.profile)
  print_parameter("\tProject Dir:    ", workflow.projectDir)
  print_parameter("\tLaunch Dir:     ", workflow.launchDir)
  print_parameter("\tWork Dir:       ", workflow.workDir)
  print_parameter("\tOut Dir:        ", params.outdir)
  print_parameter("\tTSV file:       ", tsvFile)
  print_parameter("\tdbsnp:          ",referenceMap.dbsnp)
  print_parameter("\tgenome:         ",referenceMap.genomeFile)
  print_parameter("\tintervals:      ",referenceMap.intervals)
  println LikeletUtils.print_green("-------------------------------------------------------------")
  println LikeletUtils.print_green("#                       Run analysis                          ")
  println LikeletUtils.print_green("-------------------------------------------------------------")
  checkAnalysis("\trunFacet:         ",params.runFacet)
  checkAnalysis("\trunADTex:         ",params.runADTex)
  checkAnalysis("\trunstrelka:       ",params.runstrelka)
  checkAnalysis("\trunSequenza:      ",params.runSequenza)
  checkAnalysis("\trunQualimap:      ",params.runQualimap)
  checkAnalysis("\trunFreeC:         ",params.runFreeC)
  checkAnalysis("\trunFACET:         ",params.runFACET)
  checkAnalysis("\trunMAFsummary:    ",params.runMAFsummary)
  checkAnalysis("\trunMSIsensor:     ",params.runMSIsensor)
  checkAnalysis("\trunNGScheckmate:  ",params.runNGScheckmate)
  checkAnalysis("\trun Mutect2:      ",params.runMutect2)
  checkAnalysis("\trun GATK4 based CNV analysis:      ",params.rungatk4CNV)
  checkAnalysis("\trun Delly:      ",params.runDelly)
  println LikeletUtils.print_green("-------------------------------------------------------------")
}

def exomeSeqMessage() {
  // Display Sarek message
    LikeletUtils.sysucc_ascii()
    log.info ''
    print LikeletUtils.print_yellow('==========================================================================')+"\n"
    print LikeletUtils.print_yellow('=========        SYSUCC Exome seq data processing PIPELINE      ==========')+"\n"
    print LikeletUtils.print_yellow('==========================================================================')+"\n"
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow run ExomeseqPipe.nf '
}

def startMessage() {
  // Display start message
  this.exomeSeqMessage()
  this.minimalInformationMessage()
}


def defineReferenceMap() {
  return [
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex'         : checkParamReturnFile("bwaIndex"),
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
    // common snps located in the target region
    'knownTargetSnp' : checkParamReturnFile("knownTargetSnp"),
    'knownTargetSnpIndex' : checkParamReturnFile("knownTargetSnpIndex")
  ]
}

def checkParamReturnFile(item) {
  return file(params."${item}")
}

def defineFREECref() {
  return [
    'chrLenFile'            : checkParamReturnFile("freec_chrLen"),
    'chrFile'               : checkParamReturnFile("freec_chrFile"),
    // snp file for calculating BAF
    'snpfile'               : checkParamReturnFile("freec_snpfile"),
    // bed 
    'bedfile'               : checkParamReturnFile("freec_bedfile")
  ]
}

 def defineFACETref() {
  return [
    'facetVcf'            : checkParamReturnFile("facet_vcf")
  ]
}

 def defineNGScheckMateRef() {
  return [
    'NGScheckmatePTfile'  : checkParamReturnFile("NGScheckmatePTfile")
  ]
}


// spread bam and intervals
 def generateIntervalsForVC(bams, intervals) {
    def (bamsNew, bamsForVC) = bams.into(2)
    def (intervalsNew, vcIntervals) = intervals.into(2)
    def bamsForVCNew = bamsForVC.combine(vcIntervals)
    return [bamsForVCNew, bamsNew, intervalsNew]
}
 def print_parameter(content, parameter){
    println LikeletUtils.print_cyan(content)+LikeletUtils.print_green(parameter)
}
 def checkAnalysis(software,param){
if(param) println LikeletUtils.print_yellow(software)+LikeletUtils.print_green(param)
}