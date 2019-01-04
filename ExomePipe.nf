#! /usr/bin/env nextflow

// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"

// requirement:
// - bwa
// - picard
// - samtools/sambamba
// - annovar
// - bcftools



//user options
if (params.help) {
    this.helpMessage()

    exit 1
}

//
outdir=params.outDir
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

//start message 
startMessage()
/*
// Step 1 Alignment
*/
// alignment 
process bwa_aligment{
        tag { file_tag }

        input:
        set idPatient, gender, status, idSample, file(fastqFile1), file(fastqFile2) from fastqFiles
         set file(genomeFile), file(bwaIndex) from Channel.value([
            referenceMap.genomeFile,
            referenceMap.bwaIndex
            ])


        output:
        set idPatient, status, idSample, file("${idSample}.bam") into mappedBam, mappedBamForQC
        
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
            set idPatient, status, idSample, file(bam) from mappedBamForQC

        output:
            file(idSample) into bamQCmappedReport

        when: !params.runQualimap

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
            set idPatient, status, idSample, file("${idSample}.bam") from mappedBam

        output:
            set idPatient, status, idSample, file("${idSample}_${status}.md.bam"), file("${idSample}_${status}.md.bai") into duplicateMarkedBams,mdBamToJoin
            set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai") into markDuplicatesTSV
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
            set idPatient, status, idSample, file(bam), file(bai) from duplicateMarkedBams // realignedBam
            set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
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
            set idPatient, status, idSample, file("${idSample}.recal.table") into recalibrationTable
            set idPatient, status, idSample, val("${idSample}_${status}.md.bam"), val("${idSample}_${status}.md.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

        script:
        known = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
        """
        gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        --input ${bam} \
        --output ${idSample}.recal.table \
        --tmp-dir /tmp \
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

        publishDir outdir+"/recal_bam", mode: 'link'


        input:
            set idPatient, status, idSample, file(bam), file(bai), file(recalibrationReport) from recalibrationTable
            set file(genomeFile), file(genomeIndex), file(genomeDict), file(intervals) from Channel.value([
            referenceMap.genomeFile,
            referenceMap.genomeIndex,
            referenceMap.genomeDict,
            referenceMap.intervals,
            ])

        output:
            set idPatient, status, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam, recalibratedBamForStats
            set idPatient, status, idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

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
recalibratedBam = recalibratedBam.first()
recalibratedBam.choice(bamsNormal, bamsTumor) {it[1] == 0 ? 1 : 0}
//check files
bamsNormal = bamsNormal.ifEmpty{exit 1, "No normal sample defined, check TSV file: ${tsvFile}"}
bamsTumor = bamsTumor.ifEmpty{exit 1, "No tumor sample defined, check TSV file: ${tsvFile}"}
// remove status
bamsNormal = bamsNormal.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }
bamsTumor = bamsTumor.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }



//split input bed into multiple to help accelerate mutect2 calling 
process CreateIntervalBeds {
    tag {intervals.fileName}

    input:
        file(intervals) from Channel.value(referenceMap.intervals)

    output:
        file '*.bed' into bedIntervals mode flatten

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
        """
        awk -vFS="[:-]" '{
        name = sprintf("%s_%d-%d", \$1, \$2, \$3);
        printf("%s\\t%d\\t%d\\n", \$1, \$2, \$3) > name ".bed"
        }' ${intervals}
        """
}

// process interval files 
    bedIntervals = bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
        final fields = line.split('\t')
        if (fields.size() >= 5) duration += fields[4].toFloat()
        else {
            start = fields[1].toInteger()
            end = fields[2].toInteger()
            duration += (end - start) / params.nucleotidesPerSecond
        }
        }
        [duration, intervalFile]
    }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}


    (bamsNormalTemp, bamsNormal, bedIntervals) = generateIntervalsForVC(bamsNormal, bedIntervals)
    (bamsTumorTemp, bamsTumor, bedIntervals) = generateIntervalsForVC(bamsTumor, bedIntervals)
    //combine normal and tumor 
    bamsAll = bamsNormal.combine(bamsTumor)
    // Since idPatientNormal and idPatientTumor are the same
    // It's removed from bamsAll Channel (same for genderNormal)
    // /!\ It is assumed that every sample are from the same patient
    bamsAll = bamsAll.map {
    idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
    [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
    }
    // combine interval and bamsALL
    bamsTumorNormalIntervals = bamsAll.spread(bedIntervals)


// Mutect 
process RunMutect2 {
        tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

        input:
            set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsTumorNormalIntervals
            set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
            referenceMap.genomeFile,
            referenceMap.genomeIndex,
            referenceMap.genomeDict,
            referenceMap.dbsnp,
            referenceMap.dbsnpIndex
            ])

        output:
            set val("mutect2"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output
        
        script:
        """
            gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                Mutect2 \
                -R ${genomeFile}\
                -I ${bamTumor}  -tumor ${idSampleTumor} \
                -I ${bamNormal} -normal ${idSampleNormal} \
                -L ${intervalBed} \
                -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
        """
}
mutect2Output = mutect2Output.groupTuple(by:[0,1,2,3])

// merge VCF from split vcf
process ConcatVCF {
        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

        input:
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcFiles) from mutect2Output
            file(genomeIndex) from Channel.value(referenceMap.genomeIndex)

        output:
                // we have this funny *_* pattern to avoid copying the raw calls to publishdir
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated
                // TODO DRY with ConcatVCF

        script:
        outputFile = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.vcf"

        if(params.targetBED)		// targeted
                concatOptions = "-i ${genomeIndex} -c ${task.cpus} -o ${outputFile} -t ${params.targetBED}"
            else										// WGS
                concatOptions = "-i ${genomeIndex} -c ${task.cpus} -o ${outputFile} "

            """
            concatenateVCFs.sh ${concatOptions}
        """
}

/*
Step 6  Filter Mutect output 
*/
// filtered VCF
process filtered_Mutect_VCFs{
        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

        publishDir path: {params.outdir +"/Mutect2_filtered_VCF"}, mode: "copy"

        input:
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcf), file(tbi) from vcfConcatenated
            set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                referenceMap.genomeFile,
                referenceMap.genomeIndex,
                referenceMap.genomeDict
            ])

        output:
                // we have this funny *_* pattern to avoid copying the raw calls to publishdir
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${outputFile}") into vcfHardFiltered
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("${vcf}.somatic.filter_mark.vcf.gz") into vcfMarkedFiltered

        script:
        outputFile = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.filtered.vcf"

        """
            #run filter Mutect
            gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                FilterMutectCalls \
                -V ${vcf} \
                -O ${vcf}.somatic.filter_mark.vcf.gz

            # hard filter 

            gatk --java-options "-Xmx${task.memory.toGiga()}g" \
            SelectVariants \
            -R ${genomeFile} \
            -V ${vcf}.somatic.filter_mark.vcf.gz \
            -select "vc.isNotFiltered()" \
            -O ${outputFile}

        """
}

/*
Step 7  Annotated filter variants by Annovar 
*/

process Annotate_Mutect_VCF{
        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

        publishDir path: {params.outdir +"/Annotaed_Filtered_VAFs"}, mode: "copy"

        input:
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcfFiltered) from vcfHardFiltered
            set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
                referenceMap.genomeFile,
                referenceMap.genomeIndex,
                referenceMap.genomeDict
            ])

        output:
                // we have this funny *_* pattern to avoid copying the raw calls to publishdir
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("*.anno.txt") into annovarTXT
            
        script:
        outName = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}"

        """
            table_annovar.pl example/ex2.vcf humandb/ -buildver ${params.annovarDB} \
                    -out ${outName}.somatic.anno -remove \
                    -otherinfo \
                    -protocol refGene,cosmic70,avsnp147,ALL.sites.2015_08_edit,EAS.sites.2015_08_edit,esp6500siv2_all,exac03,ljb26_all,clinvar \
                    -operation g,f,f,f,f,f,f,f,f \
                    -nastring . \
                    -vcfinput

        """
}

/*
Step 8  Collapse and convert annoted file into MAF file for futher analysis 
*/

annovarTXT.map { variantCaller,idPatient, idSampleNormal, idSampleTumor, vcfFiltered ->
  "${vcfFiltered}\t${idSampleTumor}\n"
}.collectFile(name: 'filelist.txt').set { annoFile_sample }

process Convert2MAF{
        tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

        publishDir path: {params.outdir +"/Mutect2_merged_MAF"}, mode: "copy"

        input:
            set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(annotatedTxt) from annovarTXT.collect()
            file mapFile from annoFile_sample
        output:
            set variantCaller, file("${variantCaller}_Merged.MAF") into Mutect_mergedMAF
        script:

        outName = "${variantCaller}_Merged"

        """
            perl ${baseDir}/bin/Mutect_Annovar_to_MAF.pl > ${outName}.MAF

        """
}


// Downstream analysis
// Copy number analysis
// Facet.R analysis for purity analysis
// Freec for WES

process runADTex{
    tag {file_tag}

    input:
        set val(samplename),file(tumorBam),file(normalBam)  from tumorBam_for_ADTex
        file targetBed

    output:

        set val(file_tag_new),file("${file_tag_new}_ADTex_cnv") into ADTexOutFolder

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
        set val(samplename),file(tumorBam)  from tumorBam_for_strelka
        set val(normalfilename),file(normalBam) from normalBam_for_strelka
        file targetBed
        file genome_ref_exome

    output:

        set val(file_tag_new),file("${file_tag_new}_somatic.*") into strelka_vcf_gz

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

            set val(samplename),file(tumorbam) from tumorBam_for_sequenza
            set val(normalName),file(normalbam) from normalBam_for_sequenza
            file targetBed
            file genome_gc_file
            file genome_ref_exome

        output:
            set val(samplename),file("stage2.seqz.gz") into sequenza_bin_data

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

        set val(samplename),file(bin_data) from sequenza_bin_data

        output:

        when: !params.runSequenza
        shell:
        """
        
        """

        }






//pipeline log
workflow.onComplete {

    log.info "SYSUCC bioinformatics center Exome analysis Pipeline Complete"

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




/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def helpMessage() {
  // Display help message
  this.exomeSeqMessage()
  log.info "    Usage:"
  log.info "       nextflow run ExomeSeqPipe --sample <file.tsv> "
  log.info "    --sample <file.tsv>"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : " + tsvFile
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Target BED  : " + params.targetBED
  log.info "Containers"
  if (params.repository != "") log.info "  Repository   : " + params.repository
  if (params.containerPath != "") log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  dbsnp       :\n\t" + referenceMap.dbsnp
  log.info "\t" + referenceMap.dbsnpIndex
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  intervals   :\n\t" + referenceMap.intervals
}

def exomeSeqMessage() {
  // Display Sarek message
    LikeletUtils.sysucc_ascii()
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'SYSUCC Exome seq data processing PIPELINE '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow run ExomeseqPipe.nf '
}

def startMessage() {
  // Display start message
  LikeletUtils.sysucc_ascii()
  this.exomeSeqMessage()
  this.minimalInformationMessage()
}

def checkParamReturnFile(item) {
  return file(params."${item}")
}
//adjusted  from Saret 
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
  ]
}
// spread bam and intervals
def generateIntervalsForVC(bams, intervals) {
  def (bamsNew, bamsForVC) = bams.into(2)
  def (intervalsNew, vcIntervals) = intervals.into(2)
  def bamsForVCNew = bamsForVC.combine(vcIntervals)
  return [bamsForVCNew, bamsNew, intervalsNew]
}
