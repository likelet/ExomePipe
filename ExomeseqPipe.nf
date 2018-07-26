#! /usr/bin/env nextflow

// usage : ./alignment.nf --input_folder input/ --cpu 8 --mem 32 --ref hg19.fasta --RG "PL:ILLUMINA"

// requirement:
// - bwa
// - picard
// - samtools/sambamba
// - annovar
// - bcftools


// Pipeline version
version = '0.0.1'

//user options
if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW Exome seq data processing PIPELINE v!{version}'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow ExomeseqPipe.nf '

    exit 1
}




//default values
params.help = null
params.input_folder = './'
params.out_folder = './'

// database file folder genome ref files

// database file folder genome ref files
params.genome_bwa_idex="/data/database/human/hg19/bwaIndex/genome"
params.genome_ref_exome="/data/database/human/hg19/bwaIndex/genome.fa"

params.gatkpath="/data/software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"
params.knownfile1="/data/database/human/hg19/forGATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
params.knownfile2="/data/database/human/hg19/forGATK/1000G_phase1.indels.hg19.sites.vcf"
params.knownfile3="/data/database/human/hg19/forGATK/dbsnp_138.hg19.vcf"
params.genomeref="/data/database/human/hg19/bwaIndex/genome.fa"
params.genomeindex="/data/database/human/hg19/bwaIndex/genome.fa.fai"
params.genomedict="/data/database/human/hg19/bwaIndex/genome.dict"

gatkpath=params.gatkpath
knownfile1=file(params.knownfile1)
knownfile2=file(params.knownfile2)
knownfile3=file(params.knownfile3)

genomeref=file(params.genomeref)
genomeindex=file(params.genomeindex)
genomedict=file(params.genomedict)



// software
//

// fastq file
params.fastq_ext = "*_{1,2}.fq.gz"

reads=params.fastq_ext
genome_bwa_idex = params.genome_bwa_idex
genome_ref_exome = file(params.genome_ref_exome)

bwa_index = Channel.fromPath("${params.genome_bwa_idex}*")


// read fastq
Channel.fromFilePairs(reads, size: 2)
                .ifEmpty {
            exit 1, print_red("Cannot find any reads matching: !{reads}\nNB: Path needs to be enclosed in quotes!\n")
        }
        .into {reads_for_bwa}

// alignment 
        process bwa_aligment{
            tag { file_tag }
            maxForks 5

            input:
            set val(samplename), file(pair) from reads_for_bwa
            file  bwa_id from bwa_index.collect()


            output:
            set val(file_tag_new), file("${file_tag_new}.bam") into mapped_bam
           
            shell:
            file_tag_new=samplename
            file_tag=file_tag_new
            bwa_threads = 10
	
            '''
                
 #               mkdir tmp

                bwa mem -t !{bwa_threads} genome.fa -M -R '@RG\tID:noID\tPL:ILLUMINA\tLB:noLB\tSM:'!{file_tag_new}'' !{pair[0]} !{pair[1]}  > !{file_tag_new}.sam
                sambamba  view -S -f bam -t !{bwa_threads} !{file_tag_new}.sam > !{file_tag_new}.bam                
                rm -rf !{file_tag_new}.sam
                '''
        }

process bamfile_process{
            tag { file_tag }
            maxForks 5
        input:
            set val(file_tag), file(bamfile) from mapped_bam


        output:
            set val(file_tag_new), file("${file_tag_new}_sorted_dedup.bam") into mapped_sorted_dedup_bam_from_sambamba,mapped_sorted_dedup_bam_from_sambamba2
        shell:
            file_tag_new = file_tag
            sambamba_threads = 6

            '''
                sambamba sort  -t !{sambamba_threads} -o !{file_tag_new}_sorted.bam !{file_tag_new}.bam
                sambamba markdup -t !{sambamba_threads}  !{file_tag_new}_sorted.bam !{file_tag_new}_sorted_dedup.bam
		        rm !{file_tag_new}_sorted.bam                

'''
        }


process RealignerTargetCreator {
	tag {file_tag}
    maxForks 5
 	input: 
 	file genomeref
	set val(samplename),file(bam) from mapped_sorted_dedup_bam_from_sambamba
	file knownfile1
	file knownfile2
	file genomeindex
	file genomedict
	
	output:
	set val(file_tag_new),file("${file_tag_new}_target_intervals.list") into intervalist


	shell:
	file_tag=samplename
	file_tag_new=file_tag
	'''
	sambamba index !{bam}
	java -jar !{gatkpath} -T RealignerTargetCreator -fixMisencodedQuals -R !{genomeref} -I !{bam} -known !{knownfile1} -known !{knownfile2} -o !{file_tag_new}_target_intervals.list -nt 6
	
	'''
}

//step2

process IndelRealigner {
        tag {file_tag}
        maxForks 5
        input:
            file genomeref
            set val(samplename),file(intervalfile) from intervalist
            set val(samplename),file(bamfile) from mapped_sorted_dedup_bam_from_sambamba2
            file knownfile1
            file knownfile2
            file genomeindex
            file genomedict

        output:
            set val(file_tag_new),file("${file_tag_new}_sort_dedup_realigned.bam") into realignmentBam,realignmentBam2
            set val(file_tag_new),file("${file_tag_new}_sort_dedup_realigned.bai") into realignmentBai,realignmentBai2
            shell:
            file_tag=samplename
            file_tag_new=file_tag
        '''
	sambamba index !{bamfile}
	java -jar !{gatkpath} -T IndelRealigner -R !{genomeref} -I !{bamfile} \
        -targetIntervals !{intervalfile} \
	-known !{knownfile1} -known !{knownfile2} -fixMisencodedQuals -o !{file_tag_new}_sort_dedup_realigned.bam 

        '''
}
//step3 
process BaseRecalibrator{
    tag {file_tag}
    maxForks 5
        input:
            file genomeref
            set val(samplename),file(bamfile) from realignmentBam
            set val(samplename),file(baifile) from realignmentBai
            file knownfile1
            file knownfile2
            file knownfile3
            file genomeindex
            file genomedict

        output:
            set val(file_tag_new),file("${file_tag_new}_sort_dedup_realigned_recal.grp") into GRPfile
	
	shell:
        file_tag=samplename
        file_tag_new=file_tag
	'''
	sambamba index !{bamfile}
        java -jar !{gatkpath} -T BaseRecalibrator -R !{genomeref} -I !{bamfile} \
        -knownSites !{knownfile3} -knownSites !{knownfile1} -knownSites !{knownfile2} -fixMisencodedQuals \
	-o !{file_tag_new}_sort_dedup_realigned_recal.grp  -nct 6 
#	java -jar $gatkpath -T BaseRecalibrator -R !{genomeref} -I !{bamfile} \
#	-BQSR !{file_tag_new}_sort_dedup_realigned_recal.grp \
#	-o !{file_tag_new}_sort_dedup_realigned_post_recal.grp  
        '''

}

process PrintReads{
        tag {file_tag}
        maxForks 5
        publishDir pattern: "*_recal.ba*", path: { params.out_folder + "/Result/RecalBam" }, mode: 'link', overwrite: true

    input:
        file genomeref
	    set val(samplename),file(GRPfile) from GRPfile
        set val(samplename),file(realignBam) from realignmentBam2
	    set val(samplename),file(baifile) from realignmentBai2
	    file knownfile1
        file knownfile2
        file knownfile3
	    file genomeindex
	    file genomedict
	
    output:
        set val(file_tag_new),file("${file_tag_new}_sort_dedup_realigned_recal.bam") into recalbam      

    shell:
        file_tag=samplename
        file_tag_new=file_tag
	'''
	sambamba index !{realignBam}
        java -jar !{gatkpath} -T PrintReads -R !{genomeref} -fixMisencodedQuals -I !{realignBam} -BQSR !{GRPfile} -o !{file_tag_new}_sort_dedup_realigned_recal.bam -nct 6
        '''

}


// copy number analysis
//bedtools 2.24 early
//DNAcopy R package
//ADTex
//python 2.7
process runADTex{
    tag {file_tag}

    input:
    set val(samplename),file(tumorBam)  from tumorBam_for_ADTex
    set val(normalfilename),file(normalBam) from normalBam_for_ADTex
    file targetBed

    output:

    et val(file_tag_new),file("${file_tag_new}_ADTex_cnv") into ADTexOutFolder


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


    shell:
    file_tag=samplename
    file_tag_new=file_tag
    thread=${strelka_thread}


    '''
 
        # run configuration 
        
        python !{strelkapath}/bin/configureStrelkaSomaticWorkflow.py \\
            --normalBam !{normalBam} \\
            --tumorBam !{tumorBam} \\
            --ref !{genome_ref_exome} \\
            --callRegions $targetBed \\
            --runDir ./ \\
            --exome
        python ./runWorkflow.py -m local -j !{thread}
        
        mv results/variants/somatic.snvs.vcf.gz ${file_tag_new}_somatic.snv.vcf.gz
        mv results/variants/somatic.indels.vcf.gz ${file_tag_new}_somatic.indels.vcf.gz 
        mv results/variants/somatic.snvs.vcf.gz.tbi ${file_tag_new}_somatic.snv.vcf.gz.tbi
        mv results/variants/somatic.indels.vcf.gz.tbi ${file_tag_new}_somatic.indels.vcf.gz.tbi
    '''

}

//run compelete

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
