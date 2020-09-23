

// requrement 
// samtools

// run analysis at single tumorbam file 

tumorbamPath = file(params.tumorbamPath)
// the name of bamfile should be sampleID_sort_dedup_realigned_recal.bam. This required that bamfile must be deduplicated and base ralliberated
maffile = file(params.maf)


freecFolder = params.outdir +"/Result/FREEC_OUTPUT/"

// setting annotation file path 
referenceMap = defineReferenceMap()



// samtools needed 
process mutationImputation {
    publishDir path: {params.outdir +"/Result/MWES/Imputation/"}, mode: "copy"
    input:
            file tumorbamPath
            file maffile 
            set file(genomeFile),file(genomeIndex),file(genomeDict),file(intervals) from Channel.value([
                            referenceMap.genomeFile,
                            referenceMap.genomeIndex,
                            referenceMap.genomeDict,
                            referenceMap.intervals
                            ])

    output:
        file "*.input.maf" into Impute_MAF
        file "*.input_raw.maf" into Impute_MAF
    script:
    """
        perl ${baseDir}/Mutation_imputation.pl ${maffile} ${genomeFile} ${tumorbamPath} ./
    """
}

// analysis data from multiExome seqeuencing

if(params.runpyclone && param.runFreeC){
// specify TSV file
   

    process FREEC_2_pyClone_Input {
              input:
                    file maffile from Impute_MAF
              output:
                    set idPatient, file("${idPatient}.tsv") into MutTsvFileForPyclone
              script:
              idPatient=patientID
              """
              ln -s ${baseDir}/generate_pyclone_snv_maf.pl .
              ln -s ${baseDir}/sub_pycloneInput_FREEC.R .

              sh ${baseDir}/get_PycloneInput_FREEC.sh ${idPatient} ${freecFolder} ${maffile} 

              """
            }
            


    // MutTsvFileForPyclone represents multiple mutation tsv file from single patients or same region
    process runPyclone {
      input:    
          set idPatient, file(muteTSV) from MutTsvFileForPyclone

      output:
      
      script:
      """
      mkdir ${idPatient}
        for tsv in *.tsv
        do
            pre=${tsv%%.tsv}
            PyClone build_mutations_file --in_file $muteTSV --out_file ${idPatient}/pyclone_analysis/yaml/${pre}.yaml --prior total_copy_number
        done
      """
    }


}

def helpMessage() {
  // Display help message
  log.info "    Usage:"
  log.info "       nextflow run MultiExomeProcess.nf -c nextflow_M.config --profile <your profile> --maf <maffile> --tumorbamPath <path> "
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
  print_parameter("\tMAFfile:        ", params.maf)
  print_parameter("\tTumor Bam:      ", params.tumorbam)
  print_parameter("\tdbsnp:          ",referenceMap.dbsnp)
  print_parameter("\tgenome:         ",referenceMap.genomeFile)
  print_parameter("\tintervals:      ",referenceMap.intervals)
  println LikeletUtils.print_green("-------------------------------------------------------------")
  println LikeletUtils.print_green("#                       Run analysis                          ")
  println LikeletUtils.print_green("-------------------------------------------------------------")
  checkAnalysis("\trunFreeC:         ",params.runFreeC)
  checkAnalysis("\trunPyclone:       ",params.runpyclone)
  println LikeletUtils.print_green("-------------------------------------------------------------")
}

def exomeSeqMessage() {
  // Display Sarek message
    LikeletUtils.sysucc_ascii()
    log.info ''
    print LikeletUtils.print_yellow('==========================================================================')+"\n"
    print LikeletUtils.print_yellow('=========        SYSUCC MultiExome seq data PIPELINE      ==========')+"\n"
    print LikeletUtils.print_yellow('==========================================================================')+"\n"
    log.info ''
    log.info 'Usage: '
    log.info 'Nextflow run MultiExomeProcess.nf -c nextflow_M.config  --profile <your profile> --maf <maffile> --tumorbamPath <path>'
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
  return file(params.genomes[params.genome]."${item}")
}