

// run analysis at single tumorbam file 

tumorbamPath = file(params.tumorbamPath)
// the name of bamfile should be sampleID_sort_dedup_realigned_recal.bam. This required that bamfile must be deduplicated and base ralliberated
maffile = file(maffiles)

freecFolder = params.outdir +"/Result/FREEC_OUTPUT/"

// setting annotation file path 
referenceMap = defineReferenceMap()

patientID = params.patientID


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