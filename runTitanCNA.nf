


tumorbamPath = file(params.tumorbamPath)
maffile = file(maffiles)

freecFolder = params.outdir +"/Result/FREEC_OUTPUT/"

// setting annotation file path 
referenceMap = defineReferenceMap()

patientID = params.patientID

// from https://github.com/shahcompbio/hmmcopy_utils/tree/master/data
genome_bwfile=file("map_hg19.wig")

process gcCounter {

  input:
        set file(genomeFile),file(genomeDict) from Channel.value([
                        referenceMap.genomeFile,
                        referenceMap.genomeDict
                        ])
  output:
        file ("genome.gc.seg") into GenomeGCsegfile
  script:
  """
    gcCounter -w 100000 -s ALL ${genomeFile} > genome.gc.seg
  """
}


//mapCounter is a small program for calculating average mappability for non-overlapping windows of fixed width across all sequences (i.e. chromosomes) present in a BigWig file. It is built mainly on top of an independent subset of files obtained from the UCSC Genome Browser source code (i.e. kent library) made available by Jim Kent. Generating average mappability files in 1000 base windows (default) on hg18 took 160 seconds on a 3.06 GHz Intel Core 2 Duo with 8GB RAM.

process MapCounter {
  input:
        file genome_bwfile
  output:
        file ("genome.map.seg") into GenomeMapsegFile
  script:
  """
    ./mapCounter -w 100000 -s ALL genome_bwfile > hg18.map.seg
  """
}

process readCounter  {
  input:
   file bamfile from bamfilePath
  output:
    file("readcounts.seg") into RecountSegfile
  script:
  """
  readCounter -w 100 -s all ${bamfile} > readcounts.seg
  """
}



process titanCNA {
  input:
      
  output:
  
  script:
  """
  
  """
}


// functions --------------

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