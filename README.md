# ExomePipe
## dependencies
   * fastp
   * bwa 
   * gatk
   * mutect2
   * strelka2
   * sequenza
   * tabix
   
   
# Input file 

*  `sample input file ` TSV file for sample
> adjusted from https://github.com/SciLifeLab/Sarek/blob/master/docs/INPUT.md
# TSV file for sample

Input files for ExomeSeqPipe can be specified using a tsv file given to the `--sample` parameter. The tsv file is a Tab Separated Value file with columns: `subject gender status sample lane fastq1 fastq2` or `subject gender status sample bam bai`.
The content of these columns should be quite straight-forward:

- `subject` designate the subject, it should be the ID of the Patient, or if you don't have one, it could be the Normal ID Sample.
- `gender` is the gender of the Patient, (XX or XY)
- `status` is the status of the Patient, (0 for Normal or 1 for Tumor)
- `sample` designate the Sample, it should be the ID of the Sample (it is possible to have more than one tumor sample for each patient)
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the bam file
- `bai` is the index

All examples are given for a normal/tumor pair. If no tumors are listed in the TSV file, then the workflow will proceed as if it was a single normal sample instead of a normal/tumor pair.

* Example TSV file for a normal/tumor pair with FASTQ files

In this sample for the normal case there are 3 read groups, and 2 for the tumor. It is recommended to add the absolute path of the paired FASTQ files, but relative path should work also. Note, the delimiter is the tab (\t) character:
NOTE: assume each sample has only one libraray
```
G15511    XX    0    C09DFN    pathToFiles/C09DFACXX111207.1_1.fastq.gz    pathToFiles/C09DFACXX111207.1_2.fastq.gz
G15511    XX    1    D0ENMT    pathToFiles/D0ENMACXX111207.1_1.fastq.gz    pathToFiles/D0ENMACXX111207.1_2.fastq.gz
```
