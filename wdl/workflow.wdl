version 1.0

import "tasks.wdl" as tasks

workflow wambam {
    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Compute some quick bam QC for long reads (e.g. read length, alignment identity)"
    }

    parameter_meta {
        BAM_FILE: "BAM file from running minimap2 with the --eqx flag. Either provide this file or both a FASTQ_FILE and REFERENCE_FILE."
        FASTQ_FILE: "Reads in a gzipped FASTQ file. Either provide this file and REFERENCE_FILE, or a BAM_FILE."
        REFERENCE_FILE: "FASTA file for the reference genome. Can be gzipped. Either provide this file and FASTQ_FILE, or a BAM_FILE."
    }

    input {
        File? BAM_FILE
        File? FASTQ_FILE
        File? REFERENCE_FILE
    }

    if(!defined(BAM_FILE) && defined(FASTQ_FILE) && defined(REFERENCE_FILE)){
        call tasks.runMinimap2 {
            input:
            readsFile=FASTQ_FILE,
            referenceFile=REFERENCE_FILE
        }
    }

    File cur_bam_file = select_first([BAM_FILE, runMinimap2.bam])
    
    call tasks.runWambam {
        input: bamFile=cur_bam_file
    }

    call tasks.makeWambamGraphs {
        input:
        identityCsv=runWambam.identityDist,
        lengthCsv=runWambam.lengthDist
    }

    output {
        File identity_dist_csv = runWambam.identityDist
        File length_dist_csv = runWambam.lengthDist
        File graph_pdf = makeWambamGraphs.graphsPdf
        File summary_csv = makeWambamGraphs.summaryCsv
        File? bam = runMinimap2.bam
        File? bam_index = runMinimap2.bam_index
    }
}
