version 1.0

task runWambam {
    input {
        File bamFile
        Int memSizeGB = 10
    }

    Int diskSizeGB = round(2*size(bamFile, "GB")) + 50

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        wam -i ~{bamFile} -o wambam_results
	>>>

	output {
		File identityDist = "wambam_results/identity_distribution.csv"
		File lengthDist = "wambam_results/length_distribution.csv"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/wambam:latest"
        preemptible: 1
    }
}

task makeWambamGraphs {
    input {
        File identityCsv
        File lengthCsv
        Int memSizeGB = 4
        Int diskSizeGB = 50
    }

	command <<<
        set -eux -o pipefail

        Rscript /build/wambam/scripts/make_plots.R ~{identityCsv} ~{lengthCsv} wambam-graphs.pdf wambam-summary.csv
	>>>

	output {
        File graphsPdf = "wambam-graphs.pdf"
        File summaryCsv = "wambam-summary.csv"
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/wambam:latest"
        preemptible: 1
    }
}

task runMinimap2 {
    input {
        File? readsFile
        File? referenceFile
        String preset = "map-ont"
        Int kSize = 17
        Int indexSplitSizeGb = 8
        Boolean useMd = true
        Int memSizeGB = 128
        Int threadCount = 64
    }

    Int diskSizeGB = 10 * round(size(readsFile, "GB") + size(referenceFile, "GB")) + 50
    
	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        OUTPREF=$(basename ~{readsFile} | sed -E 's/.(fastq.gz|fq.gz|fasta|fasta.gz|bam)*$//')

        INPUT_READS=~{readsFile}
        if [ "${INPUT_READS: -3}" == "bam" ]
        then
            samtools fastq ~{readsFile} | minimap2 -x ~{preset} -K 3G ~{true="--MD" false="" useMd} -I ~{indexSplitSizeGb}g -a -c --eqx -t ~{threadCount} -k ~{kSize} ~{referenceFile} - | samtools view -hb - | samtools sort -@ ~{threadCount} -o $OUTPREF.bam -
        else
            minimap2 -x ~{preset} -K 3G ~{true="--MD" false="" useMd} -I ~{indexSplitSizeGb}g -a -c --eqx -t ~{threadCount} -k ~{kSize} ~{referenceFile} ~{readsFile} | samtools view -hb - | samtools sort -@ ~{threadCount} -o $OUTPREF.bam -
        fi
        samtools index -@ ~{threadCount} -b $OUTPREF.bam
	>>>

	output {
		File bam = glob("*.bam")[0]
        File bam_idx = glob("*.bam.bai")[0]
	}

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/wambam:latest"
        preemptible: 1
    }
}
