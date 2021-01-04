#!/usr/bin/env nextflow

log.info """\
Vsearch 16S  -  N F    part 1 
================================
region   : $params.region
reads    : $params.reads
results  : $params.outdir
"""

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */

Channel
    .fromFilePairs( params.reads+"/*_{R1,R2}.fastq" )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { read_pairs_ch } 


process cleanReads {
	input:
	tuple val(pair_id), path(reads) from read_pairs_ch

	output:
	set val(pair_id), file("${pair_id}-R1-paired.fastq"), file("${pair_id}-R2-paired.fastq") into clean_read_pairs
	
	"""
	java -jar $params.trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} ${pair_id}-R1-paired.fastq ${pair_id}-R1-unpaired.fastq ${pair_id}-R2-paired.fastq ${pair_id}-R2-unpaired.fastq LEADING:15 TRAILING:15 MINLEN:200
	"""
}

process stitching {
	input:
	tuple val(pair_id), file(R1_clean), file(R2_clean) from clean_read_pairs

	output:
	file("${pair_id}-merged.fasta") into merged_read_pairs

	"""	
	vsearch --fastq_mergepairs ${R1_clean} --reverse ${R2_clean} --fastaout ${pair_id}-merged.fasta --fastqout_notmerged_fwd ${pair_id}-not-merged.fastq --fastqout_notmerged_rev ${pair_id}-not-merged.fastq --eetabbedout ${pair_id}-mergestatistics.txt --fastq_minmergelen 250 --label_suffix ';sample='${pair_id}
	"""
}


merged_read_pairs.collectFile(name: params.outdir).println{file -> "results in output file: ${file}\n" }
