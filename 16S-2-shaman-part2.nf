#!/usr/bin/env nextflow


log.info """\
Vsearch 16S  -  N F    part 2 
================================
region   : $params.region
reads    : $params.reads
results  : $params.outdir
"""

merged_fastas_file = Channel.fromPath( params.outdir )
silva_db_file = Channel.fromPath( params.dbpath )

process formatMergedFasta {
	
	input:
	file merged_fastas from merged_fastas_file

	output:
	file("${merged_fastas}") into formated_merged_fastas 

	"""
	sed "s: ::g" -i ${merged_fastas}
	"""
}

formated_merged_fastas.into { formated_merged_fastas_to_derep; formated_merged_fastas_to_map }

process dereplication {

	input:
	file(formated_merged_fasta) from formated_merged_fastas_to_derep

	output:
	file("formated_merged_fasta_derep.fasta") into formated_merged_fastas_derep 

	"""
	vsearch --derep_fulllength ${formated_merged_fasta} --output formated_merged_fasta_derep.fasta -sizeout -strand both
	"""
}

process singletonRemoval {
	
	input:
	file(formated_merged_fastas_derep) from formated_merged_fastas_derep
	
	output:
	file("derep-singleton.fasta") into formated_merged_derep_singleton

	"""
	vsearch -sortbysize ${formated_merged_fastas_derep} --minsize 10 --output derep-singleton.fasta -sizein
	"""
}

process chimera_filtering {
	
	input:
	file(formated_merged_derep_singleton) from formated_merged_derep_singleton
	
	output:
	file("derep-singleton-nochimera.fasta") into formated_merged_fastas_derep_nochimera
	
	"""	
	vsearch --uchime_denovo ${formated_merged_derep_singleton} --nonchimeras derep-singleton-nochimera.fasta --sizein
	"""
}

process clustering {
	
	input:
	file(formated_merged_fastas_derep_nochimera) from formated_merged_fastas_derep_nochimera

	output:
	file("derep-singleton-nochimera-clustered.fasta") into derep_singleton_nochimera_clustered

	"""
	vsearch --cluster_size ${formated_merged_fastas_derep_nochimera} --centroids derep-singleton-nochimera-clustered.fasta --id 0.97 --relabel "OTU_" --sizein
	"""
}

derep_singleton_nochimera_clustered.into { derep_singleton_nochimera_clustered_to_annotate; derep_singleton_nochimera_clustered_to_map }

process otu_annotation{
	publishDir params.outdir,mode: 'move'

	input:
	file(derep_singleton_nochimera_clustered) from derep_singleton_nochimera_clustered_to_annotate
	file silva_db from silva_db_file
	output:
	file("annotated_otus.tsv") into annotation

	"""
	vsearch --usearch_global ${derep_singleton_nochimera_clustered} --db ${silva_db} --id 0.75 --blast6out annotation-silva.tsv

	python $params.taxscript -i annotation-silva.tsv -d ${silva_db} -o annotated_otus.tsv -dtype silva_ssu
	"""
}

process otu_mapping{
	publishDir params.outdir, mode: 'move'
	
	input:
	file(stitched_reads) from formated_merged_fastas_to_map
	file(otus) from	derep_singleton_nochimera_clustered_to_map
	
	output:
	file("otus-countmatrix.tsv") into counts_table
	"""
	vsearch --usearch_global ${stitched_reads} --db ${otus} --id 0.97 --strand both --otutabout otus-countmatrix.tsv
	"""
}
