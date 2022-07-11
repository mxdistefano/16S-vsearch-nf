# 16S-vsearch-nf
Vsearch based targeted metagenomics data analysis pipeline.

This pipeline was developed to provide an accurate and fast pipeline for targeted metagenomics based on 16S rRNA high throughput sequencing. It was designed using Vsearch toolkit and adapted to work with Nextflow, to cover the analysis from raw data (fastq files) to the: 
- Annotation table: OTUs annotated using an available database 
- Count matrix: abundancy of each OTU per sample.

![imagen](https://user-images.githubusercontent.com/29629513/178354401-2d3b50af-0573-4b8e-9317-3bc520aaf50c.png)

Dependecies:
- Python > 2.7
- Nextflow
- trim-o-matic 
- Vsearch
- Docker

