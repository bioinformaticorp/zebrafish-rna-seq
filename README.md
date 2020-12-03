# zebrafish-rna-seq

Dataset derived from sequencing of mRNA from Danio rerio embryos in two different developmental stages.

Sequencing was performed on the Illumina platform and generated 76bp paired-end sequence data using poly-(A)+ selected RNA. 
The data files are the following:
• 2cells_1.fastq and 2cells_2.fastq: these files are based on RNA-seq data of a 2-cell zebrafish embryo
• 6h_1.fastq and 6h_2.fastq: these files are based on RNA-seq data of zebrafish embryos 6h post fertilisation.
• Danio_rerio.Zv9.66.dna.fa: chromosome sequence in fasta format of Danio rerio.
• Danio_rerio.Zv9.66.gtf: annotation from Ensembl of Danio rerio chromosome sequence.
• hisat2.sh:
1. RNA-seq data aligned to the zebrafish genome using HiSat2. 
2. Perform transcriptome reconstruction using Cufflinks.
3. Compare the gene expression between two different conditions in order to identify differentially expressed genes using Cuffdiff.
