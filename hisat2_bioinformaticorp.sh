#!/bin/bash

##########################################################################################
##########################################################################################

#PASO1


#INDEXAR un genoma

hisat2-build -p 4 Danio_rerio.Zv9.66.dna.fa Danio_rerio.Zv9.66.dna
	# -p for number of processors
	# first argument is the fasta file
	# second argument is the base name for indexed files


##########################################################################################
##########################################################################################

#PASO2

# Correr alineamiento HiSat2 para mapear cada set de lecturas

hisat2 -p 4 -x Danio_rerio.Zv9.66.dna -1 2cells_1_25000_reads.fq -2 2cells_2_25000_reads.fq -S 2cells.sam 

hisat2 -p 4 -x Danio_rerio.Zv9.66.dna -1 6h_1_25000_reads.fq -2 6h_2_25000_reads.fq -S 6h.sam


##########################################################################################
##########################################################################################

#PASO3

# convertir formato SAM a BAM 

samtools view --threads 4 -b -o 2cells.bam 2cells.sam

samtools view --threads 4 -b -o 6h.bam 6h.sam


##########################################################################################
##########################################################################################

#PASO4

# Ordenar BAM basado en coordenadas

samtools sort -m 2G -o 2cells_sorted.bam --threads 4 2cells.bam

samtools sort -m 2G -o 6h_sorted.bam --threads 4 6h.bam

#para versiones antiguas
#samtools sort 2cells.bam 2cells_sorted.bam


##########################################################################################
##########################################################################################

#PASO5

# Expresion de isoformas y ensamble de transcritos


cufflinks -o 2cells -G Danio_rerio.Zv9.66.gff -b Danio_rerio.Zv9.66.dna.fa -u 2cells_sorted.bam


cufflinks -o 6hours -G Danio_rerio.Zv9.66.gff -b Danio_rerio.Zv9.66.dna.fa -u 6h_sorted.bam


#Ordenar genes.fpkm_tracking por columna 10 en 2cells/ y 6hours/

sort -t$'\t' -r -g -k 10 genes.fpkm_tracking > genes.sorted.fpkm_tracking


#Tomar los mejores 100 genes ordenados por FPKM

head -100 genes.sorted.fpkm_tracking | cut -f 1 > top100FPKM.txt


##########################################################################################
##########################################################################################

#PASO6

# Expresion Diferencial

#Correr cuffdiff

cuffdiff -o cuffdiff_out -L ZV9_2cells,ZV9_6h -b Danio_rerio.Zv9.66.dna.fa -u Danio_rerio.Zv9.66.gff 2cells_sorted.bam 6h_sorted.bam


#Ver cuales son los genes con mayor expresion diferencial

sort -t$'\t' -g -k 13 cuffdiff_out/gene_exp.diff > cuffdiff_out/gene_exp_sorted.diff


#Cortar los 18 primeros

head -18 gene_exp_sorted.diff | cut -f 1 > top18genes.txt


##########################################################################################
##########################################################################################

#Para GRAFICAR es necesario emplear el lenguaje de programacion R a travez de un Script llamado cummeRbund.R

##########################################################################################
##########################################################################################






