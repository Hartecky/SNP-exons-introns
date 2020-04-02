#!bin/bash

# Genome partition for every chromosome

# chromosome 1
# awk '{if ($1=="1") print $0}' arabidopsis_genome.vcf >> chrom1.vcf

# chromosome 2
# awk '{if ($1=="2") print $0}' arabidopsis_genome.vcf >> chrom2.vcf

# chromosome 3
# awk '{if ($1=="3") print $0}' arabidopsis_genome.vcf >> chrom3.vcf

# chromosome 4
# awk '{if ($1=="4") print $0}' arabidopsis_genome.vcf >> chrom4.vcf

# chromosome 5
# awk '{if ($1=="5") print $0}' arabidopsis_genome.vcf >> chrom5.vcf


# Specification of introns and exons for each of the chromosomes

# Get Exons
# awk '/mRNA/ && !/intron_variant/ {print $0}' chrom1.vcf >> extracted_exons.txt

# Get Introns
# awk '/intron_variant/ {print $0}' chrom1.vcf >> extracted_introns.txt

# Get coordinates of all exons from all chromosomes
# awk '{if ($1=="1" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom_1_exonscoords.txt
# awk '{if ($1=="2" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom_2_exonscoords.txt
# awk '{if ($1=="3" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom_3_exonscoords.txt
# awk '{if ($1=="4" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom_4_exonscoords.txt
# awk '{if ($1=="5" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom_5_exonscoords.txt


# Get SNPs positions
# awk '{print $2}' extracted_exons.txt >> ex_snp_positions.txt

# Merge the files with SNP positions, and exons coordinates
# Format : 
# SNP | EXON START | EXON END
# paste ex_snp_positions.txt chrom_1_exonscoords.txt | cut -f 1,2,3 > main_exons_data.txt


# Final file with coordinates of snps, exon start and exon end\

# coords='main_snp_exons.txt'

# snp=( $( cat $coords | awk '{print $1}') )
# exon_start=( $( cat $coords | awk '{print $2}') )
# exon_end=( $( cat $coords | awk '{print $3}') )


#Working method for counting SNPs in every exon
awk ' {print $1}  $3!="" {print $2"S"; print $3"E"} ' main_exons_data.txt | sort -n | 
awk ' !/E|S/ {count++; next} /S/ {count=0; next} /E/ {print "Exon " line++": "count}' >> counted_exons_snp.txt