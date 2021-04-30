#!bin/bash

# PROJECT: 

# CHARACTERISTICS OF A SINGLE NUCLEOTIDE POLYMORPHISMS IN ARABIDOPSIS THALIANA

# ### GENOME PARTITION ---------------------------------------------------------------------------------------------
# # chromosome 1
# awk 'NR==1 || NR==23 {print $0} {if ($1=="1") print $0}' arabidopsis_genome.vcf >> chrom1.vcf

# # chromosome 2
# awk 'NR==1 || NR==23 {print $0} {if ($1=="2") print $0}' arabidopsis_genome.vcf >> chrom2.vcf

# # chromosome 3
# awk 'NR==1 || NR==23 {print $0} {if ($1=="3") print $0}' arabidopsis_genome.vcf >> chrom3.vcf

# # chromosome 4
# awk 'NR==1 || NR==23 {print $0} {if ($1=="4") print $0}' arabidopsis_genome.vcf >> chrom4.vcf

# # chromosome 5
# awk 'NR==1 || NR==23 {print $0} {if ($1=="5") print $0}' arabidopsis_genome.vcf >> chrom5.vcf



# ### EXTRACTING EXONS LENGTH FROM CHROMOSOMES -----------------------------------------------------
# chromosome 1
# awk 'BEGIN {print "exon_len"} /exon/ {if ($1=="1") print $7}' gene_features.data >> chrom1_len.txt

# # chromosome 2
# awk 'BEGIN {print "exon_len"} /exon/ {if ($1=="2") print $7}' gene_features.data >> chrom2_len.txt

# # chromosome 3
# awk 'BEGIN {print "exon_len"} /exon/ {if ($1=="3") print $7}' gene_features.data >> chrom3_len.txt

# # chromosome 4
# awk 'BEGIN {print "exon_len"} /exon/ {if ($1=="4") print $7}' gene_features.data >> chrom4_len.txt

# # chromosome 5
# awk 'BEGIN {print "exon_len"} /exon/ {if ($1=="5") print $7}' gene_features.data >> chrom5_len.txt



# ### EXTRACTING COORDINATES FOR BOTH CODING / NON-CODING SEQUENCES --------------------------------------------------

# # GENES COORDINATES
# awk '{if ($1=="1" && $4=="GENE") print $5 " " $6}' gene_features.data >> chrom1genecoords.txt
# awk '{if ($1=="2" && $4=="GENE") print $5 " " $6}' gene_features.data >> chrom2genecoords.txt
# awk '{if ($1=="3" && $4=="GENE") print $5 " " $6}' gene_features.data >> chrom3genecoords.txt
# awk '{if ($1=="4" && $4=="GENE") print $5 " " $6}' gene_features.data >> chrom4genecoords.txt
# awk '{if ($1=="5" && $4=="GENE") print $5 " " $6}' gene_features.data >> chrom5genecoords.txt

# # EXONS COORDINATES
# awk '{if ($1=="1" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom1exonscoords.txt
# awk '{if ($1=="2" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom2exonscoords.txt
# awk '{if ($1=="3" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom3exonscoords.txt
# awk '{if ($1=="4" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom4exonscoords.txt
# awk '{if ($1=="5" && $4=="exon") print $5 " " $6}' gene_features.data >> chrom5exonscoords.txt


# ### INFO COLUMN-----------------------------------------------------------------------------------------------------
# bioawk -F ";" -c vcf '{print $info}' chrom*.vcf >> info_data.vcf

# TSA - polymorphism
# awk -F ";" '{print $2}' info_data.vcf  >> info_TSA.txt

# # VE - raw variants
# awk -F "VE=" '{split($2,a," ");print FS a[1]}' info_data.vcf >> info_VE.txt


# ### UNIQ VARIANTS --------------------------------------------------------------------------------------------------

# TSA UNIQ
# awk -F "=" '{print $2}' info_TSA.txt | sort | uniq -c >> uniq_variants.txt

# VE UNIQ
# awk -F "=" '{print $2}' info_VE.txt | awk -F "|" '{print $1}' | sort | uniq -c >> uniq_VE.txt 

### GET SNPS ONLY --------------------------------------------------------------------------------------------------
awk '/TSA=SNV/ {print $0}' chrom*.vcf >> chromSNP.vcf

### EXTRACT EXONS AND INTRONS FROM A CHROM_SNP DATA ----------------------------------------------------------------

# EXONS
awk '/VE=missense_variant/ || /VE=synonymous_variant/ || /VE=stop_lost/ || /VE=stop_retained_variant/ || /VE=mature_miRNA_variant/	|| /VE=protein_altering/ {print $0}' chromSNP.vcf >> extracted_exons.vcf

# INTRONS
awk '/VE=3_prime_UTR_variant/ || /VE=5_prime_UTR_variant/ || /VE=intron_variant/ || /VE=splice_acceptor_variant/ || /VE=splice_donor_variant/ || /VE=splice_region_variant/	|| /VE=start_lost/ {print $0}' chromSNP.vcf >> extracted_introns.vcf

### SNP POSITIONS AND VARIANTS ID ----------------------------------------------------------------------------------

# EXONS
awk '{print $2, $3}' extracted_exons.vcf >> ex_snp_positions.txt

# INTRONS
awk '{print $2, $3}' extracted_introns.vcf >> in_snp_positions.txt

### MERGING SNP POSITIONS WITH COORDINATES -------------------------------------------------------------------------

# FORMAT: 
# SNP POSITION | EXON START (S) | EXON END (E)

paste ex_snp_positions.txt chrom*exonscoords.txt | cut -f 1,2,3 > main_exons_data.txt
paste in_snp_positions.txt chrom*intronscoords.txt | cut -f 1,2,3 > main_introns_data.txt


### SNP COUNTER ---------------------------------------------------------------------------------------------------

# EXONS
awk ' {print $1}  $3!="" {print $2"S"; print $3"E"} ' main_exons_data.txt | sort -n |  
awk 'BEGIN {print "SNP"} !/E|S/ {count++; next} /S/ {count=0; next} /E/ {print count}' >> counted_exons_snp.txt


# INTRONS
awk ' {print $1}  $3!="" {print $2"S"; print $3"E"} ' main_introns_data.txt | sort -n | 
awk 'BEGIN {print "SNP"} !/E|S/ {count++; next} /S/ {count=0; next} /E/ {print count}' >> counted_introns_snp.txt


### EXTRACTING EXONS FROM GENES -----------------------------------------------------------------------------------

# CHROMOSOME 1 FEATURES
awk '/GENE/ || /exon/ {if ($1=="1") print $4}' gene_features.data >> chrom1features.txt

# CHROMOSOME 1 FEATURES
awk '/GENE/ || /exon/ {if ($1=="2") print $4}' gene_features.data >> chrom2features.txt

# CHROMOSOME 1 FEATURES
awk '/GENE/ || /exon/ {if ($1=="3") print $4}' gene_features.data >> chrom3features.txt

# CHROMOSOME 1 FEATURES
awk '/GENE/ || /exon/ {if ($1=="4") print $4}' gene_features.data >> chrom4features.txt

# CHROMOSOME 1 FEATURES
awk '/GENE/ || /exon/ {if ($1=="5") print $4}' gene_features.data >> chrom5features.txt



# UNIQ METHOD
uniq -c features.txt | awk 'BEGIN {print "GENES" " " "EXONS"} /GENE/ {line++} /exon/ {print "GENE_"line " " $1}' >> results.txt

# CALCULATION METHOD
awk 'BEGIN {print "GENE" " " "Exons"} /GENE/ {line++; count=0; next} /exon/ {count++} {print line" "count}' chrom*features.txt >> counted_exons.txt

# COUNTING SNPS ONLY
awk ' {print $1}  $3!="" {print $2"S"; print $3"E"} ' main_exons_data.txt | sort -n |   
awk 'BEGIN {print "SNPs"} !/E|S/ {count++; next} /S/ {count=0; next} /E/ {print count}' >> counted_snp.txt

# MERGING DATA INTO ONE TABLE
paste -d ' ' counted_exons.txt counted_exons_snp.txt chrom*_len.txt >> chrom1_exonsresults.txt

