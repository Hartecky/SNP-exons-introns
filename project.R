### CHARACTERISTICS OF SINGLE NUCLEOTIDE POLYMORHPHISM IN ARABIDOPSIS THALIANA 

### MATERIALS | METHODS

# ----------------------------------------------------------------------------

### DATA PREPARATION

### READ DATA ABOUT SNP
chrom1 = read.csv("/home/bartlomiej/Pulpit/skn_proj/chrom1/chrom1_exonresults.txt" , sep = " ", header = TRUE)
chrom2 = read.csv("/home/bartlomiej/Pulpit/skn_proj/chrom2/chrom2_exonresults.txt" , sep = " ", header = TRUE)
chrom3 = read.csv("/home/bartlomiej/Pulpit/skn_proj/chrom3/chrom3_exonresults.txt" , sep = " ", header = TRUE)
chrom4 = read.csv("/home/bartlomiej/Pulpit/skn_proj/chrom4/chrom4_exonresults.txt" , sep = " ", header = TRUE)
chrom5 = read.csv("/home/bartlomiej/Pulpit/skn_proj/chrom5/chrom5_exonresults.txt" , sep = " ", header = TRUE)

### REDUCE NA VALUES
chrom1[is.na(chrom1)] <- 0

### READ GENE ID'S
GID1 = read.table("/home/bartlomiej/Pulpit/skn_proj/chrom1/GENE_ID1.txt", header = FALSE)
GID2 = read.table("/home/bartlomiej/Pulpit/skn_proj/chrom2/GENE_ID2.txt", header = FALSE)
GID3 = read.table("/home/bartlomiej/Pulpit/skn_proj/chrom3/GENE_ID3.txt", header = FALSE)
GID4 = read.table("/home/bartlomiej/Pulpit/skn_proj/chrom4/GENE_ID4.txt", header = FALSE)
GID5 = read.table("/home/bartlomiej/Pulpit/skn_proj/chrom5/GENE_ID5.txt", header = FALSE)

colnames(GID1) = c("GENE_ID")
colnames(GID2) = c("GENE_ID")
colnames(GID3) = c("GENE_ID")
colnames(GID4) = c("GENE_ID")
colnames(GID5) = c("GENE_ID")

### READ CHROMOSOMES IN VCF
chrom1.vcf <- read.table("/home/bartlomiej/Pulpit/skn_proj/chrom1.txt", header = TRUE)
chrom2.vcf <- read.table("/home/bartlomiej/Pulpit/skn_proj/chrom2.txt", header = TRUE)
chrom3.vcf <- read.table("/home/bartlomiej/Pulpit/skn_proj/chrom3.txt", header = TRUE)
chrom4.vcf <- read.table("/home/bartlomiej/Pulpit/skn_proj/chrom4.txt", header = TRUE)
chrom5.vcf <- read.table("/home/bartlomiej/Pulpit/skn_proj/chrom5.txt", header = TRUE)
genome.info <- read.table("/home/bartlomiej/Pulpit/skn_proj/snpmain_info.txt", header = FALSE)
### CHANGE DATA TYPES AND MAKE NECESSERY VARIABLES
chrom1.vcf$ID <- as.character(chrom1.vcf$ID)
chrom2.vcf$ID <- as.character(chrom2.vcf$ID)
chrom3.vcf$ID <- as.character(chrom3.vcf$ID)
chrom4.vcf$ID <- as.character(chrom4.vcf$ID)
chrom5.vcf$ID <- as.character(chrom5.vcf$ID)

chrom1$GENE = as.factor(chrom1$GENE)
chrom1$Exons = as.factor(chrom1$Exons)
chrom1$sum_len <- 0

chrom2$GENE = as.factor(chrom2$GENE)
chrom2$Exons = as.factor(chrom2$Exons)
chrom2$sum_len <- 0

chrom3$GENE = as.factor(chrom3$GENE)
chrom3$Exons = as.factor(chrom3$Exons)
chrom3$sum_len <- 0

chrom4$GENE = as.factor(chrom4$GENE)
chrom4$Exons = as.factor(chrom4$Exons)
chrom4$sum_len <- 0

chrom5$GENE = as.factor(chrom5$GENE)
chrom5$Exons = as.factor(chrom5$Exons)
chrom5$sum_len <- 0

### COUNT GENES FOR EACH CHROMOSOME
all.genes1 <- length(unique(chrom1[, 1]))
all.genes2 <- length(unique(chrom2[, 1]))
all.genes3 <- length(unique(chrom3[, 1]))
all.genes4 <- length(unique(chrom4[, 1]))
all.genes5 <- length(unique(chrom5[, 1]))

### COUNT EXONS LENGTH SUM FOR EACH GENE IN EACH CHROMOSOME

#Chrom 1
for (gene in 1:all.genes1){
  suma <- sum(chrom1[chrom1$GENE == gene, 4])
  chrom1[chrom1$GENE == gene, 5] <- suma
}

#Chrom 2
for (gene in 1:all.genes2){
  suma <- sum(chrom2[chrom2$GENE == gene, 4])
  chrom2[chrom2$GENE == gene, 5] <- suma
}

# Chrom 3
for (gene in 1:all.genes3){
  suma <- sum(chrom3[chrom3$GENE == gene, 4])
  chrom3[chrom3$GENE == gene, 5] <- suma
}

# Chrom 4
for (gene in 1:all.genes4){
  suma <- sum(chrom4[chrom4$GENE == gene, 4])
  chrom4[chrom4$GENE == gene, 5] <- suma
}

# Chrom 5
for (gene in 1:all.genes5){
  suma <- sum(chrom5[chrom5$GENE == gene, 4])
  chrom5[chrom5$GENE == gene, 5] <- suma
}

### ADD GENE ID TO THE DATAFRAME
GID1$GENE <- as.factor(1:nrow(GID1))
GID2$GENE <- as.factor(1:nrow(GID2))
GID3$GENE <- as.factor(1:nrow(GID3))
GID4$GENE <- as.factor(1:nrow(GID4))
GID5$GENE <- as.factor(1:nrow(GID5))

### MERGE DATA WITH GENE ID FOR EACH CHROMOSOME | SORT IT

chrom1 <- merge(chrom1, GID1, by = "GENE", by.y = "GENE", all.x = TRUE, no.dups = FALSE)
chrom1 <- chrom1[order(chrom1$GENE), ]

chrom2 <- merge(chrom2, GID2, by = "GENE", by.y = "GENE", all.x = TRUE, no.dups = FALSE)
chrom2 <- chrom2[order(chrom2$GENE), ]

chrom3 <- merge(chrom3, GID3, by = "GENE", by.y = "GENE", all.x = TRUE, no.dups = FALSE)
chrom3 <- chrom3[order(chrom1$GENE), ]

chrom4 <- merge(chrom4, GID1, by = "GENE", by.y = "GENE", all.x = TRUE, no.dups = FALSE)
chrom4 <- chrom4[order(chrom4$GENE), ]

chrom5 <- merge(chrom5, GID5, by = "GENE", by.y = "GENE", all.x = TRUE, no.dups = FALSE)
chrom5 <- chrom5[order(chrom5$GENE), ]


### COUNT SNP SCORE STATISTICS
chrom1$score = (chrom1$SNP*chrom1$exon_len) / chrom1$sum_len
chrom2$score = (chrom2$SNP*chrom2$exon_len) / chrom2$sum_len
chrom3$score = (chrom3$SNP*chrom3$exon_len) / chrom3$sum_len
chrom4$score = (chrom4$SNP*chrom4$exon_len) / chrom4$sum_len
chrom5$score = (chrom5$SNP*chrom5$exon_len) / chrom5$sum_len

# ------------------------------------------------------------------------------- 
### CHECK THE DIFFERENCES BETWEEN GROUPS IN SNP VALUES
### LOAD PACKAGES
library(lawstat)
library(ggpubr)
library(dunn.test)
library(rcompanion)
library(multcompView)
# ------------------------------------------------------------------------------

### CHROM 1

###NORMALITY TEST
variance.test <- levene.test(chrom1$score, chrom1$Exons)

### NON-PARAMETRIC ANOVA
###KRUSKAL TEST
non.param.kruskal <- kruskal.test(chrom1$SNP~chrom1$Exons)

### DUNN TEST FOR CHECK SIGNIF DIFFERENCES BEETWEN GROUPS
non.param.dunn <- dunn.test(chrom1$SNP,chrom1$Exons)
non.param.dunn$comparisons

### PAIRWISE WILCOX TEST
PT = pairwise.wilcox.test(chrom1$SNP,chrom1$Exons, p.adjust.method = "fdr")
PT = PT$p.value
PT1 = fullPTable(PT)

### LETTER CODED DIFFERENCES
letter.coded.diffs <- multcompLetters(PT,
                        compare="<",
                        threshold=0.05,
                        Letters=letters,
                        reversed = FALSE)

### SCORE STATISTICS 
variance.score.test <- levene.test(chrom1$score, chrom1$Exons)

### PAIRWISE WILCOX TEST
ST = pairwise.wilcox.test(chrom1$score, chrom1$Exons, p.adjust.method = "fdr")
ST = ST$p.value
ST1 = fullPTable(ST)

letter.score.coded.diffs <- multcompLetters(ST,
                                            compare = "<",
                                            threshold = 0.05,
                                            Letters = letters,
                                            reversed = FALSE)

# ---------------------------------------------------------------------------------

### CHROM 2

variance.test2 <- levene.test(chrom2$SNP, chrom2$Exons)

### Non parametric ANOVA
# Kruskal test
non.param.kruskal2 <- kruskal.test(chrom2$SNP~chrom2$Exons)

# Dunn test for signif differences beetwen groups
non.param.dunn2 <- dunn.test(chrom2$SNP,chrom2$Exons)
non.param.dunn2$comparisons

# Pairwise Wilcox Test for SNP
PT2 = pairwise.wilcox.test(chrom2$SNP,chrom2$Exons, p.adjust.method = "fdr")
PT2 = PT2$p.value
PT.chr2 = fullPTable(PT2)

# LETTER CODED DIFFERENCES
letter.coded.diffs2 <- multcompLetters(PT2,
                                      compare="<",
                                      threshold=0.05,
                                      Letters=letters,
                                      reversed = FALSE)

# Score statistics 
variance.score.test2 <- levene.test(chrom2$score, chrom2$Exons)

# Pairwise Wilcox Test for SNP
ST2 = pairwise.wilcox.test(chrom2$score, chrom2$Exons, p.adjust.method = "fdr")
ST2 = ST2$p.value
ST.chr2 = fullPTable(ST2)

letter.score.coded.diffs2 <- multcompLetters(ST2,
                                            compare = "<",
                                            threshold = 0.05,
                                            Letters = letters,
                                            reversed = FALSE)
# ------------------------------------------------------------------------------

# CHROM 3

# Non parametric ANOVA
# Kruskal test
non.param.kruskal3 <- kruskal.test(chrom3$SNP~chrom3$Exons)

# Dunn test for signif differences beetwen groups
non.param.dunn3 <- dunn.test(chrom3$SNP,chrom3$Exons)
non.param.dunn3$comparisons

# Pairwise Wilcox Test for SNP
PT3 = pairwise.wilcox.test(chrom3$SNP,chrom3$Exons, p.adjust.method = "fdr")
PT3 = PT3$p.value
PT.chr3 = fullPTable(PT3)

# LETTER CODED DIFFERENCES - tu nie dziala
letter.coded.diffs3 = multcompLetters(PT3, compare = "<", threshold = 0.05, Letters = letters, reversed = FALSE)

# Score statistics 
variance.score.test3 <- levene.test(chrom3$score, chrom3$Exons)

# Pairwise Wilcox Test for SNP
ST3 = pairwise.wilcox.test(chrom3$score, chrom3$Exons, p.adjust.method = "fdr")
ST3 = ST3$p.value
ST.chr3 = fullPTable(ST3)

letter.score.coded.diffs3 <- multcompLetters(ST3,
                                            compare = "<",
                                            threshold = 0.05,
                                            Letters = letters,
                                            reversed = FALSE)
# ------------------------------------------------------------------------------------------

# CHROM 4

# Non parametric ANOVA
# Kruskal test
non.param.kruskal4 <- kruskal.test(chrom4$SNP~chrom4$Exons)

# Dunn test for signif differences beetwen groups
non.param.dunn4 <- dunn.test(chrom4$SNP,chrom4$Exons)
non.param.dunn4$comparisons

# Pairwise Wilcox Test for SNP
PT4 = pairwise.wilcox.test(chrom4$SNP,chrom4$Exons, p.adjust.method = "fdr")
PT4 = PT4$p.value
PT.chr4 = fullPTable(PT4)

# LETTER CODED DIFFERENCES
letter.coded.diffs4 <- multcompLetters(PT4,
                                      compare="<",
                                      threshold=0.05,
                                      Letters=letters,
                                      reversed = FALSE)

# Score statistics 
variance.score.test4 <- levene.test(chrom4$score, chrom4$Exons)

# Pairwise Wilcox Test for SNP
ST4 = pairwise.wilcox.test(chrom4$score, chrom4$Exons, p.adjust.method = "fdr")
ST4 = ST4$p.value
ST.chr4 = fullPTable(ST4)

letter.score.coded.diffs4 <- multcompLetters(ST4,
                                            compare = "<",
                                            threshold = 0.05,
                                            Letters = letters,
                                            reversed = FALSE)

# -----------------------------------------------------------------------------------

# CHROM 5

# Non parametric ANOVA
# Kruskal test
non.param.kruskal5 <- kruskal.test(chrom5$SNP~chrom5$Exons)

# Dunn test for signif differences beetwen groups
non.param.dunn5 <- dunn.test(chrom5$SNP,chrom5$Exons)
non.param.dunn5$comparisons

# Pairwise Wilcox Test for SNP
PT5 = pairwise.wilcox.test(chrom5$SNP,chrom5$Exons, p.adjust.method = "fdr")
PT5 = PT5$p.value
PT.chr5 = fullPTable(PT5)
  
# LETTER CODED DIFFERENCES
letter.coded.diffs5 <- multcompLetters(PT5,
                                      compare="<",
                                      threshold=0.05,
                                      Letters=letters,
                                      reversed = FALSE)

# Score statistics 
variance.score.test5 <- levene.test(chrom5$score, chrom5$Exons)

# Pairwise Wilcox Test for SNP
ST5 = pairwise.wilcox.test(chrom5$score, chrom5$Exons, p.adjust.method = "fdr")
ST5 = ST5$p.value
ST.chr5 = fullPTable(ST5)

letter.score.coded.diffs5 <- multcompLetters(ST5,
                                            compare = "<",
                                            threshold = 0.05,
                                            Letters = letters,
                                            reversed = FALSE)
### ------------------------------------------------------------------------------

### DESCRIPTIVE STATISTICS
### WHOLE GENOME 

all.variants <- sum(chrom1.variants.count,chrom2.variants.count,chrom3.variants.count,chrom4.variants.count,chrom5.variants.count)

### CHROM1
chrom1.variants.count <- length(chrom1.vcf$Pos)
chrom1.info <- summary(chrom1.vcf$Info)
chrom1.snps <- summary(chrom1$SNP)

### CHROM 2
chrom2.variants.count <- length(chrom2.vcf$Pos)
### CHROM 3
chrom3.variants.count <- length(chrom3.vcf$Pos)
### CHROM 4
chrom4.variants.count <- length(chrom4.vcf$Pos)
### CHROM 5
chrom5.variants.count <- length(chrom5.vcf$Pos)
