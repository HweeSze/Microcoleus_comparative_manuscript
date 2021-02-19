library(readxl)
library(openxlsx)
library(dplyr)


read_length <- read.table("../all.read.length.txt",header = T)
contig_length = read.table('../list.scaffolds.in.all.assemblies-contig_length.txt')


S1 = read.table('../all.vs.S01-node_reads.txt', sep='\t', stringsAsFactors=F)
S2 = read.table('../all.vs.S02-node_reads.txt', sep='\t', stringsAsFactors=F)
S3 = read.table('../all.vs.S03-node_reads.txt', sep='\t', stringsAsFactors=F)
S4 = read.table('../all.vs.S04-node_reads.txt', sep='\t', stringsAsFactors=F)
S5 = read.table('../all.vs.S05-node_reads.txt', sep='\t', stringsAsFactors=F)
S6 = read.table('../all.vs.S06-node_reads.txt', sep='\t', stringsAsFactors=F)
S7 = read.table('../all.vs.S07-node_reads.txt', sep='\t', stringsAsFactors=F)
S8 = read.table('../all.vs.S08-node_reads.txt', sep='\t', stringsAsFactors=F)
S9 = read.table('../all.vs.S09-node_reads.txt', sep='\t', stringsAsFactors=F)
S10 = read.table('../all.vs.S10-node_reads.txt', sep='\t', stringsAsFactors=F)
S11 = read.table('../all.vs.S11-node_reads.txt', sep='\t', stringsAsFactors=F)


# Run for each sample
S1 = read.table('../all.vs.S1-node_reads.txt', sep='\t', stringsAsFactors=F)
colnames(S1) = c("Contig", "Number_of_reads")
S1$Contig_length = contig_length$V1
S1$Total_read_length = S1$Number_of_reads * read_length$S1

# Import the two data tables
cov_table_S1 = S1
bin_table = read.table('../../7b-Refined-bins/all-refined-derep-bins-divided-per-sample/all.scaffolds2bin.txt')
colnames(bin_table) = c('Contig', 'Bin')

# Join the two tables together based on the value of the Contig column, remove unbinned contigs, then drop the columns we don't care about
cov_table_S1 = left_join(cov_table_S1, bin_table, on=Contig) %>% filter(., !is.na(Bin) ) %>% select(., Bin, Contig, Contig_length, Total_read_length)

# For each Bin in the table, calculate the coverage then cast the result as a dataframe
bin_coverage_S1 = group_by(cov_table_S1, Bin) %>% summarize( cov_table_S1 = sum(Total_read_length) / sum(Contig_length) ) %>% as.data.frame()


#REPLACE USING RSCRIPT


# Make final table
final_cov_table = bin_coverage_S1
final_cov_table$bin_coverage_S2 = bin_coverage_S2$cov_table_S2
final_cov_table$bin_coverage_S3 = bin_coverage_S3$cov_table_S3
final_cov_table$bin_coverage_S4 = bin_coverage_S4$cov_table_S4
final_cov_table$bin_coverage_S5 = bin_coverage_S5$cov_table_S5
final_cov_table$bin_coverage_S6 = bin_coverage_S6$cov_table_S6
final_cov_table$bin_coverage_S7 = bin_coverage_S7$cov_table_S7
final_cov_table$bin_coverage_S8 = bin_coverage_S8$cov_table_S8
final_cov_table$bin_coverage_S9 = bin_coverage_S9$cov_table_S9
final_cov_table$bin_coverage_S10 = bin_coverage_S10$cov_table_S10
final_cov_table$bin_coverage_S11 = bin_coverage_S11$cov_table_S11

write.xlsx(x = final_cov_table, file = "../genome-coverage-table.xlsx", col.names = T, row.names = F)
