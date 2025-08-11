
# ALIGNMENT OF FASTQ FILES


# Load libraries

library(Rsubread); library(Rsamtools)

fastq.files.R1 <- list.files(path = "fastq/", pattern = "R1_001")
fastq.files.R2 <- list.files(path = "fastq/", pattern = "R2_001")
bam.files <- gsub(pattern = "_R1_001.fastq.gz", replacement = ".bam", fastq.files.R1)

for (i in 40:1) {
  bam.files <- gsub(pattern = paste0("_S",i), replacement = "", bam.files)
}
bam.files <- gsub(pattern = "__", replacement = "_", bam.files)
bam.files <- gsub(pattern = "-", replacement = ".", bam.files)
bam.files <- gsub(pattern = "R1", replacement = "rep1", bam.files) 
bam.files <- gsub(pattern = "R2", replacement = "rep2", bam.files) 


fastq.files.R1 <- paste0("fastq/",fastq.files.R1)
fastq.files.R2 <- paste0("fastq/",fastq.files.R2) 
bam.files <- paste0("bam/", bam.files)

# Aligning FASTQ files to the human genome (hg38) using Rsubread.
# First, build the genome index (only needed once).
# Then, perform paired-end alignment and generate BAM files.

##Index --> hg38
#setwd("~/genomes/hg38/")
#buildindex(basename="hg38_index_Rsubread",reference="fasta/hg38.fa")


# Alignment

align(index="hg38_index_Rsubread/hg38_index_Rsubread",
      readfile1=fastq.files.R1,
      readfile2 = fastq.files.R2,
      TH1=2,
      type=1,
      input_format="gzFASTQ",
      output_file=bam.files,
      unique = TRUE,
      nthreads = 20) 

# Sorting

for (bam in bam.files) {
  out <- suppressWarnings(sortBam(bam, "temporal"))
  file.rename(out, bam)
}

# Alignment Statistics

diagnostics <- list()
for (bam in bam.files) {
  total <- countBam(bam)$records
  mapped <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE)))$records
  marked <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
  diagnostics[[bam]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics)) 
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats
write.table(diag.stats, "bam/mapped_statistics.tsv", sep="\t", dec = ".")

# Index

indexBam(bam.files)


# Aligning FASTQ files to the Saccharomyces cerevisiae (sacCer2) genome for spike-in normalization.
# First, build the genome index (only once).
# Then, align paired-end reads, sort and index BAM files, and extract mapping statistics.

##Index --> Saccharomyces cerevisiae (sacCer2)
#setwd("~/genomes/sacCer2/")
#buildindex(basename="sacCer2_index_Rsubread",reference="sacCer2.fa")


# Spike-in

fastq.files.R1 <- list.files(path = "fastq/", pattern = "R1_001")
fastq.files.R2 <- list.files(path = "fastq/", pattern = "R2_001")
bam.files.sacCer <- gsub(pattern = "_R1_001.fastq.gz",replacement = "_spike-in.bam", fastq.files.R1)
bam.files.sacCer <- gsub(pattern = "_001", replacement = "", bam.files.sacCer)
for (i in 40:1) {
  bam.files.sacCer <- gsub(pattern = paste0("_S",i), replacement = "", bam.files.sacCer)
}
bam.files.sacCer <- gsub(pattern = "__", replacement = "_", bam.files.sacCer)
bam.files.sacCer <- gsub(pattern = "-", replacement = ".", bam.files.sacCer)
bam.files.sacCer <- gsub(pattern = "spike.in", replacement = "spike-in", bam.files.sacCer)
bam.files.sacCer <- gsub(pattern = "R1", replacement = "rep1", bam.files.sacCer)
bam.files.sacCer <- gsub(pattern = "R2", replacement = "rep2", bam.files.sacCer)



fastq.files.R1 <- paste0("fastq/",fastq.files.R1)
fastq.files.R2 <- paste0("fastq/",fastq.files.R2) 
bam.files.sacCer <- paste0("bam/", bam.files.sacCer)

# Alignment

align(index="sacCer2_index_Rsubread", readfile1=fastq.files.R1, readfile2 = fastq.files.R2, TH1=2, type=1, input_format="gzFASTQ", output_file=bam.files.sacCer, unique = TRUE, nthreads = 20) 

# Sorting

for (bam in bam.files.sacCer) {
  out <- suppressWarnings(sortBam(bam, "temporal"))
  file.rename(out, bam)
}

# Alignment Statistics

diagnostics <- list()
for (bam in bam.files.sacCer) {
  total <- countBam(bam)$records
  mapped <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE)))$records
  marked <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
  diagnostics[[bam]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics)) 
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats
write.table(diag.stats, "mapped_statistics_spike-in.tsv", sep="\t")

# Index

indexBam(bam.files.sacCer)



# Convert BAM files to BigWig format for genome browser visualization.
# File names are cleaned and standardized before running bamCoverage (deepTools).


bam_dir <- "bam"
bw_dir <- "bw"

fastq.files.R1 <- list.files(path = "fastq/", pattern = "R1_001")
bam.files <- gsub(pattern = "_R1_001.fastq.gz", replacement = ".bam", fastq.files.R1)
bam.files <- gsub(pattern = "_001", replacement = "", bam.files)
for (i in 40:1) {
  bam.files <- gsub(pattern = paste0("_S",i), replacement = "", bam.files)
}
bam.files <- gsub(pattern = "__", replacement = "_", bam.files)
bam.files <- gsub(pattern = "-", replacement = ".", bam.files)
bam.files <- gsub(pattern = "R1", replacement = "rep1", bam.files)
bam.files <- gsub(pattern = "R2", replacement = "rep2", bam.files)

for (bam_file in bam.files) {
  bw_file <- file.path(bw_dir, paste0(basename(tools::file_path_sans_ext(bam_file)), ".bw"))
  
  command <- paste("bamCoverage -p 20 --normalizeUsing RPKM -b", paste0(bam_dir, "/", bam_file), "-o", bw_file)
  
  system(command)
}

