
# Starting with BAM files obtained from read alignment, we generate a raw count matrix
# by quantifying reads assigned to genomic features (genes).

# Load necessary libraries 

library (Rsubread); library (Rsamtools)


# Generate the 'mapped_statistics' file. 
# This step creates summary statistics tables based on the alignment data from .bam files

# setwd('F:/')

bam_files <- list.files("bam/", pattern= "bam")
bam_files <- paste0("bam/", bam_files)
bam_files


# Sort each BAM file by chromosome and start coordinate 
# to ensure all reads are ordered correctly for processing.

for (bam in bam_files) {
    out <- suppressWarnings(sortBam(bam, "temporal"))
  
    file.rename(out, bam)
}

# Alignment statistics: generate a summary table of mapping metrics for each .bam file

diagnostics <- list()
for (bam in bam_files) {
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
write.csv(diag.stats, "bam/mapped_statistics.csv")

# Index the .bam files by creating corresponding .bam.bai index files

indexBam(bam_files)


# Quantification of reads per genomic feature

diag.stats <- read.csv("bam/mapped_statistics.csv", header = T)
lib.sizes <- as.vector(diag.stats[,2])
filenames <- bam_files

fc <- featureCounts(files = filenames,
                    annot.ext = "gtf/gencode.v38.annotation.gtf",
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="gene",
                    GTF.attrType="gene_name",
                    isPairedEnd = TRUE,
                    strandSpecific = 2,
                    nthreads = 20)

write.table(fc$counts, file = "raw_counts.tsv", sep = "\t")
