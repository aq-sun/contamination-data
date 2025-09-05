#!/usr/bin/env Rscript

library("vegan")
library("ecodist")
library(ade4)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

dir.create(output_dir, showWarnings = FALSE)

labeled_files <- list.files(input_dir, pattern = "*.pcoa.txt", full.names = TRUE)

for (labeled_file in labeled_files) {
  # Get name of file (which tells us the site)
  name_in_list <- strsplit(labeled_file, split='/|\\.')
  name <- tail(name_in_list[[1]], 3)[1]

  # Read file as a dataframe
  labeled_df <- read.table(labeled_file, header=TRUE)
  unlabeled <- subset(labeled_df, select=-c(site))

  # Calculate a distance matrix
  dist <- vegdist(unlabeled, method = "bray")

  pcoa <- cmdscale(dist, eig = TRUE)
  pcoa_coords <- as.data.frame(pcoa$points) # Extract PCoA coordinates
  
  # Make and save plot as png
  png(file=file.path(paste0(output_dir, name, "_pcoa.png")), height=2000, width=2000, res=300)

  print(labeled_df$site)
  # Plot PCoA using ade4
  s.class(pcoa_coords, fac = as.factor(labeled_df$site), 
          col = c("red", "olivedrab3"),
          addaxes = FALSE)
  
  # Calculate the proportion of variance explained by PCoA1 and PCoA2
  pcoa1 <- 100 * (pcoa$eig[1] / sum(pcoa$eig))
  pcoa2 <- 100 * (pcoa$eig[2] / sum(pcoa$eig))
  mtext(paste0("PC1 (",round(pcoa1, digits=1), "%)"), side = 1, line = 2, cex = 1.2)
  mtext(paste0("PC2 (",round(pcoa2, digits=1), "%)"), side = 2, line = 2, cex = 1.2)
  dev.off() #call off png 

  # Run PERMANOVA
  permanova_result <- adonis2(dist ~ site, data = labeled_df)
  
  write.table(permanova_result, file.path(paste0(output_dir, name, "_permanova.txt")), sep="\t", quote = FALSE, row.names = FALSE) # Write tab delimiter text file
}