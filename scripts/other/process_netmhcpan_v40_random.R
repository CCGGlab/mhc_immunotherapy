# MHC-I
# Note: Run this script from the terminal with "Rscript"
# In RStudio a "Too many open files error" might occur
# while processing thousands of output files
# Takes about 8 GB of RAM and 50s to calculate
library(tidyverse)
library(multidplyr)

# For parallel processing
cluster <- multidplyr::new_cluster(15)
cluster_library(cluster, "dplyr")

# Where the NetMHCIIpan output files were saved
output_folder <- "temp/random_nmp/netMHCpan-4.0_output"

# Create a data frame with the path to the NetMHCIIpan output files
# and the corresponding alleles
output_files <- list.files(output_folder, full.names = T) %>%
  # Call the column with the file paths "path"
  as_tibble_col("path") %>%
  # Obtain the allele name from the file path
  mutate(allele = str_remove(basename(path), '_[^_]*\\.xls$'))

result <- output_files %>%
  group_by(allele) %>%
  # Process result files in parallel
  partition(cluster) %>%
  # Steps to do per allele
  summarise({
    # Read the output TSV files
    vroom::vroom(
      path,
      delim = "\t",
      quote = "",
      skip = 1,
      num_threads = 8,
      show_col_types = F
    ) %>%
      # Save memory:
      # Only select the columns we need
      select(ID, nM, Rank)
  }) %>%
  # Perform the actual calculation (in parallel)
  collect()

mhc1_Kd_best <- result %>%
  pivot_wider(id_cols = "ID",
              names_from = "allele",
              # values_from = "min(nM)",
              values_from = "nM"
  ) %>%
  # Make order of rows match with the original GPPM_rand matrix
  # (allows comparing both)
  mutate(ID = as.numeric(ID)) %>%
  arrange(ID) %>%
  select(-ID) %>%
  as.matrix

saveRDS(mhc1_Kd_best, "temp/mhc1_rand_matrix.rds")