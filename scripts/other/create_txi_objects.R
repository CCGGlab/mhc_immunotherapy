library(tidyverse)

# RSEM output folders per study
folders <- c(
  "Riaz_2017",
  "Liu_2019",
  "Hugo_2016",
  "Gide_2019"
)

parse_rsem <- function(study_name) {
  folder_basename <- study_name
  input_folder <- paste0("temp/RSEM_output/", folder_basename)
  input_files <- list.files(input_folder, pattern = "\\.genes.results$", full.names = TRUE)
  
  # Exclude problematic sample for Hugo
  if (study_name == "Hugo_2016") {
    # For Pt27 there are two samples available
    # (one "upper arm", the other "upper back", exclude the second sample)
    input_files <- purrr::discard(input_files, str_detect, "SRR3184299")
  }
  
  # Obtain accession numbers from file names 
  samples_lst <- basename(input_files) %>%
    str_remove("\\.genes\\.results") %>%
    as_tibble_col
  
  # Map accession numbers to patient IDs according to metadata
  # - Different for Riaz (as we need to exclude "On" treatment samples)
  if (study_name == "Riaz_2017") {
    # PRJNA356761
    metadata_rna <- read_tsv("downloads/pub/riaz_2017/metadata_rna_conv.tsv", col_names = F) %>%
      rename(Run := 1, title := 2) %>%
      transmute(Run, `Sample Name` = str_match(title, "^[^_]+_[^_]+")[,1])

    pt <- samples_lst %>%
      left_join(metadata_rna, by=c("value" = "Run")) %>%
      mutate(sample_id = `Sample Name`, `Sample Name`=NULL) %>%
      separate(sample_id, c("patient", "sample_type"), sep="_") %>%
      transmute(sample_id = value, patient, is_pre = (sample_type == "Pre"))

    # Only include samples before treatment (Pre)
    input_files <- input_files[pt$is_pre]
    names(input_files) <- pt[pt$is_pre,] %>%
      select(sample_id, patient) %>%
      deframe
  } else {
    # For other studies
    source(paste0("scripts/functions/study_sample_info/", study_name, ".R"), local = T)
    pt <- get_sample_info(samples_lst)
    names(input_files) <- deframe(select(pt, sample_id, patient))
  }

  # Read RSEM files
  txi.rsem <- tximport::tximport(input_files, type = "rsem", txIn = FALSE, txOut = FALSE)
}

# Apply parse_rsem for all studies (output folders)
txi_objects <- map(folders, parse_rsem)

saveRDS(txi_objects, "data/txi_objects_rsem.rds")
