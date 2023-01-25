library(tidyverse)
library(biokit)
library(org.Hs.eg.db)
library(decoupleR)
library(igraph)


#### OMNIPATH ####

# load prior knowledge network from Omnipath
pkn_initial <- OmnipathR::import_omnipath_interactions() %>%
  dplyr::transmute(source = source_genesymbol, 
                   target = target_genesymbol,
                   mor = case_when(consensus_stimulation == 1 & consensus_inhibition == 0 ~ 1,
                                   consensus_stimulation == 0 & consensus_inhibition == 1 ~ -1,
                                   TRUE ~ 0)) %>%
  dplyr::filter(source != target & mor != 0) %>%
  mutate(int_id = ifelse(mor == 1, paste0(source, "-->", target), paste0(source, "--|", target))) %>%
  distinct()

# retain the giant connected component from the graph
g <- igraph::graph_from_data_frame(pkn_initial)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g_filtered <- igraph::induced_subgraph(g, vert_ids)
pkn_filtered <- as_data_frame(g_filtered) %>%
  dplyr::rename(source = from, target = to)

#### DOROTHEA ####

# load dorothea and filter to TFs inside the network
dorothea_abc <- decoupleR::get_dorothea() %>%
  filter(source %in% pkn_filtered$source | source %in% pkn_filtered$target)

#### TIME-COURSE DATA PREPROCESSING ####

### data downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127988

input_rpm_file <- "GSE127988_RPM_Gerosa_et_al.csv.gz"

# read metadata, data, and format it
metadata <- read_csv(gzfile(input_rpm_file), n_max = 5, col_names  = FALSE) %>%
  dplyr::select(-X1) %>%
  column_to_rownames("X2") %>%
  t() %>%
  as.data.frame() %>%
  mutate(col_id = paste0("C", 1:nrow(.)))
colnames(metadata) <- make.names(colnames(metadata))
write_tsv(as.data.frame(metadata), 'experimental_metadata.tsv')

rpm_data <- read_csv(gzfile(input_rpm_file), skip = 7, col_names = FALSE ) %>%
  dplyr::select(-X2) %>%
  column_to_rownames("X1") %>%
  as.matrix()

colnames(rpm_data) <- metadata$col_id

# translate rpm data to gene symbols
rpm_symbols <- translateMatrixWithDb(mat = rpm_data, db = org.Hs.eg.db, sourceKey = "ENSEMBL", targetKey = "SYMBOL")

# select a single time-course dataset
int_metadata <- metadata %>%
  dplyr::filter(Type == "Individual" &Vemurafenib.24h..uM. == 1 & Cobimetinib.24h..uM. == 1) %>%
  dplyr::rename(egf_time = EGF.100ng.mL..h.)
int_rpm <- rpm_data[, int_metadata$col_id]
time_points <- unique(int_metadata$egf_time)

write_tsv(as.data.frame(int_metadata), 'metadata.tsv')


# define contrasts as consecutive temporal contrasts
contrast_df <- data.frame()
for( i in 1:(length(time_points)-1) ) {
  
  to_add_df <- data.frame(
                          num_time = time_points[i+1],
                          den_time = time_points[i]
                          )
  contrast_df <- bind_rows(contrast_df, to_add_df)
  
}

# compute t values from RPM data
t_df_list <- lapply(1:nrow(contrast_df), function(int_row) {
  
  int_num <- as.character(contrast_df[int_row, "num_time"])
  int_den <- as.character(contrast_df[int_row, "den_time"])
  message("Comparing ", int_num, " vs ", int_den)
  
  # get col indexes
  num_ids <- int_metadata %>%
    dplyr::filter(egf_time == int_num) %>%
    dplyr::pull(col_id)
  den_ids <- int_metadata %>%
    dplyr::filter(egf_time == int_den) %>%
    dplyr::pull(col_id)
  
  # subset expression matrix to interesting samples and perform pairwise comparisons
  int_mat <- rpm_symbols[, c(num_ids, den_ids)]
  t_values <- apply(int_mat, 1, function(x) {
    
    t_value <- t.test(x[num_ids], x[den_ids])$statistic
    return(t_value)
    
  })
  
  out_df <- data.frame(row.names = names(t_values), t = t_values)
  colnames(out_df) <- paste0("t", int_num, "_t", int_den)
  
  return(out_df)
  
})
t_matrix <- Reduce( bind_cols, t_df_list) %>%
  as.matrix()

# remove rows with NAs
keep <- apply(t_matrix, 1, function(x) !any(is.na(x)))
t_matrix <- t_matrix[keep,]

write_tsv(as.data.frame(t_matrix) %>% rownames_to_column(), 't_matrix_de.tsv')

#### PERFORM DECOUPLER ANALYSIS ####
decoupler_results <- decoupleR::run_wmean(mat = t_matrix, network = dorothea_abc, seed = 149) %>%
  dplyr::filter(statistic == "norm_wmean")

write_tsv(as.data.frame(decoupler_results), 'tfs.tsv')

##### PREPARE CORNETO INPUT #####
corneto_meas <- decoupler_results %>%
  group_by(condition) %>%
  arrange(desc(abs(score))) %>%
  #slice_head(n = 50) %>%
  ungroup() %>%
  transmute(condition, type = "measurement", feature = source, score) %>%
  distinct()
corneto_pert <- data.frame(condition = unique(corneto_meas$condition),
                           type = "perturbation",
                           feature = "EGF",
                           score = 1)
corneto_input <- bind_rows(corneto_pert, corneto_meas)
corneto_pkn <- pkn_filtered %>%
  dplyr::select(-int_id) %>%
  transmute(source, mor, target)

# save the network and the corneto input
write_tsv(corneto_pkn, "pkn.tsv")
write_tsv(corneto_input, "data.tsv")


