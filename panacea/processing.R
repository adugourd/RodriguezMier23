# libraries
library(tidyverse)
library(readxl)
library(limma)
library(edgeR)
library(biokit)
library(org.Hs.eg.db)
library(decoupleR)
library(OmnipathR)
library(igraph)
library(jsonlite)


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

#### PANACEA DATA PREPROCESSING ####

## DREAM COUNTS FILE WERE DOWNLOADED FROM THE FOLLOWING URL:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186341

# read from csv files
counts_file <- "GSE186341_DU145_dream_counts.csv.gz"
raw_content <- read_csv(counts_file)
count_matrix <- raw_content %>%
  column_to_rownames(var = "...1") %>%
  as.matrix()

# extract metadata from column names for every dataset
metadata <- data.frame(id = colnames(count_matrix)) %>%
    separate(id, into = c("treatment", "dose", "time"), sep = "_", remove = FALSE) %>%
    mutate(time = str_replace(time, "\\.\\..*", ""))

# apply processing
translated_matrix <- biokit::translateMatrixWithDb(count_matrix, db = org.Hs.eg.db, sourceKey = "ENTREZID", targetKey = "SYMBOL")
  
# prune low expressed genes using edgeR
design_matrix <- model.matrix(~ 0 + treatment, data = metadata)
colnames(design_matrix) <- str_replace(colnames(design_matrix), "treatment", "")
keep <- edgeR::filterByExpr(y = translated_matrix, design = design_matrix)
filtered_matrix <- translated_matrix[keep,]
  
# normalize using TMMs + CPMs
norm_matrix <- biokit::countsToTmm(filtered_matrix, log = TRUE, prior.count = 1)

#### PANACEA DATA DIFFERENTIAL EXPRESSION ANALYSIS ####

# define contrasts
int_contrasts <- unique(metadata$treatment) %>%
  .[. != "DMSO" & . != "UNTREATED"] %>%
  paste0(., "-DMSO")
contrast_matrix <- makeContrasts(contrasts = int_contrasts, levels = design_matrix)

# limma differential expression analysis
fit <- lmFit(norm_matrix, design = design_matrix) %>%
  contrasts.fit(., contrasts = contrast_matrix) %>%
  eBayes(.) 
de_table <- lapply(int_contrasts, function(c) {
  limma::topTable(fit, coef = c, number = Inf) %>%
    rownames_to_column(var = "gene") %>%
    mutate(contrast = c)
}) %>%
  bind_rows()

#### PANACEA DATA TF ANALYSIS WITH DECOUPLER ####

# input
decoupler_input <- de_table %>%
  dplyr::transmute(gene, id = str_replace(contrast, "-", "_"), t) %>%
  pivot_wider(values_from = t, names_from = id) %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()
decoupler_results <- decoupleR::run_wmean(mat = decoupler_input, network = dorothea_abc, seed = 149) %>%
  dplyr::filter(statistic == "norm_wmean")

#### EXCEL METADATA FILE DOWNLOADED FROM PAPER SUPPLEMENTARY MATERIALS
#### https://ars.els-cdn.com/content/image/1-s2.0-S2666379121003694-mmc3.xlsx

# read excel metadata
targets_file <- "1-s2.0-S2666379121003694-mmc3.xlsx"
targets_table <- read_xlsx(path = targets_file)

# fix typos in the resulting data
typos_tofix <- c("ERRB2" = "ERBB2",
                "AKK1" = "ATK1")
for(i in 1:length(typos_tofix)) {
  targets_table$`Drug-Targets` <- str_replace(string = targets_table$`Drug-Targets`, 
                                              pattern = names(typos_tofix)[i], replacement = typos_tofix[i])
}

# add targets to the decoupler table
combined_df <- decoupler_results %>%
  mutate(compound = str_replace(condition, "_DMSO", "")) %>%
  left_join(., transmute(targets_table, compound = toupper(...1), target = `Drug-Targets`), by = "compound")

# filter tfs that are not in the prior knowledge
filt_combined_df <- combined_df %>%
  subset(source %in% pkn_filtered$source | source %in% pkn_filtered$target)

##### PREPARE INPUT DATA #####

# format data for corneto
corneto_meas <- filt_combined_df %>%
  group_by(compound) %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  transmute( compound, type = "measurement", feature = source, score, target) %>%
  distinct()
corneto_pert <- corneto_meas %>%
  dplyr::select(compound, target) %>%
  mutate(target = str_split(target, ",")) %>%
  unnest(target) %>%
  transmute(compound, type = "perturbation", feature = target, score = -1)
corneto_input <- corneto_meas %>%
  dplyr::select(-target) %>%
  bind_rows(., corneto_pert)
corneto_pkn <- pkn_filtered %>%
  dplyr::select(-int_id) %>%
  transmute(source, mor, target)

# visualize input
p <- ggplot(corneto_input, aes(x = compound, y = feature, fill = score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
  facet_grid(rows = vars(type), space = "free", scales = "free") +
  guides(x = guide_axis(angle = 60)) +
  ylab("Feature") +
  theme(axis.text.y = element_blank(), axis.ticks.y= element_blank())
ggsave("input_viz.png", p, width = 10, height = 6, dpi = 150)

# save the network and the corneto input
write_tsv(corneto_pkn, "pkn.tsv")
write_tsv(corneto_input, "data.tsv")