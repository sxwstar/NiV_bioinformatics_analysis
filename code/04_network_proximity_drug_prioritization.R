rm(list = ls())
gc()

library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(igraph)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(parallel)

set.seed(123)

setwd("/users/hzhang1/project/network")

data_dir <- "data"
out_dir  <- "result"
if (!dir.exists(out_dir)) dir.create(out_dir)

topk_dir <- file.path(out_dir, "topk_distance")
if (!dir.exists(topk_dir)) dir.create(topk_dir)

############################################################
# 新增：网络作图专用输出目录
############################################################
network_table_dir <- file.path(out_dir, "network_plot_tables")
if (!dir.exists(network_table_dir)) dir.create(network_table_dir)

############################################################
# 1. NiV disease module
############################################################
niv <- read_excel(file.path(data_dir, "Nipah_PPI.xlsx"))

niv_genes <- niv %>%
  dplyr::mutate(Human_GeneID = as.character(Human_GeneID)) %>%
  tidyr::separate_rows(Human_GeneID, sep = ";") %>%
  dplyr::filter(!is.na(Human_GeneID), Human_GeneID != "") %>%
  dplyr::mutate(Human_GeneID = trimws(Human_GeneID)) %>%
  dplyr::distinct(Human_GeneID) %>%
  dplyr::pull(Human_GeneID)

############################################################
# 2. PPI network
############################################################
ppi <- read_excel(file.path(data_dir, "human.xlsx"))

edges <- ppi %>%
  dplyr::select(x_id, y_id) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%
  dplyr::distinct()

g <- igraph::graph_from_data_frame(edges, directed = FALSE)

niv_genes <- intersect(niv_genes, igraph::V(g)$name)
deg <- igraph::degree(g)

############################################################
# ⭐ degree阈值
############################################################
small_deg_threshold <- as.numeric(stats::quantile(deg, 0.1))

############################################################
# 3. drug-target + SYMBOL
############################################################
drug <- fread(file.path(data_dir, "drug_target.csv"))

drug_clean <- drug %>%
  dplyr::mutate(Entrezid = as.character(Entrezid)) %>%
  tidyr::separate_rows(Entrezid, sep = ";") %>%
  dplyr::filter(!is.na(Entrezid), Entrezid != "") %>%
  dplyr::mutate(Entrezid = trimws(Entrezid)) %>%
  dplyr::filter(grepl("^[0-9]+$", Entrezid))

id2symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(drug_clean$Entrezid),
  columns = "SYMBOL",
  keytype = "ENTREZID"
) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

symbol_map <- setNames(id2symbol$SYMBOL, id2symbol$ENTREZID)

drug_clean <- drug_clean %>%
  dplyr::left_join(id2symbol, by = c("Entrezid" = "ENTREZID"))

drug_list <- drug_clean %>%
  dplyr::group_by(drug_id) %>%
  dplyr::summarise(
    targets = list(unique(Entrezid)),
    target_symbols = list(stats::na.omit(unique(SYMBOL))),
    .groups = "drop"
  )

############################################################
# 新增：输出网络作图所需的 clean node/edge tables
# 只输出，不改变后续 proximity 计算逻辑
############################################################

# 3.1 NiV module genes clean
niv_symbol_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(niv_genes),
  columns = "SYMBOL",
  keytype = "ENTREZID"
) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

NiV_module_genes_clean <- data.frame(
  Entrezid = niv_genes,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(niv_symbol_df, by = c("Entrezid" = "ENTREZID")) %>%
  dplyr::mutate(
    node_type = "NiV_module_gene"
  ) %>%
  dplyr::arrange(Entrezid)

write.csv(
  NiV_module_genes_clean,
  file.path(network_table_dir, "NiV_module_genes_clean.csv"),
  row.names = FALSE
)

# 3.2 drug targets clean
drug_targets_clean <- drug_clean %>%
  dplyr::select(drug_id, Entrezid, SYMBOL) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    node_type = "drug_target"
  ) %>%
  dplyr::arrange(drug_id, Entrezid)

write.csv(
  drug_targets_clean,
  file.path(network_table_dir, "drug_targets_clean.csv"),
  row.names = FALSE
)

# 3.3 module internal PPI edges
module_internal_edges <- edges %>%
  dplyr::filter(x_id %in% niv_genes, y_id %in% niv_genes) %>%
  dplyr::mutate(edge_type = "module_module_ppi") %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(x_id = Entrezid, x_symbol = SYMBOL),
    by = "x_id"
  ) %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(y_id = Entrezid, y_symbol = SYMBOL),
    by = "y_id"
  ) %>%
  dplyr::select(
    x_id, x_symbol,
    y_id, y_symbol,
    edge_type
  ) %>%
  dplyr::distinct()

write.csv(
  module_internal_edges,
  file.path(network_table_dir, "module_internal_edges.csv"),
  row.names = FALSE
)

# 3.4 drug target to NiV module direct PPI edges
target_genes <- unique(drug_targets_clean$Entrezid)

target_module_edges_1 <- edges %>%
  dplyr::filter(x_id %in% target_genes, y_id %in% niv_genes) %>%
  dplyr::transmute(
    target_id = x_id,
    module_id = y_id,
    edge_type = "target_module_ppi"
  )

target_module_edges_2 <- edges %>%
  dplyr::filter(y_id %in% target_genes, x_id %in% niv_genes) %>%
  dplyr::transmute(
    target_id = y_id,
    module_id = x_id,
    edge_type = "target_module_ppi"
  )

target_module_edges_gene <- dplyr::bind_rows(
  target_module_edges_1,
  target_module_edges_2
) %>%
  dplyr::distinct()

drug_target_to_module_edges <- target_module_edges_gene %>%
  dplyr::left_join(
    drug_targets_clean %>%
      dplyr::select(drug_id, target_id = Entrezid, target_symbol = SYMBOL),
    by = "target_id"
  ) %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(module_id = Entrezid, module_symbol = SYMBOL),
    by = "module_id"
  ) %>%
  dplyr::select(
    drug_id,
    target_id, target_symbol,
    module_id, module_symbol,
    edge_type
  ) %>%
  dplyr::arrange(drug_id, target_id, module_id) %>%
  dplyr::distinct()

write.csv(
  drug_target_to_module_edges,
  file.path(network_table_dir, "drug_target_to_module_edges.csv"),
  row.names = FALSE
)

# 3.5 subnetwork edges for plotting
drug_target_edges_for_plot <- drug_targets_clean %>%
  dplyr::filter(Entrezid %in% unique(drug_target_to_module_edges$target_id)) %>%
  dplyr::transmute(
    from = paste0("Drug:", drug_id),
    to = paste0("Gene:", Entrezid),
    from_type = "drug",
    to_type = "drug_target",
    edge_type = "drug_target",
    drug_id = drug_id,
    target_id = Entrezid,
    module_id = NA_character_
  ) %>%
  dplyr::distinct()

target_module_edges_for_plot <- drug_target_to_module_edges %>%
  dplyr::transmute(
    from = paste0("Gene:", target_id),
    to = paste0("Gene:", module_id),
    from_type = "drug_target",
    to_type = "NiV_module_gene",
    edge_type = "target_module_ppi",
    drug_id = drug_id,
    target_id = target_id,
    module_id = module_id
  ) %>%
  dplyr::distinct()

module_internal_edges_for_plot <- module_internal_edges %>%
  dplyr::transmute(
    from = paste0("Gene:", x_id),
    to = paste0("Gene:", y_id),
    from_type = "NiV_module_gene",
    to_type = "NiV_module_gene",
    edge_type = "module_module_ppi",
    drug_id = NA_character_,
    target_id = NA_character_,
    module_id = NA_character_
  ) %>%
  dplyr::distinct()

drug_target_module_subnetwork_edges <- dplyr::bind_rows(
  drug_target_edges_for_plot,
  target_module_edges_for_plot,
  module_internal_edges_for_plot
) %>%
  dplyr::distinct()

write.csv(
  drug_target_module_subnetwork_edges,
  file.path(network_table_dir, "drug_target_module_subnetwork_edges.csv"),
  row.names = FALSE
)

# 3.6 subnetwork nodes for plotting
drug_nodes <- drug_target_module_subnetwork_edges %>%
  dplyr::filter(from_type == "drug") %>%
  dplyr::transmute(
    node_id = from,
    node_label = gsub("^Drug:", "", from),
    node_type = "drug"
  ) %>%
  dplyr::distinct()

gene_ids_in_edges <- unique(gsub(
  "^Gene:",
  "",
  c(
    drug_target_module_subnetwork_edges$from[
      grepl("^Gene:", drug_target_module_subnetwork_edges$from)
    ],
    drug_target_module_subnetwork_edges$to[
      grepl("^Gene:", drug_target_module_subnetwork_edges$to)
    ]
  )
))

gene_nodes <- data.frame(
  Entrezid = gene_ids_in_edges,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    dplyr::bind_rows(
      drug_targets_clean %>%
        dplyr::select(Entrezid, SYMBOL),
      NiV_module_genes_clean %>%
        dplyr::select(Entrezid, SYMBOL)
    ) %>%
      dplyr::distinct(Entrezid, .keep_all = TRUE),
    by = "Entrezid"
  ) %>%
  dplyr::mutate(
    node_id = paste0("Gene:", Entrezid),
    node_label = ifelse(is.na(SYMBOL) | SYMBOL == "", Entrezid, SYMBOL),
    is_drug_target = Entrezid %in% drug_targets_clean$Entrezid,
    is_NiV_module = Entrezid %in% niv_genes,
    node_type = dplyr::case_when(
      is_drug_target & is_NiV_module ~ "drug_target_and_NiV_module",
      is_drug_target ~ "drug_target",
      is_NiV_module ~ "NiV_module_gene",
      TRUE ~ "other_gene"
    )
  ) %>%
  dplyr::select(
    node_id,
    node_label,
    node_type,
    Entrezid,
    SYMBOL,
    is_drug_target,
    is_NiV_module
  ) %>%
  dplyr::distinct()

drug_target_module_subnetwork_nodes <- dplyr::bind_rows(
  drug_nodes %>%
    dplyr::mutate(
      Entrezid = NA_character_,
      SYMBOL = NA_character_,
      is_drug_target = FALSE,
      is_NiV_module = FALSE
    ) %>%
    dplyr::select(
      node_id,
      node_label,
      node_type,
      Entrezid,
      SYMBOL,
      is_drug_target,
      is_NiV_module
    ),
  gene_nodes
) %>%
  dplyr::distinct()

write.csv(
  drug_target_module_subnetwork_nodes,
  file.path(network_table_dir, "drug_target_module_subnetwork_nodes.csv"),
  row.names = FALSE
)

cat("Saved network plot tables to:", network_table_dir, "\n")

############################################################
# 4. distance（只保留 topK）
############################################################
calc_distance_topk <- function(targets, disease, g) {
  
  targets <- intersect(targets, igraph::V(g)$name)
  disease <- intersect(disease, igraph::V(g)$name)
  
  if (length(targets) == 0 || length(disease) == 0) return(NA_real_)
  
  d <- igraph::distances(g, v = targets, to = disease)
  d[d == Inf] <- NA
  
  d1 <- apply(d, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    min(x)
  })
  
  d1 <- d1[!is.na(d1)]
  if (length(d1) == 0) return(NA_real_)
  
  k <- min(5, length(d1))
  mean(sort(d1)[1:k])
}

############################################################
# 5. degree matching
# 注意：按你的要求，fallback 不做修改
############################################################
get_degree_matched <- function(targets, deg, tol = 0.2) {
  
  sampled <- character(0)
  
  for (t in targets) {
    
    if (!(t %in% names(deg))) next
    
    d <- deg[[t]]
    
    if (d <= small_deg_threshold) {
      candidates <- names(deg)[abs(deg - d) <= 2]
    } else {
      candidates <- names(deg)[abs(deg - d) / (d + 1) <= tol]
    }
    
    candidates <- setdiff(candidates, c(sampled, t))
    
    if (length(candidates) < 5) {
      candidates <- names(deg)[order(abs(deg - d))[1:50]]
    }
    
    if (length(candidates) == 0) next
    
    sampled <- c(sampled, sample(candidates, 1))
  }
  
  sampled
}

############################################################
# 6. zscore
############################################################
calc_zscore <- function(targets, disease, g, deg,
                        tol = 0.2, n = 1000) {
  
  targets <- intersect(targets, igraph::V(g)$name)
  disease <- intersect(disease, igraph::V(g)$name)
  
  if (length(targets) < 1 || length(disease) < 1) {
    return(list(distance = NA, z = NA, p = NA))
  }
  
  d_obs <- calc_distance_topk(targets, disease, g)
  if (is.na(d_obs)) {
    return(list(distance = NA, z = NA, p = NA))
  }
  
  rand <- numeric(n)
  
  for (i in seq_len(n)) {
    rand_targets <- get_degree_matched(targets, deg, tol)
    rand[i] <- calc_distance_topk(rand_targets, disease, g)
  }
  
  rand <- rand[!is.na(rand)]
  
  if (length(rand) == 0) {
    return(list(distance = d_obs, z = NA, p = NA))
  }
  
  m <- mean(rand)
  s <- sd(rand)
  
  z <- ifelse(s == 0, NA, (d_obs - m) / s)
  p <- (sum(rand <= d_obs) + 1) / (length(rand) + 1)
  
  list(distance = d_obs, z = z, p = p)
}

############################################################
# 7. 并行运行（5核）
############################################################
run_analysis <- function() {
  
  cat("Running: topk\n")
  
  cl <- makeCluster(5)
  
  ############################################################
  # 新增：固定并行随机数流，提高 permutation 可复现性
  ############################################################
  parallel::clusterSetRNGStream(cl, 123)
  
  clusterExport(cl, c(
    "drug_list", "niv_genes", "g", "deg",
    "calc_zscore", "get_degree_matched",
    "symbol_map", "small_deg_threshold",
    "calc_distance_topk"
  ), envir = environment())
  
  clusterEvalQ(cl, library(igraph))
  
  results_list <- parLapply(cl, seq_len(nrow(drug_list)), function(i) {
    
    targets <- drug_list$targets[[i]]
    
    res <- calc_zscore(targets, niv_genes, g, deg)
    
    overlap <- intersect(targets, niv_genes)
    
    overlap_symbols <- symbol_map[overlap]
    overlap_symbols <- unique(na.omit(overlap_symbols))
    
    data.frame(
      drug_id = drug_list$drug_id[i],
      method = "topk",
      distance = res$distance,
      zscore = res$z,
      p = res$p,
      target_n = length(targets),
      overlap_n = length(overlap),
      target_symbols = paste(drug_list$target_symbols[[i]], collapse = ";"),
      overlap_symbols = paste(overlap_symbols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  stopCluster(cl)
  
  results <- dplyr::bind_rows(results_list) %>%
    dplyr::arrange(zscore)
  
  write.csv(results, file.path(topk_dir, "final_results.csv"), row.names = FALSE)
  
  cat("Saved:", file.path(topk_dir, "final_results.csv"), "\n")
}

############################################################
# 8. 运行
############################################################
run_analysis()

cat("DONE\n")



############################################################
# 16_prepare_network_plot_tables_from_existing_results_FIXED.R
#
# Purpose:
#   Generate network-ready node/edge tables for plotting
#   drug-target-NiV module network from existing proximity results.
#
# Important:
#   This script DOES NOT recalculate network proximity.
#   It only prepares plotting tables.
#
# Input:
#   /users/hzhang1/project/network/data/Nipah_PPI.xlsx
#   /users/hzhang1/project/network/data/human.xlsx
#   /users/hzhang1/project/network/data/drug_target.csv
#   /users/hzhang1/project/network/result/topk_distance/final_results_with_info.csv
#   or
#   /users/hzhang1/project/network/result/topk_distance/final_results.csv
#
# Output:
#   /users/hzhang1/project/network/result/topk_distance/network_plot_tables/
############################################################

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

############################################################
# 0. Paths
############################################################
project_dir <- "/users/hzhang1/project/network"
setwd(project_dir)

data_dir <- file.path(project_dir, "data")

topk_dir <- file.path(
  project_dir,
  "result",
  "topk_distance"
)

network_table_dir <- file.path(
  topk_dir,
  "network_plot_tables"
)

dir.create(network_table_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:\n", project_dir, "\n\n")
cat("TopK result directory:\n", topk_dir, "\n\n")
cat("Network table output directory:\n", network_table_dir, "\n\n")

############################################################
# 1. User settings
############################################################
# selected_mode:
#   "zscore_lt_cutoff" : select drugs with zscore < z_cutoff
#   "top_n"            : select top N drugs with the lowest zscore
#   "manual"           : manually specify drug_id

selected_mode <- "zscore_lt_cutoff"

z_cutoff <- -1
top_n <- 20

manual_drug_ids <- c(
  # Example:
  # "DB01083", "DB12095"
)

############################################################
# 2. Read existing network proximity results
############################################################
result_with_info_file <- file.path(
  topk_dir,
  "final_results_with_info.csv"
)

result_file <- file.path(
  topk_dir,
  "final_results.csv"
)

if (file.exists(result_with_info_file)) {
  cat("Reading:\n", result_with_info_file, "\n\n")
  proximity_res <- as.data.frame(data.table::fread(result_with_info_file))
} else if (file.exists(result_file)) {
  cat("Reading:\n", result_file, "\n\n")
  proximity_res <- as.data.frame(data.table::fread(result_file))
} else {
  stop(
    "Cannot find final_results_with_info.csv or final_results.csv in: ",
    topk_dir
  )
}

required_result_cols <- c(
  "drug_id",
  "distance",
  "zscore",
  "p",
  "target_n",
  "overlap_n"
)

missing_result_cols <- setdiff(required_result_cols, colnames(proximity_res))

if (length(missing_result_cols) > 0) {
  stop(
    "Missing required columns in proximity result: ",
    paste(missing_result_cols, collapse = ", ")
  )
}

proximity_res <- proximity_res %>%
  dplyr::mutate(
    drug_id = as.character(drug_id),
    zscore = as.numeric(zscore),
    p = as.numeric(p),
    distance = as.numeric(distance),
    target_n = as.numeric(target_n),
    overlap_n = as.numeric(overlap_n)
  ) %>%
  dplyr::arrange(zscore)

if (!"name" %in% colnames(proximity_res)) {
  proximity_res$name <- proximity_res$drug_id
}

write.csv(
  proximity_res,
  file.path(network_table_dir, "network_proximity_results_used.csv"),
  row.names = FALSE
)

############################################################
# 3. Select drugs for subnetwork plotting
############################################################
if (selected_mode == "manual") {
  
  selected_drugs <- proximity_res %>%
    dplyr::filter(drug_id %in% manual_drug_ids)
  
} else if (selected_mode == "top_n") {
  
  selected_drugs <- proximity_res %>%
    dplyr::filter(!is.na(zscore)) %>%
    dplyr::arrange(zscore) %>%
    dplyr::slice_head(n = top_n)
  
} else if (selected_mode == "zscore_lt_cutoff") {
  
  selected_drugs <- proximity_res %>%
    dplyr::filter(!is.na(zscore), zscore < z_cutoff) %>%
    dplyr::arrange(zscore)
  
  if (nrow(selected_drugs) == 0) {
    selected_drugs <- proximity_res %>%
      dplyr::filter(!is.na(zscore)) %>%
      dplyr::arrange(zscore) %>%
      dplyr::slice_head(n = top_n)
  }
  
} else {
  stop("Unknown selected_mode: ", selected_mode)
}

selected_drugs <- selected_drugs %>%
  dplyr::distinct(drug_id, .keep_all = TRUE)

write.csv(
  selected_drugs,
  file.path(network_table_dir, "selected_drugs_for_subnetwork.csv"),
  row.names = FALSE
)

cat("Selected drugs for subnetwork:", nrow(selected_drugs), "\n\n")

############################################################
# 4. NiV disease module genes
############################################################
niv_file <- file.path(data_dir, "Nipah_PPI.xlsx")

if (!file.exists(niv_file)) {
  stop("Missing file: ", niv_file)
}

niv <- readxl::read_excel(niv_file)

if (!"Human_GeneID" %in% colnames(niv)) {
  stop("Nipah_PPI.xlsx must contain column: Human_GeneID")
}

niv_genes <- niv %>%
  dplyr::mutate(Human_GeneID = as.character(Human_GeneID)) %>%
  tidyr::separate_rows(Human_GeneID, sep = ";") %>%
  dplyr::mutate(Human_GeneID = trimws(Human_GeneID)) %>%
  dplyr::filter(!is.na(Human_GeneID), Human_GeneID != "") %>%
  dplyr::filter(grepl("^[0-9]+$", Human_GeneID)) %>%
  dplyr::distinct(Human_GeneID) %>%
  dplyr::pull(Human_GeneID)

niv_symbol_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(niv_genes),
  columns = "SYMBOL",
  keytype = "ENTREZID"
) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

NiV_module_genes_clean <- data.frame(
  Entrezid = niv_genes,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    niv_symbol_df,
    by = c("Entrezid" = "ENTREZID")
  ) %>%
  dplyr::mutate(
    node_type = "NiV_module_gene"
  ) %>%
  dplyr::arrange(Entrezid)

write.csv(
  NiV_module_genes_clean,
  file.path(network_table_dir, "NiV_module_genes_clean.csv"),
  row.names = FALSE
)

############################################################
# 5. Human PPI edges
############################################################
ppi_file <- file.path(data_dir, "human.xlsx")

if (!file.exists(ppi_file)) {
  stop("Missing file: ", ppi_file)
}

ppi <- readxl::read_excel(ppi_file)

required_ppi_cols <- c("x_id", "y_id")
missing_ppi_cols <- setdiff(required_ppi_cols, colnames(ppi))

if (length(missing_ppi_cols) > 0) {
  stop(
    "human.xlsx must contain columns: ",
    paste(required_ppi_cols, collapse = ", ")
  )
}

edges <- ppi %>%
  dplyr::select(x_id, y_id) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
  dplyr::mutate(
    x_id = trimws(x_id),
    y_id = trimws(y_id)
  ) %>%
  dplyr::filter(!is.na(x_id), !is.na(y_id), x_id != "", y_id != "") %>%
  dplyr::filter(x_id != y_id) %>%
  dplyr::distinct()

ppi_nodes <- unique(c(edges$x_id, edges$y_id))
niv_genes_in_ppi <- intersect(niv_genes, ppi_nodes)

NiV_module_genes_clean <- NiV_module_genes_clean %>%
  dplyr::mutate(
    in_human_ppi = Entrezid %in% niv_genes_in_ppi
  )

write.csv(
  NiV_module_genes_clean,
  file.path(network_table_dir, "NiV_module_genes_clean.csv"),
  row.names = FALSE
)

############################################################
# 6. Drug-target clean table
############################################################
drug_target_file <- file.path(data_dir, "drug_target.csv")

if (!file.exists(drug_target_file)) {
  stop("Missing file: ", drug_target_file)
}

drug <- data.table::fread(drug_target_file)

required_drug_cols <- c("drug_id", "Entrezid")
missing_drug_cols <- setdiff(required_drug_cols, colnames(drug))

if (length(missing_drug_cols) > 0) {
  stop(
    "drug_target.csv must contain columns: ",
    paste(required_drug_cols, collapse = ", ")
  )
}

drug_clean <- drug %>%
  dplyr::mutate(
    drug_id = as.character(drug_id),
    Entrezid = as.character(Entrezid)
  ) %>%
  tidyr::separate_rows(Entrezid, sep = ";") %>%
  dplyr::mutate(
    Entrezid = trimws(Entrezid)
  ) %>%
  dplyr::filter(!is.na(Entrezid), Entrezid != "") %>%
  dplyr::filter(grepl("^[0-9]+$", Entrezid))

id2symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(drug_clean$Entrezid),
  columns = "SYMBOL",
  keytype = "ENTREZID"
) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

drug_clean <- drug_clean %>%
  dplyr::left_join(
    id2symbol,
    by = c("Entrezid" = "ENTREZID")
  )

drug_targets_clean <- drug_clean %>%
  dplyr::select(drug_id, Entrezid, SYMBOL) %>%
  dplyr::distinct() %>%
  dplyr::left_join(
    proximity_res %>%
      dplyr::select(
        drug_id,
        drug_name = name,
        distance,
        zscore,
        p,
        target_n,
        overlap_n
      ) %>%
      dplyr::distinct(drug_id, .keep_all = TRUE),
    by = "drug_id"
  ) %>%
  dplyr::mutate(
    node_type = "drug_target",
    in_human_ppi = Entrezid %in% ppi_nodes,
    is_selected_drug = drug_id %in% selected_drugs$drug_id
  ) %>%
  dplyr::arrange(drug_id, Entrezid)

write.csv(
  drug_targets_clean,
  file.path(network_table_dir, "drug_targets_clean.csv"),
  row.names = FALSE
)

############################################################
# 7. NiV module internal PPI edges
############################################################
module_internal_edges <- edges %>%
  dplyr::filter(
    x_id %in% niv_genes_in_ppi,
    y_id %in% niv_genes_in_ppi
  ) %>%
  dplyr::mutate(
    edge_type = "module_module_ppi"
  ) %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(x_id = Entrezid, x_symbol = SYMBOL),
    by = "x_id"
  ) %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(y_id = Entrezid, y_symbol = SYMBOL),
    by = "y_id"
  ) %>%
  dplyr::select(
    x_id,
    x_symbol,
    y_id,
    y_symbol,
    edge_type
  ) %>%
  dplyr::distinct()

write.csv(
  module_internal_edges,
  file.path(network_table_dir, "module_internal_edges.csv"),
  row.names = FALSE
)

############################################################
# 8. Drug target to NiV module direct PPI edges
############################################################
target_genes <- unique(drug_targets_clean$Entrezid)

target_module_edges_1 <- edges %>%
  dplyr::filter(
    x_id %in% target_genes,
    y_id %in% niv_genes_in_ppi
  ) %>%
  dplyr::transmute(
    target_id = x_id,
    module_id = y_id,
    edge_type = "target_module_ppi"
  )

target_module_edges_2 <- edges %>%
  dplyr::filter(
    y_id %in% target_genes,
    x_id %in% niv_genes_in_ppi
  ) %>%
  dplyr::transmute(
    target_id = y_id,
    module_id = x_id,
    edge_type = "target_module_ppi"
  )

target_module_edges_gene <- dplyr::bind_rows(
  target_module_edges_1,
  target_module_edges_2
) %>%
  dplyr::distinct()

drug_target_to_module_edges <- target_module_edges_gene %>%
  dplyr::left_join(
    drug_targets_clean %>%
      dplyr::select(
        drug_id,
        drug_name,
        target_id = Entrezid,
        target_symbol = SYMBOL,
        distance,
        zscore,
        p,
        target_n,
        overlap_n,
        is_selected_drug
      ),
    by = "target_id"
  ) %>%
  dplyr::left_join(
    NiV_module_genes_clean %>%
      dplyr::select(
        module_id = Entrezid,
        module_symbol = SYMBOL
      ),
    by = "module_id"
  ) %>%
  dplyr::filter(!is.na(drug_id)) %>%
  dplyr::select(
    drug_id,
    drug_name,
    distance,
    zscore,
    p,
    target_n,
    overlap_n,
    is_selected_drug,
    target_id,
    target_symbol,
    module_id,
    module_symbol,
    edge_type
  ) %>%
  dplyr::arrange(zscore, drug_id, target_id, module_id) %>%
  dplyr::distinct()

write.csv(
  drug_target_to_module_edges,
  file.path(network_table_dir, "drug_target_to_module_edges.csv"),
  row.names = FALSE
)

############################################################
# 9. Build selected drug-target-NiV module subnetwork edges
############################################################
selected_target_module_edges <- drug_target_to_module_edges %>%
  dplyr::filter(drug_id %in% selected_drugs$drug_id)

selected_connected_targets <- unique(selected_target_module_edges$target_id)
selected_connected_modules <- unique(selected_target_module_edges$module_id)

drug_target_edges_for_plot <- drug_targets_clean %>%
  dplyr::filter(
    drug_id %in% selected_drugs$drug_id,
    Entrezid %in% selected_connected_targets
  ) %>%
  dplyr::transmute(
    from = paste0("Drug:", drug_id),
    to = paste0("Gene:", Entrezid),
    from_type = "drug",
    to_type = "drug_target",
    edge_type = "drug_target",
    drug_id = drug_id,
    drug_name = drug_name,
    target_id = Entrezid,
    target_symbol = SYMBOL,
    module_id = NA_character_,
    module_symbol = NA_character_,
    zscore = zscore,
    p = p
  ) %>%
  dplyr::distinct()

target_module_edges_for_plot <- selected_target_module_edges %>%
  dplyr::transmute(
    from = paste0("Gene:", target_id),
    to = paste0("Gene:", module_id),
    from_type = "drug_target",
    to_type = "NiV_module_gene",
    edge_type = "target_module_ppi",
    drug_id = drug_id,
    drug_name = drug_name,
    target_id = target_id,
    target_symbol = target_symbol,
    module_id = module_id,
    module_symbol = module_symbol,
    zscore = zscore,
    p = p
  ) %>%
  dplyr::distinct()

module_internal_edges_for_plot <- module_internal_edges %>%
  dplyr::filter(
    x_id %in% selected_connected_modules,
    y_id %in% selected_connected_modules
  ) %>%
  dplyr::transmute(
    from = paste0("Gene:", x_id),
    to = paste0("Gene:", y_id),
    from_type = "NiV_module_gene",
    to_type = "NiV_module_gene",
    edge_type = "module_module_ppi",
    drug_id = NA_character_,
    drug_name = NA_character_,
    target_id = NA_character_,
    target_symbol = NA_character_,
    module_id = NA_character_,
    module_symbol = NA_character_,
    zscore = NA_real_,
    p = NA_real_
  ) %>%
  dplyr::distinct()

drug_target_module_subnetwork_edges <- dplyr::bind_rows(
  drug_target_edges_for_plot,
  target_module_edges_for_plot,
  module_internal_edges_for_plot
) %>%
  dplyr::distinct()

write.csv(
  drug_target_module_subnetwork_edges,
  file.path(network_table_dir, "drug_target_module_subnetwork_edges.csv"),
  row.names = FALSE
)

############################################################
# 10. Build selected drug-target-NiV module subnetwork nodes
############################################################
drug_nodes <- selected_drugs %>%
  dplyr::filter(drug_id %in% unique(selected_target_module_edges$drug_id)) %>%
  dplyr::transmute(
    node_id = paste0("Drug:", drug_id),
    node_label = ifelse(is.na(name) | name == "", drug_id, name),
    node_type = "drug",
    drug_id = drug_id,
    drug_name = name,
    Entrezid = NA_character_,
    SYMBOL = NA_character_,
    zscore = zscore,
    p = p,
    is_drug_target = FALSE,
    is_NiV_module = FALSE
  ) %>%
  dplyr::distinct()

gene_ids_in_edges <- unique(gsub(
  "^Gene:",
  "",
  c(
    drug_target_module_subnetwork_edges$from[
      grepl("^Gene:", drug_target_module_subnetwork_edges$from)
    ],
    drug_target_module_subnetwork_edges$to[
      grepl("^Gene:", drug_target_module_subnetwork_edges$to)
    ]
  )
))

gene_symbol_source <- dplyr::bind_rows(
  drug_targets_clean %>%
    dplyr::select(Entrezid, SYMBOL),
  NiV_module_genes_clean %>%
    dplyr::select(Entrezid, SYMBOL)
) %>%
  dplyr::distinct(Entrezid, .keep_all = TRUE)

gene_nodes <- data.frame(
  Entrezid = gene_ids_in_edges,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_symbol_source, by = "Entrezid") %>%
  dplyr::mutate(
    node_id = paste0("Gene:", Entrezid),
    node_label = ifelse(is.na(SYMBOL) | SYMBOL == "", Entrezid, SYMBOL),
    is_drug_target = Entrezid %in% selected_connected_targets,
    is_NiV_module = Entrezid %in% niv_genes_in_ppi,
    node_type = dplyr::case_when(
      is_drug_target & is_NiV_module ~ "drug_target_and_NiV_module",
      is_drug_target ~ "drug_target",
      is_NiV_module ~ "NiV_module_gene",
      TRUE ~ "other_gene"
    ),
    drug_id = NA_character_,
    drug_name = NA_character_,
    zscore = NA_real_,
    p = NA_real_
  ) %>%
  dplyr::select(
    node_id,
    node_label,
    node_type,
    drug_id,
    drug_name,
    Entrezid,
    SYMBOL,
    zscore,
    p,
    is_drug_target,
    is_NiV_module
  ) %>%
  dplyr::distinct()

drug_target_module_subnetwork_nodes <- dplyr::bind_rows(
  drug_nodes,
  gene_nodes
) %>%
  dplyr::distinct()

write.csv(
  drug_target_module_subnetwork_nodes,
  file.path(network_table_dir, "drug_target_module_subnetwork_nodes.csv"),
  row.names = FALSE
)

############################################################
# 11. QC summary
############################################################
qc_lines <- c(
  "Network plot table preparation finished.",
  "",
  paste0("Selected mode: ", selected_mode),
  paste0("Z-score cutoff: ", z_cutoff),
  paste0("Top N: ", top_n),
  "",
  paste0("NiV module genes raw: ", length(niv_genes)),
  paste0("NiV module genes in PPI: ", length(niv_genes_in_ppi)),
  paste0("Drug targets clean rows: ", nrow(drug_targets_clean)),
  paste0("Unique drug target genes: ", length(unique(drug_targets_clean$Entrezid))),
  paste0("Module internal PPI edges: ", nrow(module_internal_edges)),
  paste0("Drug-target-to-module direct PPI edges: ", nrow(drug_target_to_module_edges)),
  "",
  paste0("Selected drugs: ", nrow(selected_drugs)),
  paste0(
    "Selected drugs with direct target-module PPI edges: ",
    length(unique(selected_target_module_edges$drug_id))
  ),
  paste0("Subnetwork edges: ", nrow(drug_target_module_subnetwork_edges)),
  paste0("Subnetwork nodes: ", nrow(drug_target_module_subnetwork_nodes))
)

writeLines(
  qc_lines,
  con = file.path(network_table_dir, "network_plot_table_QC_summary.txt"),
  useBytes = TRUE
)

############################################################
# 12. Finish
############################################################
cat("\n============================================================\n")
cat("Network plot tables generated successfully.\n\n")

cat("Input result directory:\n")
cat(normalizePath(topk_dir), "\n\n")

cat("Output directory:\n")
cat(normalizePath(network_table_dir), "\n\n")

cat("Generated files:\n")
cat("1) NiV_module_genes_clean.csv\n")
cat("2) drug_targets_clean.csv\n")
cat("3) module_internal_edges.csv\n")
cat("4) drug_target_to_module_edges.csv\n")
cat("5) drug_target_module_subnetwork_edges.csv\n")
cat("6) drug_target_module_subnetwork_nodes.csv\n")
cat("7) selected_drugs_for_subnetwork.csv\n")
cat("8) network_plot_table_QC_summary.txt\n")
cat("============================================================\n")