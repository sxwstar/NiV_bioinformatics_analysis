# 02_temporal_analysis_and_cluster_enrichment_analysis.R
#
# 
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# 0) Paths
############################################################
work_dir <- "D:/WeiRuan/NiV/GSE166707/results_FINAL_timepoint_FULL_SCI/08_maSigPro"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

expr_file <- "D:/WeiRuan/NiV/GSE166707/results_FINAL_timepoint_FULL_SCI/00_QC/expr_gene_matrix.csv"
sample_file <- "D:/WeiRuan/NiV/GSE166707/results_FINAL_timepoint_FULL_SCI/00_QC/sample_info_all.csv"

subdirs <- c(
  "00_input_check",
  "01_maSigPro_tables",
  "02_cluster_results",
  "03_cluster_plots",
  "04_enrichment_results/GO",
  "04_enrichment_results/KEGG",
  "04_enrichment_results/summary",
  "05_enrichment_plots/GO",
  "05_enrichment_plots/KEGG",
  "05_enrichment_plots/selected_summary",
  "06_combined_plots"
)
invisible(lapply(file.path(work_dir, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE))

cat("============================================================\n")
cat("NiV maSigPro temporal analysis and cluster enrichment START\n")
cat("Final temporal cluster number: k = 6\n")
cat("Working directory:", work_dir, "\n")
cat("Expression file:", expr_file, "\n")
cat("Sample info file:", sample_file, "\n")
cat("============================================================\n\n")

############################################################
# 1) Packages
############################################################
need_pkgs <- c(
  "maSigPro", "limma", "ggplot2", "dplyr", "tidyr", "tibble",
  "stringr", "pheatmap", "RColorBrewer", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "patchwork"
)

install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(miss) == 0) return(invisible(TRUE))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc_avail <- BiocManager::available()
  bioc_to_install <- intersect(miss, bioc_avail)
  cran_to_install <- setdiff(miss, bioc_to_install)
  if (length(bioc_to_install) > 0) BiocManager::install(bioc_to_install, ask = FALSE, update = FALSE)
  if (length(cran_to_install) > 0) install.packages(cran_to_install)
  invisible(TRUE)
}
install_if_missing(need_pkgs)

suppressPackageStartupMessages({
  library(maSigPro)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(patchwork)
})

############################################################
# 2) Fixed parameters
############################################################
degree_use <- 2
maSigPro_Q <- 0.05
tfit_alpha <- 0.05
rsq_cutoff <- 0.6

# Final cluster number used in the manuscript
cluster_k <- 6

mean_cutoff <- 1
sd_cutoff <- 0.2
do_quantile_normalization <- FALSE

group_levels <- c("Uninfected", "Infection")
time_levels <- c("4h", "8h", "12h", "16h")

group_cols <- c("Uninfected" = "#4C72B0", "Infection" = "#D94F3D")
time_cols <- c("4h" = "#1B9E77", "8h" = "#D95F02", "12h" = "#7570B3", "16h" = "#E7298A")

############################################################
# 3) Helper functions
############################################################
theme_sci <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1, color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey96", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", color = "black")
    )
}

write_table_dual <- function(x, file_csv, row.names = FALSE, quote = FALSE, sep_txt = "\t", na = "") {
  utils::write.csv(x, file_csv, row.names = row.names, quote = quote, na = na)
  file_txt <- sub("\\.csv$", ".txt", file_csv, ignore.case = TRUE)
  utils::write.table(
    x,
    file = file_txt,
    sep = sep_txt,
    row.names = row.names,
    col.names = TRUE,
    quote = FALSE,
    na = na
  )
  invisible(TRUE)
}

safe_save_plot <- function(p, filename, width = 8, height = 6, dpi = 320) {
  ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = dpi, bg = "white")
}

zscore_rows <- function(mat) {
  mat2 <- t(scale(t(mat)))
  mat2[!is.finite(mat2)] <- 0
  mat2
}

convert_ids <- function(genes) {
  genes <- unique(genes)
  ids <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
  ids <- unique(ids[!is.na(ids)])
  ids
}

parse_gene_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) != 2) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  })
}

wrap_term <- function(x, width = 38) {
  stringr::str_wrap(x, width = width)
}

collect_tables_recursive <- function(x, path = "root", max_depth = 6, depth = 1) {
  out <- list()
  if (depth > max_depth || is.null(x)) return(out)

  if (is.data.frame(x) || is.matrix(x)) {
    out[[path]] <- as.data.frame(x)
    return(out)
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      child_path <- paste0(path, "/", nms[i])
      out <- c(out, collect_tables_recursive(x[[i]], child_path, max_depth, depth + 1))
    }
  }
  out
}

collect_raw_candidates_recursive <- function(x, path = "root", max_depth = 6, depth = 1) {
  out <- list()
  if (depth > max_depth || is.null(x)) return(out)

  if (is.list(x) || is.data.frame(x) || is.matrix(x)) {
    out[[path]] <- x
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) {
      child_path <- paste0(path, "/", nms[i])
      out <- c(out, collect_raw_candidates_recursive(x[[i]], child_path, max_depth, depth + 1))
    }
  }
  out
}

pick_export_sig_table <- function(sigs_obj, fallback_genes = NULL) {
  tabs <- collect_tables_recursive(sigs_obj, path = "sigs")
  if (length(tabs) == 0) return(NULL)

  score_df <- data.frame(
    path = names(tabs),
    nrow = vapply(tabs, nrow, integer(1)),
    ncol = vapply(tabs, ncol, integer(1)),
    overlap = 0L,
    stringsAsFactors = FALSE
  )

  if (!is.null(fallback_genes)) {
    score_df$overlap <- vapply(tabs, function(tb) {
      rn <- rownames(tb)
      if (is.null(rn)) return(0L)
      length(intersect(rn, fallback_genes))
    }, integer(1))
  }

  score_df <- score_df[score_df$nrow > 0 & score_df$ncol >= 1, , drop = FALSE]
  if (nrow(score_df) == 0) return(NULL)

  score_df <- score_df[order(-score_df$overlap, -score_df$nrow, -score_df$ncol), , drop = FALSE]
  best_path <- score_df$path[1]
  best_tab <- tabs[[best_path]]
  attr(best_tab, "selected_path") <- best_path
  best_tab
}

find_working_raw_for_see_genes <- function(sigs_obj, dis, cluster_k = 6) {
  candidates <- collect_raw_candidates_recursive(sigs_obj, path = "sigs")
  if (length(candidates) == 0) return(NULL)

  cand_names <- names(candidates)
  log_df <- data.frame(
    candidate = cand_names,
    ok = FALSE,
    n_genes = NA_integer_,
    message = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(candidates)) {
    obj <- candidates[[i]]
    tmp_pdf <- tempfile(fileext = ".pdf")

    res <- tryCatch({
      pdf(tmp_pdf, width = 6, height = 4)
      tmp <- see.genes(
        obj,
        show.fit = FALSE,
        dis = dis,
        cluster.method = "hclust",
        cluster.data = 1,
        k = cluster_k
      )
      dev.off()

      if (!is.null(tmp$cut) && length(tmp$cut) > 0) {
        list(ok = TRUE, n_genes = length(tmp$cut), msg = "success", result = tmp)
      } else {
        list(ok = FALSE, n_genes = NA_integer_, msg = "empty cut", result = NULL)
      }
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      list(ok = FALSE, n_genes = NA_integer_, msg = conditionMessage(e), result = NULL)
    })

    log_df$ok[i] <- res$ok
    log_df$n_genes[i] <- res$n_genes
    log_df$message[i] <- res$msg

    if (isTRUE(res$ok)) {
      attr(obj, "selected_path") <- cand_names[i]
      attr(obj, "trial_log") <- log_df
      return(obj)
    }
  }

  attr(candidates, "trial_log") <- log_df
  NULL
}

############################################################
# 4) Read inputs
############################################################
if (!file.exists(expr_file)) stop("Missing expr file: ", expr_file)
if (!file.exists(sample_file)) stop("Missing sample file: ", sample_file)

expr_df <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
sample_info <- read.csv(sample_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "group", "time") %in% colnames(sample_info))) {
  stop("sample_info_all.csv must contain columns: sample, group, time")
}

expr0 <- as.matrix(expr_df)
storage.mode(expr0) <- "numeric"
stopifnot(!anyNA(expr0))

sample_info <- sample_info %>%
  mutate(
    group = factor(group, levels = group_levels),
    time = factor(time, levels = time_levels)
  ) %>%
  as.data.frame()

common_samples <- intersect(sample_info$sample, colnames(expr0))
if (length(common_samples) == 0) stop("No overlapping samples between expr matrix and sample_info.")

sample_info <- sample_info[match(common_samples, sample_info$sample), , drop = FALSE]
expr0 <- expr0[, common_samples, drop = FALSE]
rownames(sample_info) <- sample_info$sample

write_table_dual(
  sample_info,
  file.path(work_dir, "00_input_check", "sample_info_used.csv"),
  row.names = FALSE
)

cat("Expression matrix dim:", nrow(expr0), "genes x", ncol(expr0), "samples\n")

############################################################
# 5) Preprocessing
############################################################
expr_range <- range(expr0, na.rm = TRUE)
cat("Expression range before transform:", paste(round(expr_range, 3), collapse = " ~ "), "\n")

if (max(expr0, na.rm = TRUE) > 50) {
  expr0 <- log2(expr0 + 1)
  cat("Applied log2(x+1) transform.\n")
}

if (do_quantile_normalization) {
  expr0 <- normalizeBetweenArrays(expr0, method = "quantile")
  cat("Applied quantile normalization.\n")
}

expr <- expr0
expr <- expr[rowMeans(expr) > mean_cutoff, , drop = FALSE]
expr <- expr[apply(expr, 1, sd) > sd_cutoff, , drop = FALSE]

writeLines(
  c(
    paste0("Original genes: ", nrow(expr0)),
    paste0("After filtering: ", nrow(expr))
  ),
  file.path(work_dir, "00_input_check", "filtering_summary.txt")
)

############################################################
# 6) Build maSigPro design
############################################################
time_num <- as.numeric(sub("h", "", as.character(sample_info$time)))

Replicate <- ave(
  seq_along(time_num),
  interaction(sample_info$group, sample_info$time, drop = TRUE),
  FUN = seq_along
)

edesign <- data.frame(
  Time = time_num,
  Replicate = Replicate,
  check.names = FALSE
)

group_factor <- factor(sample_info$group, levels = group_levels)
group_dummies <- model.matrix(~ 0 + group_factor)
colnames(group_dummies) <- sub("^group_factor", "", colnames(group_dummies))

edesign <- cbind(edesign, group_dummies)
rownames(edesign) <- sample_info$sample

write_table_dual(
  edesign,
  file.path(work_dir, "00_input_check", "maSigPro_edesign.csv"),
  row.names = TRUE
)

design <- make.design.matrix(edesign, degree = degree_use)

write_table_dual(
  design$dis,
  file.path(work_dir, "00_input_check", "maSigPro_design_matrix.csv"),
  row.names = TRUE
)

############################################################
# 7) Run maSigPro
############################################################
cat("=== Running p.vector ===\n")
fit <- p.vector(
  expr,
  design,
  Q = maSigPro_Q,
  MT.adjust = "BH"
)

if (is.null(fit$SELEC) || nrow(fit$SELEC) == 0) {
  stop("No significant genes identified by p.vector.")
}

write_table_dual(
  fit$SELEC,
  file.path(work_dir, "01_maSigPro_tables", "pvector_selected_genes.csv"),
  row.names = TRUE
)

cat("=== Running T.fit ===\n")
tstep <- T.fit(
  fit,
  step.method = "backward",
  alfa = tfit_alpha
)

cat("=== Running get.siggenes ===\n")
sigs <- get.siggenes(
  tstep,
  rsq = rsq_cutoff,
  vars = "groups"
)

capture.output({
  cat("===== names(sigs) =====\n")
  print(names(sigs))
  cat("===== str(sigs, max.level = 3) =====\n")
  str(sigs, max.level = 3)
}, file = file.path(work_dir, "01_maSigPro_tables", "get_siggenes_structure.txt"))

############################################################
# 8) Stable export table for downstream analysis
############################################################
sig_table <- pick_export_sig_table(
  sigs_obj = sigs,
  fallback_genes = rownames(fit$SELEC)
)

if (is.null(sig_table) || nrow(sig_table) == 0) {
  stop("Could not extract a stable sig table from get.siggenes(). Check get_siggenes_structure.txt")
}

selected_export_path <- attr(sig_table, "selected_path")
if (is.null(selected_export_path)) selected_export_path <- "unknown"

writeLines(
  paste0("Selected export sig table path: ", selected_export_path),
  file.path(work_dir, "01_maSigPro_tables", "selected_sig_component.txt")
)

if (is.null(rownames(sig_table)) || all(rownames(sig_table) == seq_len(nrow(sig_table)))) {
  if (nrow(sig_table) == nrow(fit$SELEC)) {
    rownames(sig_table) <- rownames(fit$SELEC)
  } else {
    stop("Selected sig table has no reliable gene rownames.")
  }
}

write_table_dual(
  as.data.frame(sig_table),
  file.path(work_dir, "01_maSigPro_tables", "siggenes_primary_component.csv"),
  row.names = TRUE
)

############################################################
# 9) Raw object for see.genes() clustering
############################################################
cat("=== Finding working raw object for see.genes() ===\n")

sig_component_raw <- find_working_raw_for_see_genes(
  sigs_obj = sigs,
  dis = design$dis,
  cluster_k = cluster_k
)

if (is.null(sig_component_raw)) {
  candidates <- collect_raw_candidates_recursive(sigs, path = "sigs")
  write.csv(
    data.frame(candidate = names(candidates), stringsAsFactors = FALSE),
    file.path(work_dir, "01_maSigPro_tables", "see_genes_candidate_paths.csv"),
    row.names = FALSE
  )
  stop("Could not find a working raw object for see.genes(). Check get_siggenes_structure.txt")
}

selected_raw_path <- attr(sig_component_raw, "selected_path")
if (is.null(selected_raw_path)) selected_raw_path <- "unknown"

writeLines(
  paste0("Selected raw object for see.genes(): ", selected_raw_path),
  file.path(work_dir, "01_maSigPro_tables", "selected_raw_component_for_see_genes.txt")
)

trial_log <- attr(sig_component_raw, "trial_log")
if (!is.null(trial_log)) {
  write.csv(
    trial_log,
    file.path(work_dir, "01_maSigPro_tables", "see_genes_trial_log.csv"),
    row.names = FALSE
  )
}

############################################################
# 10) Final k = 6 temporal clustering
############################################################
cat("=== Running see.genes() for final k = 6 temporal clusters ===\n")

pdf(file.path(work_dir, "02_cluster_results", "sig_genes_clusters_k6_raw.pdf"), width = 15, height = 12)
see_result <- see.genes(
  sig_component_raw,
  show.fit = TRUE,
  dis = design$dis,
  cluster.method = "hclust",
  cluster.data = 1,
  k = cluster_k
)
dev.off()

if (is.null(see_result$cut) || length(see_result$cut) == 0) {
  stop("see.genes() returned empty cut vector.")
}

cluster_gene_df <- data.frame(
  gene = names(see_result$cut),
  cluster = as.integer(see_result$cut),
  Cluster = as.integer(see_result$cut),
  Cluster_label = paste0("Cluster_", as.integer(see_result$cut)),
  stringsAsFactors = FALSE
) %>%
  arrange(cluster, gene)

write_table_dual(
  cluster_gene_df,
  file.path(work_dir, "02_cluster_results", "cluster_gene_mapping.csv"),
  row.names = FALSE
)
write_table_dual(
  cluster_gene_df,
  file.path(work_dir, "02_cluster_results", "cluster_gene_mapping_k6.csv"),
  row.names = FALSE
)

cluster_gene_list <- split(cluster_gene_df$gene, cluster_gene_df$cluster)

cluster_size_df <- data.frame(
  cluster = names(cluster_gene_list),
  Cluster = as.integer(names(cluster_gene_list)),
  Cluster_label = paste0("Cluster_", names(cluster_gene_list)),
  n_genes = sapply(cluster_gene_list, length),
  stringsAsFactors = FALSE
) %>%
  arrange(as.numeric(cluster))

write_table_dual(
  cluster_size_df,
  file.path(work_dir, "02_cluster_results", "cluster_gene_counts.csv"),
  row.names = FALSE
)
write_table_dual(
  cluster_size_df,
  file.path(work_dir, "02_cluster_results", "cluster_gene_counts_k6.csv"),
  row.names = FALSE
)

############################################################
# 11) Export downstream gene list
############################################################
downstream_genes <- unique(rownames(as.data.frame(sig_table)))
write_table_dual(
  data.frame(gene_symbol = downstream_genes),
  file.path(work_dir, "01_maSigPro_tables", "maSigPro_downstream_gene_list.csv"),
  row.names = FALSE
)

############################################################
# 12) Cluster size plot
############################################################
p_cluster_size <- ggplot(cluster_size_df, aes(x = factor(Cluster, levels = 1:cluster_k), y = n_genes)) +
  geom_col(width = 0.72, fill = "#6A9BCB", color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_genes), vjust = -0.25, size = 3.6, fontface = "bold") +
  labs(
    title = "Gene counts across maSigPro temporal clusters",
    x = "Cluster",
    y = "Gene count"
  ) +
  coord_cartesian(ylim = c(0, max(cluster_size_df$n_genes) * 1.12)) +
  theme_sci(12)

safe_save_plot(
  p_cluster_size,
  file.path(work_dir, "03_cluster_plots", "Cluster_size_barplot_k6_SCI.png"),
  width = 6.8, height = 5.2, dpi = 320
)
safe_save_plot(
  p_cluster_size,
  file.path(work_dir, "03_cluster_plots", "Cluster_size_barplot_k6_SCI.pdf"),
  width = 6.8, height = 5.2, dpi = 320
)

############################################################
# 13) SCI-style cluster trajectory plot
############################################################
sig_genes_plot <- unique(cluster_gene_df$gene)
expr_sig <- expr[sig_genes_plot, rownames(sample_info), drop = FALSE]
expr_sig_z <- zscore_rows(expr_sig)

sample_info_plot <- sample_info[, c("sample", "group", "time"), drop = FALSE]

expr_long <- as.data.frame(expr_sig_z) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "z_expr") %>%
  dplyr::left_join(sample_info_plot, by = "sample") %>%
  dplyr::left_join(cluster_gene_df, by = "gene") %>%
  dplyr::mutate(
    cluster = factor(cluster, levels = sort(unique(cluster_gene_df$cluster))),
    time = factor(time, levels = time_levels),
    time_num = as.numeric(sub("h", "", as.character(time))),
    group = factor(group, levels = group_levels)
  )

cluster_summary <- expr_long %>%
  dplyr::group_by(cluster, group, time, time_num) %>%
  dplyr::summarise(
    mean_z = mean(z_expr, na.rm = TRUE),
    se_z = sd(z_expr, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

p_cluster <- ggplot(cluster_summary, aes(x = time_num, y = mean_z, color = group, group = group)) +
  geom_ribbon(
    aes(ymin = mean_z - se_z, ymax = mean_z + se_z, fill = group),
    alpha = 0.14, colour = NA
  ) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  scale_color_manual(
    values = group_cols,
    breaks = c("Uninfected", "Infection"),
    labels = c("Control", "Infection")
  ) +
  scale_fill_manual(
    values = group_cols,
    breaks = c("Uninfected", "Infection"),
    labels = c("Control", "Infection")
  ) +
  scale_x_continuous(breaks = c(4, 8, 12, 16), labels = c("4h", "8h", "12h", "16h")) +
  facet_wrap(~ cluster, ncol = 3, scales = "free_y") +
  labs(
    title = "Temporal expression trajectories of maSigPro clusters",
    x = "Time post infection",
    y = "Mean scaled expression",
    color = NULL,
    fill = NULL
  ) +
  theme_sci(12) +
  theme(legend.position = "top")

safe_save_plot(
  p_cluster,
  file.path(work_dir, "03_cluster_plots", "Cluster_trajectory_plot_k6_SCI.png"),
  width = 11.5, height = 8.8, dpi = 320
)
safe_save_plot(
  p_cluster,
  file.path(work_dir, "03_cluster_plots", "Cluster_trajectory_plot_k6_SCI.pdf"),
  width = 11.5, height = 8.8, dpi = 320
)

############################################################
# 14) SCI-style heatmap
############################################################
sample_order_plot <- sample_info %>%
  mutate(time_num = as.numeric(sub("h", "", as.character(time)))) %>%
  arrange(time_num, group, sample) %>%
  pull(sample)

cluster_order <- cluster_gene_df %>%
  arrange(cluster, gene) %>%
  pull(gene)

hm_mat <- expr_sig_z[cluster_order, sample_order_plot, drop = FALSE]

ann_col <- data.frame(
  Time = sample_info[sample_order_plot, "time"],
  Group = sample_info[sample_order_plot, "group"]
)
rownames(ann_col) <- sample_order_plot

ann_colors <- list(
  Time = time_cols,
  Group = group_cols
)

gaps_row <- cumsum(as.numeric(table(cluster_gene_df$cluster)))
gaps_row <- gaps_row[-length(gaps_row)]

png(file.path(work_dir, "03_cluster_plots", "Cluster_heatmap_k6_SCI.png"),
    width = 2200, height = 2600, res = 320)
pheatmap(
  hm_mat,
  color = colorRampPalette(c("#2C5AA0", "white", "#C73E4D"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  main = "maSigPro significant genes ordered by temporal clusters",
  fontsize = 10
)
dev.off()

pdf(file.path(work_dir, "03_cluster_plots", "Cluster_heatmap_k6_SCI.pdf"),
    width = 8.5, height = 10)
pheatmap(
  hm_mat,
  color = colorRampPalette(c("#2C5AA0", "white", "#C73E4D"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  main = "maSigPro significant genes ordered by temporal clusters",
  fontsize = 10
)
dev.off()

############################################################
# 15) Enrichment analysis by cluster
############################################################
cat("=== Enrichment analysis for each k = 6 temporal cluster ===\n")

id_stats <- data.frame()
gene_entrez_list <- list()
go_summary <- data.frame()
kegg_summary <- data.frame()

for (cid in names(cluster_gene_list)) {
  genes_symbol <- unique(cluster_gene_list[[cid]])
  genes_entrez <- convert_ids(genes_symbol)
  gene_entrez_list[[cid]] <- genes_entrez

  id_stats <- rbind(
    id_stats,
    data.frame(
      Cluster = as.integer(cid),
      Cluster_label = paste0("Cluster_", cid),
      symbol_input = length(genes_symbol),
      entrez_mapped = length(genes_entrez),
      map_rate = round(length(genes_entrez) / length(genes_symbol), 3),
      stringsAsFactors = FALSE
    )
  )

  if (length(genes_entrez) < 3) {
    writeLines(
      paste0("Cluster ", cid, ": fewer than 3 mapped ENTREZ genes; enrichment skipped."),
      file.path(work_dir, "04_enrichment_results", "summary", paste0("cluster_", cid, "_SKIPPED_TOO_FEW_GENES.txt"))
    )
    next
  }

  ##########################################################
  # GO enrichment
  ##########################################################
  ego <- tryCatch(
    enrichGO(
      gene = genes_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      minGSSize = 2,
      readable = TRUE
    ),
    error = function(e) {
      warning("GO enrichment failed for cluster ", cid, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_simplified <- tryCatch(
      simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) ego
    )
    go_df <- as.data.frame(ego_simplified)

    if (nrow(go_df) > 0) {
      go_df$Cluster <- as.integer(cid)
      go_df$Cluster_label <- paste0("Cluster_", cid)
      go_df$cluster_id <- cid
      if ("GeneRatio" %in% colnames(go_df)) {
        go_df$GeneRatio_num <- parse_gene_ratio(go_df$GeneRatio)
      }

      go_summary <- rbind(go_summary, go_df)

      write_table_dual(
        go_df,
        file.path(work_dir, "04_enrichment_results", "GO", paste0("cluster_", cid, "_GO.csv")),
        row.names = FALSE
      )

      p_go <- dotplot(ego_simplified, showCategory = 18, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY ~ ., scales = "free_y") +
        theme_sci(11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste0("Cluster ", cid, " GO enrichment"))

      safe_save_plot(
        p_go,
        file.path(work_dir, "05_enrichment_plots", "GO", paste0("cluster_", cid, "_GO_dotplot.png")),
        width = 11, height = 9, dpi = 320
      )
      safe_save_plot(
        p_go,
        file.path(work_dir, "05_enrichment_plots", "GO", paste0("cluster_", cid, "_GO_dotplot.pdf")),
        width = 11, height = 9, dpi = 320
      )
    }
  }

  ##########################################################
  # KEGG enrichment
  ##########################################################
  ekegg <- tryCatch(
    enrichKEGG(
      gene = genes_entrez,
      organism = "hsa",
      pvalueCutoff = 0.1,
      pAdjustMethod = "BH",
      minGSSize = 2
    ),
    error = function(e) {
      warning("KEGG enrichment failed for cluster ", cid, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    kegg_df <- as.data.frame(ekegg)

    if (nrow(kegg_df) > 0) {
      kegg_df$Cluster <- as.integer(cid)
      kegg_df$Cluster_label <- paste0("Cluster_", cid)
      kegg_df$cluster_id <- cid
      if ("GeneRatio" %in% colnames(kegg_df)) {
        kegg_df$GeneRatio_num <- parse_gene_ratio(kegg_df$GeneRatio)
      }

      kegg_summary <- rbind(kegg_summary, kegg_df)

      write_table_dual(
        kegg_df,
        file.path(work_dir, "04_enrichment_results", "KEGG", paste0("cluster_", cid, "_KEGG.csv")),
        row.names = FALSE
      )

      p_kegg <- dotplot(ekegg, showCategory = 18) +
        theme_sci(11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste0("Cluster ", cid, " KEGG enrichment"))

      safe_save_plot(
        p_kegg,
        file.path(work_dir, "05_enrichment_plots", "KEGG", paste0("cluster_", cid, "_KEGG_dotplot.png")),
        width = 10, height = 7.5, dpi = 320
      )
      safe_save_plot(
        p_kegg,
        file.path(work_dir, "05_enrichment_plots", "KEGG", paste0("cluster_", cid, "_KEGG_dotplot.pdf")),
        width = 10, height = 7.5, dpi = 320
      )
    }
  }
}

write_table_dual(
  id_stats,
  file.path(work_dir, "04_enrichment_results", "summary", "All_clusters_ID_mapping_summary.csv"),
  row.names = FALSE
)

write_table_dual(
  go_summary,
  file.path(work_dir, "04_enrichment_results", "summary", "All_clusters_GO_summary.csv"),
  row.names = FALSE
)

write_table_dual(
  kegg_summary,
  file.path(work_dir, "04_enrichment_results", "summary", "All_clusters_KEGG_summary.csv"),
  row.names = FALSE
)

############################################################
# 16) Final selected GO / KEGG summary plots
############################################################
cat("=== Plotting selected GO / KEGG cluster summary bubbles ===\n")

selected_plot_dir <- file.path(work_dir, "05_enrichment_plots", "selected_summary")
dir.create(selected_plot_dir, recursive = TRUE, showWarnings = FALSE)

p_cutoff <- 0.05
cluster_levels <- paste0("Cluster_", 1:cluster_k)
cluster_labels_clean <- paste0("Cluster", 1:cluster_k)

fig_width <- 9.8
fig_height <- 7.4

clean_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pretty_cluster <- function(x) {
  x <- as.character(x)
  ifelse(str_detect(x, "^Cluster_"), x, paste0("Cluster_", x))
}

cluster_clean_label <- function(x) {
  x <- pretty_cluster(x)
  str_replace(x, "Cluster_", "Cluster")
}

standardize_term <- function(x) {
  x <- str_squish(as.character(x))
  x <- str_replace_all(x, "–", "-")
  x <- str_replace_all(x, "—", "-")
  x
}

pretty_term <- function(x) {
  x0 <- str_to_lower(str_squish(as.character(x)))

  case_when(
    x0 == "response to oxidative stress" ~ "Response to oxidative stress",
    x0 == "cellular response to oxidative stress" ~ "Cellular response to oxidative stress",
    x0 == "2-oxoglutarate metabolic process" ~ "2-Oxoglutarate metabolic process",
    x0 == "2-oxocarboxylic acid metabolic process" ~ "2-Oxocarboxylic acid metabolic process",

    x0 == "pattern specification process" ~ "Pattern specification process",
    x0 == "regionalization" ~ "Regionalization",
    x0 == "anterior/posterior pattern specification" ~ "Anterior/posterior pattern specification",

    x0 == "myeloid leukocyte activation" ~ "Myeloid leukocyte activation",
    x0 == "regulation of presynaptic membrane potential" ~ "Regulation of presynaptic membrane potential",
    x0 == "regulation of myeloid leukocyte differentiation" ~ "Regulation of myeloid leukocyte differentiation",

    x0 == "dna replication" ~ "DNA replication",
    x0 == "dna-templated dna replication" ~ "DNA-templated DNA replication",
    x0 == "dna strand elongation involved in dna replication" ~ "DNA strand elongation in DNA replication",

    x0 == "response to endoplasmic reticulum stress" ~ "Response to ER stress",
    x0 == "endoplasmic reticulum unfolded protein response" ~ "ER unfolded protein response",
    x0 == "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress" ~ "Intrinsic apoptotic signaling pathway in response to ER stress",

    x0 == "protein folding in endoplasmic reticulum" ~ "Protein folding in ER",
    x0 == "nucleocytoplasmic transport" ~ "Nucleocytoplasmic transport",
    x0 == "nuclear transport" ~ "Nuclear transport",
    x0 == "organelle fission" ~ "Organelle fission",

    x0 == "2-oxocarboxylic acid metabolism" ~ "2-Oxocarboxylic acid metabolism",
    x0 == "cholesterol metabolism" ~ "Cholesterol metabolism",
    x0 == "retinol metabolism" ~ "Retinol metabolism",

    x0 == "proteoglycans in cancer" ~ "Proteoglycans in cancer",
    x0 == "cytoskeleton in muscle cells" ~ "Cytoskeleton in muscle cells",
    x0 == "focal adhesion" ~ "Focal adhesion",

    x0 == "neuroactive ligand signaling" ~ "Neuroactive ligand signaling",
    x0 == "synaptic vesicle cycle" ~ "Synaptic vesicle cycle",
    x0 == "neuroactive ligand-receptor interaction" ~ "Neuroactive ligand-receptor interaction",

    x0 == "mismatch repair" ~ "Mismatch repair",
    x0 == "cell cycle" ~ "Cell cycle",

    x0 == "il-17 signaling pathway" ~ "IL-17 signaling pathway",
    x0 == "apoptosis" ~ "Apoptosis",
    x0 == "apoptosis - multiple species" ~ "Apoptosis - multiple species",

    x0 == "human immunodeficiency virus 1 infection" ~ "HIV-1 infection",
    x0 == "foxo signaling pathway" ~ "FoxO signaling pathway",
    x0 == "human t-cell leukemia virus 1 infection" ~ "HTLV-1 infection",
    x0 == "protein processing in endoplasmic reticulum" ~ "Protein processing in ER",

    TRUE ~ str_to_sentence(x0)
  )
}

selected_go_terms <- data.frame(
  Cluster_label = c(
    rep("Cluster_1", 3),
    rep("Cluster_2", 3),
    rep("Cluster_3", 3),
    rep("Cluster_4", 3),
    rep("Cluster_5", 3),
    rep("Cluster_6", 4)
  ),
  Description_selected = c(
    "response to oxidative stress",
    "2-oxoglutarate metabolic process",
    "cellular response to oxidative stress",

    "pattern specification process",
    "regionalization",
    "anterior/posterior pattern specification",

    "myeloid leukocyte activation",
    "regulation of presynaptic membrane potential",
    "regulation of myeloid leukocyte differentiation",

    "DNA replication",
    "DNA-templated DNA replication",
    "DNA strand elongation involved in DNA replication",

    "response to endoplasmic reticulum stress",
    "endoplasmic reticulum unfolded protein response",
    "intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress",

    "protein folding in endoplasmic reticulum",
    "nucleocytoplasmic transport",
    "nuclear transport",
    "organelle fission"
  ),
  stringsAsFactors = FALSE
)

selected_kegg_terms <- data.frame(
  Cluster_label = c(
    rep("Cluster_1", 3),
    rep("Cluster_2", 3),
    rep("Cluster_3", 3),
    rep("Cluster_4", 3),
    rep("Cluster_5", 3),
    rep("Cluster_6", 4)
  ),
  Description_selected = c(
    "2-Oxocarboxylic acid metabolism",
    "Cholesterol metabolism",
    "Retinol metabolism",

    "Proteoglycans in cancer",
    "Cytoskeleton in muscle cells",
    "Focal adhesion",

    "Neuroactive ligand signaling",
    "Synaptic vesicle cycle",
    "Neuroactive ligand-receptor interaction",

    "DNA replication",
    "Mismatch repair",
    "Cell cycle",

    "IL-17 signaling pathway",
    "Apoptosis",
    "Apoptosis - multiple species",

    "Human immunodeficiency virus 1 infection",
    "FoxO signaling pathway",
    "Human T-cell leukemia virus 1 infection",
    "Protein processing in endoplasmic reticulum"
  ),
  stringsAsFactors = FALSE
)

clean_enrich_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  df <- df %>%
    mutate(
      Cluster = dplyr::case_when(
        "Cluster" %in% colnames(.) ~ as.character(.data[["Cluster"]]),
        "cluster_id" %in% colnames(.) ~ as.character(.data[["cluster_id"]]),
        "cluster" %in% colnames(.) ~ as.character(.data[["cluster"]]),
        TRUE ~ NA_character_
      ),
      Cluster_label = dplyr::case_when(
        "Cluster_label" %in% colnames(.) ~ pretty_cluster(.data[["Cluster_label"]]),
        TRUE ~ pretty_cluster(Cluster)
      ),
      Cluster_label = factor(Cluster_label, levels = cluster_levels),
      Description = standardize_term(Description),
      desc_std = str_to_lower(str_squish(Description)),
      Count = clean_num(Count),
      pvalue = clean_num(pvalue),
      p.adjust = if ("p.adjust" %in% colnames(.)) clean_num(`p.adjust`) else NA_real_,
      GeneRatio_num = if ("GeneRatio" %in% colnames(.)) parse_gene_ratio(GeneRatio) else NA_real_,
      log10p = -log10(pvalue)
    ) %>%
    filter(
      !is.na(Cluster_label),
      !is.na(Description),
      Description != "",
      !is.na(Count),
      !is.na(pvalue),
      pvalue < p_cutoff
    )

  df
}

select_manual_terms <- function(df, selected_terms, out_dir, missing_file) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  selected_terms2 <- selected_terms %>%
    mutate(
      Cluster_label = factor(Cluster_label, levels = cluster_levels),
      selected_std = str_to_lower(str_squish(standardize_term(Description_selected))),
      Selected_order = row_number()
    )

  df2 <- df %>%
    mutate(desc_std = str_to_lower(str_squish(standardize_term(Description))))

  plot_df <- selected_terms2 %>%
    left_join(
      df2,
      by = c("Cluster_label", "selected_std" = "desc_std")
    )

  missing_terms <- plot_df %>%
    filter(is.na(Description)) %>%
    select(Cluster_label, Description_selected)

  if (nrow(missing_terms) > 0) {
    write.csv(
      missing_terms,
      file.path(out_dir, missing_file),
      row.names = FALSE,
      quote = FALSE
    )
    warning("Some selected terms were not found. Check: ", missing_file)
  }

  plot_df %>%
    filter(!is.na(Description)) %>%
    mutate(
      Display = pretty_term(Description),
      Display_wrapped = wrap_term(Display, width = 38),
      Cluster_label = factor(as.character(Cluster_label), levels = cluster_levels),
      Cluster_display = factor(
        cluster_clean_label(Cluster_label),
        levels = cluster_labels_clean
      )
    )
}

theme_selected_bubble <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.45),
      panel.grid.major.y = element_line(color = "grey95", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.55),

      axis.title.x = element_text(size = base_size + 1, color = "black"),
      axis.title.y = element_blank(),

      axis.text.x = element_text(
        size = base_size,
        color = "black",
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.text.y = element_text(
        size = base_size - 1,
        color = "black",
        lineheight = 0.95
      ),

      plot.title = element_blank(),

      legend.title = element_text(size = base_size, color = "black"),
      legend.text = element_text(size = base_size - 1, color = "black"),
      legend.key = element_blank(),

      plot.margin = margin(10, 18, 10, 14)
    )
}

make_cluster_bubble_plot <- function(plot_df, count_breaks, fill_limits) {
  if (is.null(plot_df) || nrow(plot_df) == 0) return(NULL)

  plot_df <- plot_df %>%
    group_by(Cluster_label) %>%
    arrange(Count, pvalue, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(
      Display_wrapped = factor(Display_wrapped, levels = rev(unique(Display_wrapped)))
    )

  my_cols <- c("#253494", "#7A1FA2", "#C51B8A", "#E34A33", "#F03B20")

  ggplot(plot_df, aes(x = Cluster_display, y = Display_wrapped)) +
    geom_point(
      aes(size = Count, fill = log10p),
      shape = 21,
      color = "grey22",
      stroke = 0.32,
      alpha = 0.96
    ) +
    scale_size_continuous(
      range = c(3.8, 10.5),
      breaks = count_breaks,
      limits = range(plot_all_df$Count, na.rm = TRUE),
      name = "Gene count",
      guide = guide_legend(
        override.aes = list(
          shape = 16,
          colour = "black",
          fill = "black",
          alpha = 1,
          stroke = 0
        )
      )
    ) +
    scale_fill_gradientn(
      colours = my_cols,
      limits = fill_limits,
      name = expression(-log[10]("P value"))
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    labs(
      x = "maSigPro cluster",
      y = NULL
    ) +
    theme_selected_bubble(base_size = 13)
}

go_df_for_plot <- clean_enrich_df(go_summary)
if ("ONTOLOGY" %in% colnames(go_df_for_plot)) {
  go_df_for_plot <- go_df_for_plot %>% filter(ONTOLOGY == "BP")
}
if ("Ontology" %in% colnames(go_df_for_plot)) {
  go_df_for_plot <- go_df_for_plot %>% filter(Ontology == "BP")
}

kegg_df_for_plot <- clean_enrich_df(kegg_summary)

go_plot_df <- select_manual_terms(
  go_df_for_plot,
  selected_go_terms,
  selected_plot_dir,
  "missing_selected_GO_terms.csv"
)

kegg_plot_df <- select_manual_terms(
  kegg_df_for_plot,
  selected_kegg_terms,
  selected_plot_dir,
  "missing_selected_KEGG_terms.csv"
)

write.csv(
  go_plot_df,
  file.path(selected_plot_dir, "GO_BP_selected_terms_for_plot.csv"),
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  kegg_plot_df,
  file.path(selected_plot_dir, "KEGG_selected_terms_for_plot.csv"),
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)

plot_all_df <- bind_rows(
  go_plot_df %>%
    transmute(Source = "GO", Count = as.numeric(Count), log10p = as.numeric(log10p)),
  kegg_plot_df %>%
    transmute(Source = "KEGG", Count = as.numeric(Count), log10p = as.numeric(log10p))
) %>%
  filter(!is.na(Count), !is.na(log10p))

if (nrow(plot_all_df) > 0) {
  count_breaks <- pretty(plot_all_df$Count, n = 4)
  count_breaks <- count_breaks[
    count_breaks >= min(plot_all_df$Count, na.rm = TRUE) &
      count_breaks <= max(plot_all_df$Count, na.rm = TRUE)
  ]

  fill_limits <- range(plot_all_df$log10p, na.rm = TRUE)

  p_go_selected <- make_cluster_bubble_plot(go_plot_df, count_breaks, fill_limits)
  p_kegg_selected <- make_cluster_bubble_plot(kegg_plot_df, count_breaks, fill_limits)

  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

  if (!is.null(p_go_selected)) {
    ggsave(
      filename = file.path(selected_plot_dir, "GO_BP_selected_cluster_summary.pdf"),
      plot = p_go_selected,
      width = fig_width,
      height = fig_height,
      units = "in",
      device = pdf_device
    )

    ggsave(
      filename = file.path(selected_plot_dir, "GO_BP_selected_cluster_summary.png"),
      plot = p_go_selected,
      width = fig_width,
      height = fig_height,
      units = "in",
      dpi = 600
    )
  }

  if (!is.null(p_kegg_selected)) {
    ggsave(
      filename = file.path(selected_plot_dir, "KEGG_selected_cluster_summary.pdf"),
      plot = p_kegg_selected,
      width = fig_width,
      height = fig_height,
      units = "in",
      device = pdf_device
    )

    ggsave(
      filename = file.path(selected_plot_dir, "KEGG_selected_cluster_summary.png"),
      plot = p_kegg_selected,
      width = fig_width,
      height = fig_height,
      units = "in",
      dpi = 600
    )
  }

  ############################################################
  # 17) Combined enrichment figure
  ############################################################
  if (!is.null(p_go_selected) && !is.null(p_kegg_selected)) {
    p_enrichment_combined <- (p_go_selected / p_kegg_selected) +
      plot_layout(heights = c(1, 1), guides = "collect") &
      theme(legend.position = "right")

    ggsave(
      filename = file.path(work_dir, "06_combined_plots", "GO_KEGG_selected_cluster_summary_combined.pdf"),
      plot = p_enrichment_combined,
      width = 10.2,
      height = 14.6,
      units = "in",
      device = pdf_device
    )

    ggsave(
      filename = file.path(work_dir, "06_combined_plots", "GO_KEGG_selected_cluster_summary_combined.png"),
      plot = p_enrichment_combined,
      width = 10.2,
      height = 14.6,
      units = "in",
      dpi = 600
    )
  }
} else {
  warning("No selected GO/KEGG terms available for plotting. Selected summary plots were skipped.")
}

############################################################
# 18) Finish
############################################################
cat("\n============================================================\n")
cat("NiV maSigPro temporal analysis and cluster enrichment FINISHED\n")
cat("Final temporal cluster number: k = 6\n")
cat("Selected export sig table:", selected_export_path, "\n")
cat("Selected raw see.genes object:", selected_raw_path, "\n")
cat("Interpretation direction: Infection vs Control\n")
cat("Output root:", work_dir, "\n")
cat("Main output directories:\n")
cat("  01_maSigPro_tables/\n")
cat("  02_cluster_results/\n")
cat("  03_cluster_plots/\n")
cat("  04_enrichment_results/\n")
cat("  05_enrichment_plots/\n")
cat("  06_combined_plots/\n")
cat("============================================================\n")
