############################################################
# GSE46263_Metformin_GSE166707_516_signature_validation.R
#
# Purpose:
#   Use external metformin perturbation data from GSE46263
#   to validate GSE166707-derived 516 temporal core signatures.
#
# Primary comparison:
#   5.5 mM Glu + Met_Normal
#   vs
#   5.5 mM Glu_No Met_Normal
#
# Three validation layers:
#
# 1. Overall 516 expression-direction reversal
#    Core516_UP   expected NES < 0
#    Core516_DOWN expected NES > 0
#
# 2. Residual non-viral expression-direction reversal
#    Residual_UP   expected NES < 0
#    Residual_DOWN expected NES > 0
#
# 3. Virus-host factor functional direction validation
#    HRF_desired_up     expected NES > 0
#    HDF_desired_down   expected NES < 0
#
# Input:
#   D:/WeiRuan/NiV/二甲双胍GSE46263/probe.xlsx
#   D:/WeiRuan/NiV/二甲双胍GSE46263/count.xlsx
#   D:/WeiRuan/NiV/二甲双胍GSE46263/516带表达方向带病毒标签.csv
#
# Output:
#   D:/WeiRuan/NiV/二甲双胍GSE46263/metformin_GSE166707_516_validation_FINAL
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 600)
set.seed(123)

############################################################
# 0. Paths and parameters
############################################################

work_dir <- "D:/WeiRuan/NiV/二甲双胍GSE46263"
setwd(work_dir)

probe_file <- "D:/WeiRuan/NiV/二甲双胍GSE46263/probe.xlsx"
expr_file  <- "D:/WeiRuan/NiV/二甲双胍GSE46263/count.xlsx"
signature_file <- "D:/WeiRuan/NiV/二甲双胍GSE46263/516带表达方向带病毒标签.csv"

out_root <- file.path(work_dir, "metformin_GSE166707_516_validation_FINAL")

ABSENT_FLAG <- -9999

min_gs_size <- 5
max_gs_size <- 5000

top_n_heatmap <- 80

dir_list <- c(
  out_root,
  file.path(out_root, "00_preprocess_QC"),
  file.path(out_root, "01_limma_metformin"),
  file.path(out_root, "02_ranked_gene_list"),
  file.path(out_root, "03_signature_sets"),
  file.path(out_root, "04_GSEA_results"),
  file.path(out_root, "05_score_tables"),
  file.path(out_root, "06_leading_edge"),
  file.path(out_root, "07_plots"),
  file.path(out_root, "08_summary")
)

invisible(lapply(dir_list, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!file.exists(probe_file)) stop("Cannot find probe file: ", probe_file)
if (!file.exists(expr_file)) stop("Cannot find expression file: ", expr_file)
if (!file.exists(signature_file)) stop("Cannot find signature file: ", signature_file)

############################################################
# 1. Packages
############################################################

required_packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "limma",
  "clusterProfiler"
)

missing_packages <- required_packages[
  !sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  bioc_avail <- BiocManager::available()
  bioc_to_install <- intersect(missing_packages, bioc_avail)
  cran_to_install <- setdiff(missing_packages, bioc_to_install)
  
  if (length(bioc_to_install) > 0) {
    BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
  }
  
  if (length(cran_to_install) > 0) {
    install.packages(cran_to_install)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

############################################################
# 2. Helper functions
############################################################

write_csv_safe <- function(x, file) {
  write.csv(x, file, row.names = FALSE, quote = FALSE, na = "")
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  x <- stringr::str_split(
    x,
    "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*",
    simplify = TRUE
  )[, 1]
  
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x[x == "NA"] <- NA_character_
  x[x == "NULL"] <- NA_character_
  x
}

theme_sci <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black")
    )
}

get_quantiles <- function(mat, label) {
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  
  data.frame(
    step = label,
    min = as.numeric(quantile(v, 0, na.rm = TRUE)),
    q01 = as.numeric(quantile(v, 0.01, na.rm = TRUE)),
    q25 = as.numeric(quantile(v, 0.25, na.rm = TRUE)),
    median = as.numeric(quantile(v, 0.50, na.rm = TRUE)),
    q75 = as.numeric(quantile(v, 0.75, na.rm = TRUE)),
    q99 = as.numeric(quantile(v, 0.99, na.rm = TRUE)),
    max = as.numeric(quantile(v, 1, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

impute_row_median <- function(m) {
  m2 <- m
  for (i in seq_len(nrow(m2))) {
    x <- m2[i, ]
    if (anyNA(x)) {
      med <- median(x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      x[is.na(x)] <- med
      m2[i, ] <- x
    }
  }
  m2
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

find_gene_column <- function(df) {
  cn <- colnames(df)
  
  candidates <- c(
    "gene_symbol",
    "GeneSymbol",
    "Gene Symbol",
    "Gene symbol",
    "GENE_SYMBOL",
    "SYMBOL",
    "Symbol",
    "symbol",
    "Gene",
    "gene",
    "Genes",
    "genes",
    "gene_assignment",
    "Gene Assignment"
  )
  
  hit <- candidates[candidates %in% cn]
  if (length(hit) > 0) return(hit[1])
  
  hit2 <- cn[stringr::str_detect(tolower(cn), "symbol|gene")]
  if (length(hit2) > 0) return(hit2[1])
  
  stop("Cannot detect gene symbol column.")
}

break_ties_deterministic <- function(x) {
  x + seq_along(x) * 1e-12
}

format_p_label <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

get_expected_direction <- function(analysis_level, gene_set) {
  
  if (analysis_level %in% c("Overall_516_reversal", "Residual_nonviral_reversal")) {
    if (stringr::str_detect(gene_set, "_UP$")) return("NES < 0")
    if (stringr::str_detect(gene_set, "_DOWN$")) return("NES > 0")
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    if (gene_set == "HRF_desired_up") return("NES > 0")
    if (gene_set == "HDF_desired_down") return("NES < 0")
  }
  
  NA_character_
}

is_direction_matched <- function(analysis_level, gene_set, NES) {
  
  if (is.na(NES)) return(FALSE)
  
  if (analysis_level %in% c("Overall_516_reversal", "Residual_nonviral_reversal")) {
    if (stringr::str_detect(gene_set, "_UP$")) return(NES < 0)
    if (stringr::str_detect(gene_set, "_DOWN$")) return(NES > 0)
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    if (gene_set == "HRF_desired_up") return(NES > 0)
    if (gene_set == "HDF_desired_down") return(NES < 0)
  }
  
  FALSE
}

run_one_gsea <- function(geneList, term2gene, analysis_level, comparison) {
  
  term2gene_use <- term2gene %>%
    dplyr::filter(gene %in% names(geneList))
  
  if (nrow(term2gene_use) == 0) {
    return(NULL)
  }
  
  gsea_res <- tryCatch(
    clusterProfiler::GSEA(
      geneList = geneList,
      TERM2GENE = term2gene_use,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) {
      message("GSEA failed: ", analysis_level, " | ", comparison, " | ", e$message)
      NULL
    }
  )
  
  if (is.null(gsea_res)) return(NULL)
  
  gsea_df <- as.data.frame(gsea_res)
  
  if (nrow(gsea_df) == 0) return(NULL)
  
  gsea_df$analysis_level <- analysis_level
  gsea_df$comparison <- comparison
  
  gsea_df
}

score_one_level <- function(df_one) {
  
  analysis_level <- unique(df_one$analysis_level)
  comparison <- unique(df_one$comparison)
  
  if (length(analysis_level) != 1 || length(comparison) != 1) {
    stop("score_one_level requires one analysis_level and one comparison.")
  }
  
  get_nes <- function(term) {
    x <- df_one$NES[df_one$ID == term]
    if (length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }
  
  get_padj <- function(term) {
    x <- df_one$p.adjust[df_one$ID == term]
    if (length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }
  
  if (analysis_level == "Overall_516_reversal") {
    
    nes_up <- get_nes("Core516_UP")
    nes_down <- get_nes("Core516_DOWN")
    
    score <- (nes_down - nes_up) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "Core516_DOWN",
      term_negative_expected = "Core516_UP",
      NES_positive_expected = nes_down,
      NES_negative_expected = nes_up,
      FDR_positive_expected = get_padj("Core516_DOWN"),
      FDR_negative_expected = get_padj("Core516_UP"),
      score = score,
      both_direction_matched = (!is.na(nes_down) && nes_down > 0) &&
        (!is.na(nes_up) && nes_up < 0),
      score_name = "overall_reversal_score",
      interpretation = ifelse(
        (!is.na(nes_down) && nes_down > 0) && (!is.na(nes_up) && nes_up < 0),
        "Directional reversal supported",
        "Directional reversal not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  if (analysis_level == "Residual_nonviral_reversal") {
    
    nes_up <- get_nes("Residual_UP")
    nes_down <- get_nes("Residual_DOWN")
    
    score <- (nes_down - nes_up) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "Residual_DOWN",
      term_negative_expected = "Residual_UP",
      NES_positive_expected = nes_down,
      NES_negative_expected = nes_up,
      FDR_positive_expected = get_padj("Residual_DOWN"),
      FDR_negative_expected = get_padj("Residual_UP"),
      score = score,
      both_direction_matched = (!is.na(nes_down) && nes_down > 0) &&
        (!is.na(nes_up) && nes_up < 0),
      score_name = "residual_reversal_score",
      interpretation = ifelse(
        (!is.na(nes_down) && nes_down > 0) && (!is.na(nes_up) && nes_up < 0),
        "Residual signature reversal supported",
        "Residual signature reversal not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    
    nes_hrf <- get_nes("HRF_desired_up")
    nes_hdf <- get_nes("HDF_desired_down")
    
    score <- (nes_hrf - nes_hdf) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "HRF_desired_up",
      term_negative_expected = "HDF_desired_down",
      NES_positive_expected = nes_hrf,
      NES_negative_expected = nes_hdf,
      FDR_positive_expected = get_padj("HRF_desired_up"),
      FDR_negative_expected = get_padj("HDF_desired_down"),
      score = score,
      both_direction_matched = (!is.na(nes_hrf) && nes_hrf > 0) &&
        (!is.na(nes_hdf) && nes_hdf < 0),
      score_name = "virus_factor_modulation_score",
      interpretation = ifelse(
        (!is.na(nes_hrf) && nes_hrf > 0) && (!is.na(nes_hdf) && nes_hdf < 0),
        "Virus-host factor functional modulation supported",
        "Virus-host factor functional modulation not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  stop("Unknown analysis_level: ", analysis_level)
}

############################################################
# 3. GSE46263 sample information
############################################################

cat("============================================================\n")
cat("Building GSE46263 sample information\n")
cat("============================================================\n")

primary_no_met <- c("GSM1127793", "GSM1127817", "GSM1127839")
primary_met    <- c("GSM1127806", "GSM1127828", "GSM1127850")

sample_info <- data.frame(
  sample_id = c(primary_no_met, primary_met),
  group = c(rep("NoMet_NormalGlu", 3), rep("Met_NormalGlu", 3)),
  comparison_group = c(rep("NoMet", 3), rep("Met", 3)),
  glucose = "5.5mM_normal",
  oxygen = "normoxia",
  condition_note = c(
    "5.5mmM Glu_No Met_Normal_rep1",
    "5.5mmM Glu_No Met_Normal_rep2",
    "5.5mmM Glu_No Met_Normal_rep3",
    "5.5mmM Glu + Met_Normal_rep1",
    "5.5mmM Glu + Met_Normal_rep2",
    "5.5mmM Glu + Met_Normal_rep3"
  ),
  stringsAsFactors = FALSE
)

sample_info$comparison_group <- factor(
  sample_info$comparison_group,
  levels = c("NoMet", "Met")
)

rownames(sample_info) <- sample_info$sample_id

write_csv_safe(
  sample_info,
  file.path(out_root, "00_preprocess_QC", "sample_info_primary_normal_glucose.csv")
)

print(sample_info)

############################################################
# 4. Read expression matrix
############################################################

cat("============================================================\n")
cat("Reading expression matrix\n")
cat("File: ", expr_file, "\n", sep = "")
cat("============================================================\n")

expr_df <- readxl::read_excel(expr_file)

colnames(expr_df)[1] <- "ID"

expected_samples <- sample_info$sample_id

missing_samples <- setdiff(expected_samples, colnames(expr_df))

if (length(missing_samples) > 0) {
  stop(
    "count.xlsx missing expected sample columns: ",
    paste(missing_samples, collapse = ", ")
  )
}

expr_df <- expr_df %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::filter(!is.na(ID), ID != "") %>%
  dplyr::filter(!stringr::str_starts(ID, "!"))

for (cc in expected_samples) {
  expr_df[[cc]] <- as.character(expr_df[[cc]])
  expr_df[[cc]][expr_df[[cc]] %in% c("NULL", "null", "NA", "NaN", "")] <- NA
  expr_df[[cc]] <- suppressWarnings(as.numeric(expr_df[[cc]]))
}

expr_mat_raw <- as.matrix(expr_df[, expected_samples])
rownames(expr_mat_raw) <- expr_df$ID
storage.mode(expr_mat_raw) <- "numeric"

cat("Raw expression matrix for primary comparison: ",
    nrow(expr_mat_raw), " probes/features x ",
    ncol(expr_mat_raw), " samples\n", sep = "")

write_csv_safe(
  data.frame(
    item = c("expr_file", "probe_file", "signature_file", "n_features_raw", "n_samples_used"),
    value = c(expr_file, probe_file, signature_file, nrow(expr_mat_raw), ncol(expr_mat_raw))
  ),
  file.path(out_root, "00_preprocess_QC", "input_file_summary.csv")
)

write_csv_safe(
  get_quantiles(expr_mat_raw, "raw_primary_samples"),
  file.path(out_root, "00_preprocess_QC", "expression_quantiles_raw_primary_samples.csv")
)

############################################################
# 5. Read probe annotation
############################################################

cat("============================================================\n")
cat("Reading probe annotation\n")
cat("File: ", probe_file, "\n", sep = "")
cat("============================================================\n")

probe_df_raw <- readxl::read_excel(probe_file)

colnames(probe_df_raw)[1] <- "ID"

gene_symbol_col <- find_gene_column(probe_df_raw)

probe_df <- probe_df_raw %>%
  dplyr::transmute(
    ID = as.character(ID),
    gene_symbol = clean_symbol(.data[[gene_symbol_col]])
  ) %>%
  dplyr::filter(
    !is.na(ID),
    ID != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  )

cat("Probe annotation rows after cleaning: ", nrow(probe_df), "\n", sep = "")
cat("Gene symbol column used: ", gene_symbol_col, "\n", sep = "")

write_csv_safe(
  data.frame(
    gene_symbol_col_used = gene_symbol_col,
    probe_annotation_rows_cleaned = nrow(probe_df),
    unique_gene_symbols = length(unique(probe_df$gene_symbol))
  ),
  file.path(out_root, "00_preprocess_QC", "probe_annotation_summary.csv")
)

############################################################
# 6. Preprocess expression matrix
############################################################

cat("============================================================\n")
cat("Preprocessing expression matrix\n")
cat("============================================================\n")

expr_mat <- expr_mat_raw

absent_n <- sum(expr_mat == ABSENT_FLAG, na.rm = TRUE)
expr_mat[expr_mat == ABSENT_FLAG] <- NA
na_after_absent <- sum(is.na(expr_mat))

# GCRMA data are usually log2-like.
# Still check distribution and avoid unnecessary log2 transformation.
q_before <- get_quantiles(expr_mat, "after_absent_to_NA_before_log2_decision")

x_nonNA <- as.numeric(expr_mat[!is.na(expr_mat)])
q99_before <- as.numeric(quantile(x_nonNA, 0.99, names = FALSE, na.rm = TRUE))
xmax_before <- max(x_nonNA, na.rm = TRUE)
xmin_before <- min(x_nonNA, na.rm = TRUE)

is_log2_like <- (q99_before < 50 && xmax_before < 100)

expr_processed <- expr_mat
offset_used <- 0

if (!is_log2_like) {
  min_nonNA <- suppressWarnings(min(expr_processed, na.rm = TRUE))
  offset_used <- ifelse(
    is.finite(min_nonNA) && min_nonNA <= 0,
    abs(min_nonNA) + 0.001,
    0
  )
  
  expr_shifted <- expr_processed
  expr_shifted[!is.na(expr_shifted)] <- expr_shifted[!is.na(expr_shifted)] + offset_used
  expr_processed <- log2(expr_shifted)
  
  transform_note <- "Applied offset + log2 because expression was not log2-like."
} else {
  transform_note <- "Expression was log2-like; kept as-is. GEO processing reports GCRMA normalization."
}

q_after <- get_quantiles(expr_processed, "final_expression_after_log2_decision")

# Keep features with enough data in the 6 primary samples
keep_feature <- rowSums(!is.na(expr_processed)) >= 4
expr_processed1 <- expr_processed[keep_feature, , drop = FALSE]

# Match annotation
common_ids <- intersect(rownames(expr_processed1), probe_df$ID)

expr_anno <- expr_processed1[common_ids, , drop = FALSE]
probe_sub <- probe_df[match(common_ids, probe_df$ID), , drop = FALSE]

idx_anno <- !is.na(probe_sub$gene_symbol) &
  probe_sub$gene_symbol != "" &
  probe_sub$gene_symbol != "NA" &
  probe_sub$gene_symbol != "NULL"

expr_anno <- expr_anno[idx_anno, , drop = FALSE]
probe_sub <- probe_sub[idx_anno, , drop = FALSE]

qc_filter_stat <- data.frame(
  total_rows_in_counts = nrow(expr_mat_raw),
  absent_flag_count = absent_n,
  NA_after_absent_to_NA = na_after_absent,
  rows_after_at_least_4_of_6_present = nrow(expr_processed1),
  rows_after_id_match = length(common_ids),
  rows_after_annotation_filter = nrow(expr_anno),
  unique_gene_symbols_after_annotation = length(unique(probe_sub$gene_symbol)),
  q99_before_log2 = q99_before,
  max_before_log2 = xmax_before,
  min_before_log2 = xmin_before,
  is_log2_like = is_log2_like,
  log2_applied = !is_log2_like,
  offset_used = offset_used,
  additional_quantile_normalization = FALSE,
  normalizeBetweenArrays = FALSE,
  stringsAsFactors = FALSE
)

write_csv_safe(
  rbind(q_before, q_after),
  file.path(out_root, "00_preprocess_QC", "expression_quantiles_before_after_log2_decision.csv")
)

write_csv_safe(
  qc_filter_stat,
  file.path(out_root, "00_preprocess_QC", "QC_filter_transform_summary.csv")
)

writeLines(
  c(
    "GSE46263 preprocessing notes",
    "============================================================",
    "Primary comparison: 5.5 mM Glu + Met_Normal vs 5.5 mM Glu_No Met_Normal",
    "GEO data processing reports: Partek GS 6.5 and GCRMA normalization.",
    paste0("Absent flag: ", ABSENT_FLAG, " -> NA"),
    paste0("Absent flag count: ", absent_n),
    paste0("NA after absent-to-NA: ", na_after_absent),
    paste0("q99 before log2 decision: ", q99_before),
    paste0("max before log2 decision: ", xmax_before),
    paste0("min before log2 decision: ", xmin_before),
    paste0("is_log2_like: ", is_log2_like),
    paste0("log2 applied: ", !is_log2_like),
    paste0("offset used: ", offset_used),
    transform_note,
    "Additional quantile normalization: FALSE",
    "normalizeBetweenArrays: FALSE"
  ),
  file.path(out_root, "00_preprocess_QC", "preprocess_notes.txt")
)

print(qc_filter_stat)

############################################################
# 7. Probe-to-gene collapse
############################################################

cat("============================================================\n")
cat("Probe-to-gene collapse\n")
cat("============================================================\n")

tmp <- data.frame(
  probe_id = rownames(expr_anno),
  gene_symbol = probe_sub$gene_symbol,
  expr_anno,
  check.names = FALSE
)

tmp$mean_expr <- rowMeans(
  tmp[, 3:(ncol(tmp) - 1)],
  na.rm = TRUE
)

gene_unique <- tmp %>%
  dplyr::filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    is.finite(mean_expr)
  ) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

expr_gene <- as.matrix(gene_unique[, 3:(ncol(gene_unique) - 1)])
rownames(expr_gene) <- gene_unique$gene_symbol
storage.mode(expr_gene) <- "numeric"

expr_gene <- expr_gene[, rownames(sample_info), drop = FALSE]

cat("Gene-level matrix: ",
    nrow(expr_gene), " genes x ",
    ncol(expr_gene), " samples\n", sep = "")

write_csv_safe(
  gene_unique %>%
    dplyr::select(probe_id, gene_symbol, mean_expr),
  file.path(out_root, "00_preprocess_QC", "probe_to_gene_collapse_mapping.csv")
)

write_csv_safe(
  data.frame(gene_symbol = rownames(expr_gene), expr_gene, check.names = FALSE),
  file.path(out_root, "00_preprocess_QC", "expr_gene_matrix_primary_normal_glucose.csv")
)

############################################################
# 8. QC plots
############################################################

cat("============================================================\n")
cat("QC plots\n")
cat("============================================================\n")

box_df <- data.frame(
  value = as.vector(expr_gene),
  sample = rep(colnames(expr_gene), each = nrow(expr_gene))
) %>%
  dplyr::left_join(sample_info, by = c("sample" = "sample_id")) %>%
  dplyr::filter(is.finite(value))

p_box <- ggplot(box_df, aes(x = sample, y = value, fill = comparison_group)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("NoMet" = "grey70", "Met" = "#6A9BCB")) +
  labs(
    title = "GSE46263 expression distribution",
    subtitle = "Primary comparison: normal glucose, metformin vs no metformin",
    x = NULL,
    y = "Expression",
    fill = NULL
  ) +
  theme_sci(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top"
  )

ggsave(
  file.path(out_root, "00_preprocess_QC", "Boxplot_gene_level_primary_normal_glucose_SCI.png"),
  p_box,
  width = 8.5,
  height = 5.3,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "00_preprocess_QC", "Boxplot_gene_level_primary_normal_glucose_SCI.pdf"),
  p_box,
  width = 8.5,
  height = 5.3
)

expr_pca <- impute_row_median(expr_gene)
expr_pca <- expr_pca[
  apply(expr_pca, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

pca <- prcomp(t(expr_pca), scale. = FALSE)
var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  comparison_group = sample_info[rownames(pca$x), "comparison_group"],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

write_csv_safe(
  pca_df,
  file.path(out_root, "00_preprocess_QC", "PCA_primary_normal_glucose.csv")
)

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = comparison_group, label = sample_id)
) +
  geom_point(size = 4.2, alpha = 0.95) +
  ggrepel::geom_text_repel(size = 3.3, max.overlaps = Inf) +
  scale_color_manual(values = c("NoMet" = "grey40", "Met" = "#6A9BCB")) +
  labs(
    title = "PCA of GSE46263 primary samples",
    subtitle = "Normal glucose condition",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = NULL
  ) +
  theme_sci(12) +
  theme(legend.position = "top")

ggsave(
  file.path(out_root, "00_preprocess_QC", "PCA_primary_normal_glucose_SCI.png"),
  p_pca,
  width = 7,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "00_preprocess_QC", "PCA_primary_normal_glucose_SCI.pdf"),
  p_pca,
  width = 7,
  height = 6
)

############################################################
# 9. limma: Metformin vs No Metformin
############################################################

cat("============================================================\n")
cat("limma: Metformin vs No Metformin under normal glucose\n")
cat("============================================================\n")

design <- model.matrix(~ 0 + comparison_group, data = sample_info)
colnames(design) <- gsub("^comparison_group", "", colnames(design))

contrast_matrix <- limma::makeContrasts(
  Met_vs_NoMet_NormalGlu = Met - NoMet,
  levels = design
)

write_csv_safe(
  data.frame(sample_id = rownames(design), design, check.names = FALSE),
  file.path(out_root, "01_limma_metformin", "design_matrix.csv")
)

write_csv_safe(
  data.frame(contrast = rownames(contrast_matrix), contrast_matrix, check.names = FALSE),
  file.path(out_root, "01_limma_metformin", "contrast_matrix.csv")
)

fit <- limma::lmFit(expr_gene, design)
fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

coef_name <- "Met_vs_NoMet_NormalGlu"

deg_all <- limma::topTable(
  fit2,
  coef = coef_name,
  adjust.method = "fdr",
  number = Inf,
  sort.by = "P"
)

deg_all$gene_symbol <- rownames(deg_all)

deg_all <- deg_all %>%
  dplyr::select(
    gene_symbol,
    logFC,
    AveExpr,
    t,
    P.Value,
    adj.P.Val,
    B
  ) %>%
  dplyr::mutate(
    comparison = coef_name,
    direction = dplyr::case_when(
      logFC > 0 ~ "Up_in_metformin",
      logFC < 0 ~ "Down_in_metformin",
      TRUE ~ "Neutral"
    )
  )

write_csv_safe(
  deg_all,
  file.path(out_root, "01_limma_metformin", "DEG_all_Met_vs_NoMet_NormalGlu.csv")
)

rank_df <- deg_all %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "", !is.na(t), is.finite(t)) %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(t)) %>%
  dplyr::select(gene_symbol, t, logFC, P.Value, adj.P.Val, comparison)

rank_df$t_ranked <- break_ties_deterministic(rank_df$t)

write_csv_safe(
  rank_df,
  file.path(out_root, "02_ranked_gene_list", "ranked_gene_list_t_Met_vs_NoMet_NormalGlu.csv")
)

geneList <- rank_df$t_ranked
names(geneList) <- rank_df$gene_symbol
geneList <- sort(geneList, decreasing = TRUE)

############################################################
# 10. Load 516 signature table
############################################################

cat("============================================================\n")
cat("Loading GSE166707 516 signature table\n")
cat("File: ", signature_file, "\n", sep = "")
cat("============================================================\n")

sig_df <- read.csv(signature_file, check.names = FALSE)

required_cols <- c("Genes", "direction_type", "Host factor type")
missing_cols <- setdiff(required_cols, colnames(sig_df))

if (length(missing_cols) > 0) {
  stop(
    "Signature file missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

sig_df <- sig_df %>%
  dplyr::mutate(
    Genes = clean_symbol(Genes),
    direction_type = as.character(direction_type),
    direction_type = stringr::str_trim(direction_type),
    host_factor_type_raw = as.character(`Host factor type`),
    host_factor_type_raw = stringr::str_trim(host_factor_type_raw),
    host_factor_type = dplyr::case_when(
      toupper(host_factor_type_raw) == "HDF" ~ "HDF",
      toupper(host_factor_type_raw) == "HRF" ~ "HRF",
      TRUE ~ "Non_viral"
    )
  ) %>%
  dplyr::filter(!is.na(Genes), Genes != "") %>%
  dplyr::distinct(Genes, .keep_all = TRUE)

############################################################
# 11. Build three layers of gene sets
############################################################

# Layer 1: Overall 516 expression-direction reversal
core_up <- sig_df %>%
  dplyr::filter(direction_type == "all_up") %>%
  dplyr::pull(Genes) %>%
  unique()

core_down <- sig_df %>%
  dplyr::filter(direction_type == "all_down") %>%
  dplyr::pull(Genes) %>%
  unique()

# Layer 2: Residual non-viral expression-direction reversal
residual_up <- sig_df %>%
  dplyr::filter(host_factor_type == "Non_viral", direction_type == "all_up") %>%
  dplyr::pull(Genes) %>%
  unique()

residual_down <- sig_df %>%
  dplyr::filter(host_factor_type == "Non_viral", direction_type == "all_down") %>%
  dplyr::pull(Genes) %>%
  unique()

# Layer 3: Virus-host factor functional direction
# HRF: desired up
# HDF: desired down
hrf_desired_up <- sig_df %>%
  dplyr::filter(host_factor_type == "HRF") %>%
  dplyr::pull(Genes) %>%
  unique()

hdf_desired_down <- sig_df %>%
  dplyr::filter(host_factor_type == "HDF") %>%
  dplyr::pull(Genes) %>%
  unique()

signature_summary <- data.frame(
  analysis_level = c(
    "Overall_516_reversal",
    "Overall_516_reversal",
    "Residual_nonviral_reversal",
    "Residual_nonviral_reversal",
    "Virus_factor_functional_direction",
    "Virus_factor_functional_direction"
  ),
  gene_set = c(
    "Core516_UP",
    "Core516_DOWN",
    "Residual_UP",
    "Residual_DOWN",
    "HRF_desired_up",
    "HDF_desired_down"
  ),
  expected_direction = c(
    "NES < 0",
    "NES > 0",
    "NES < 0",
    "NES > 0",
    "NES > 0",
    "NES < 0"
  ),
  biological_meaning = c(
    "Metformin suppresses GSE166707 core genes upregulated during NiV infection",
    "Metformin restores GSE166707 core genes downregulated during NiV infection",
    "Metformin suppresses non-viral-label residual upregulated host-response genes",
    "Metformin restores non-viral-label residual downregulated host-response genes",
    "Metformin upregulates host restriction / antiviral factors",
    "Metformin downregulates host dependency / proviral factors"
  ),
  n_genes = c(
    length(core_up),
    length(core_down),
    length(residual_up),
    length(residual_down),
    length(hrf_desired_up),
    length(hdf_desired_down)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  signature_summary,
  file.path(out_root, "03_signature_sets", "GSE166707_516_signature_set_summary.csv")
)

write_csv_safe(data.frame(gene_symbol = core_up),
               file.path(out_root, "03_signature_sets", "Core516_UP.csv"))
write_csv_safe(data.frame(gene_symbol = core_down),
               file.path(out_root, "03_signature_sets", "Core516_DOWN.csv"))
write_csv_safe(data.frame(gene_symbol = residual_up),
               file.path(out_root, "03_signature_sets", "Residual_UP.csv"))
write_csv_safe(data.frame(gene_symbol = residual_down),
               file.path(out_root, "03_signature_sets", "Residual_DOWN.csv"))
write_csv_safe(data.frame(gene_symbol = hrf_desired_up),
               file.path(out_root, "03_signature_sets", "HRF_desired_up.csv"))
write_csv_safe(data.frame(gene_symbol = hdf_desired_down),
               file.path(out_root, "03_signature_sets", "HDF_desired_down.csv"))

print(signature_summary)

term2gene_overall <- data.frame(
  term = c(
    rep("Core516_UP", length(core_up)),
    rep("Core516_DOWN", length(core_down))
  ),
  gene = c(core_up, core_down),
  analysis_level = "Overall_516_reversal",
  stringsAsFactors = FALSE
)

term2gene_residual <- data.frame(
  term = c(
    rep("Residual_UP", length(residual_up)),
    rep("Residual_DOWN", length(residual_down))
  ),
  gene = c(residual_up, residual_down),
  analysis_level = "Residual_nonviral_reversal",
  stringsAsFactors = FALSE
)

term2gene_virus <- data.frame(
  term = c(
    rep("HRF_desired_up", length(hrf_desired_up)),
    rep("HDF_desired_down", length(hdf_desired_down))
  ),
  gene = c(hrf_desired_up, hdf_desired_down),
  analysis_level = "Virus_factor_functional_direction",
  stringsAsFactors = FALSE
)

term2gene_all <- dplyr::bind_rows(
  term2gene_overall,
  term2gene_residual,
  term2gene_virus
)

write_csv_safe(
  term2gene_all,
  file.path(out_root, "03_signature_sets", "TERM2GENE_all_signature_layers.csv")
)

mapped_summary <- term2gene_all %>%
  dplyr::group_by(analysis_level, term) %>%
  dplyr::summarise(
    n_original = dplyr::n_distinct(gene),
    n_mapped_to_GSE46263_ranked_list = length(intersect(unique(gene), names(geneList))),
    .groups = "drop"
  )

write_csv_safe(
  mapped_summary,
  file.path(out_root, "03_signature_sets", "signature_mapping_to_GSE46263_ranked_list_summary.csv")
)

print(mapped_summary)

############################################################
# 12. Run GSEA for three layers
############################################################

cat("============================================================\n")
cat("Running GSEA validation\n")
cat("============================================================\n")

gsea_overall <- run_one_gsea(
  geneList = geneList,
  term2gene = term2gene_overall,
  analysis_level = "Overall_516_reversal",
  comparison = coef_name
)

gsea_residual <- run_one_gsea(
  geneList = geneList,
  term2gene = term2gene_residual,
  analysis_level = "Residual_nonviral_reversal",
  comparison = coef_name
)

gsea_virus <- run_one_gsea(
  geneList = geneList,
  term2gene = term2gene_virus,
  analysis_level = "Virus_factor_functional_direction",
  comparison = coef_name
)

gsea_all <- dplyr::bind_rows(
  gsea_overall,
  gsea_residual,
  gsea_virus
)

if (is.null(gsea_all) || nrow(gsea_all) == 0) {
  stop("No GSEA results generated. Please check gene sets and ranked list.")
}

gsea_all <- gsea_all %>%
  dplyr::mutate(
    expected_direction = mapply(get_expected_direction, analysis_level, ID),
    direction_matched = mapply(is_direction_matched, analysis_level, ID, NES),
    direction_status = ifelse(direction_matched, "Matched", "Not matched")
  )

write_csv_safe(
  gsea_all,
  file.path(out_root, "04_GSEA_results", "GSEA_all_layers_Met_vs_NoMet_NormalGlu_long.csv")
)

score_table <- gsea_all %>%
  dplyr::group_by(analysis_level, comparison) %>%
  dplyr::group_split() %>%
  lapply(score_one_level) %>%
  dplyr::bind_rows()

score_table <- score_table %>%
  dplyr::mutate(
    analysis_level = factor(
      analysis_level,
      levels = c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )
    )
  ) %>%
  dplyr::arrange(analysis_level)

write_csv_safe(
  score_table,
  file.path(out_root, "05_score_tables", "Metformin_GSE166707_516_validation_score_table.csv")
)

print(score_table)

############################################################
# 13. Leading-edge genes
############################################################

leading_edge_df <- gsea_all %>%
  dplyr::select(
    analysis_level,
    comparison,
    ID,
    NES,
    pvalue,
    p.adjust,
    expected_direction,
    direction_matched,
    core_enrichment
  ) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(gene_symbol = core_enrichment) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

write_csv_safe(
  leading_edge_df,
  file.path(out_root, "06_leading_edge", "GSE166707_516_metformin_GSEA_leading_edge_genes.csv")
)

leading_edge_count <- leading_edge_df %>%
  dplyr::count(analysis_level, ID, gene_symbol, sort = TRUE)

write_csv_safe(
  leading_edge_count,
  file.path(out_root, "06_leading_edge", "GSE166707_516_metformin_leading_edge_gene_frequency.csv")
)

############################################################
# 14. SCI-style figures
############################################################

cat("============================================================\n")
cat("Plotting SCI-style figures\n")
cat("============================================================\n")

analysis_label_map <- c(
  "Overall_516_reversal" = "Overall 516 reversal",
  "Residual_nonviral_reversal" = "Residual non-viral reversal",
  "Virus_factor_functional_direction" = "Virus-host factor direction"
)

gene_set_label_map <- c(
  "Core516_UP" = "Core516 UP",
  "Core516_DOWN" = "Core516 DOWN",
  "Residual_UP" = "Residual UP",
  "Residual_DOWN" = "Residual DOWN",
  "HRF_desired_up" = "HRF desired-up",
  "HDF_desired_down" = "HDF desired-down"
)

plot_df <- gsea_all %>%
  dplyr::mutate(
    analysis_label = factor(
      analysis_label_map[analysis_level],
      levels = analysis_label_map[c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )]
    ),
    gene_set_label = gene_set_label_map[ID],
    neg_log10_fdr = -log10(pmax(p.adjust, 1e-300)),
    direction_status = ifelse(direction_matched, "Matched", "Not matched"),
    label_text = paste0(
      "NES=", sprintf("%.2f", NES),
      "\nFDR=", format_p_label(p.adjust)
    )
  )

# Main NES dotplot
p_dot <- ggplot(
  plot_df,
  aes(x = NES, y = gene_set_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_point(
    aes(
      color = gene_set_label,
      shape = direction_status,
      size = neg_log10_fdr
    ),
    alpha = 0.92,
    stroke = 0.9
  ) +
  geom_text(
    aes(label = label_text),
    nudge_y = 0.23,
    size = 3.1,
    color = "black"
  ) +
  facet_grid(analysis_label ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c("Matched" = 16, "Not matched" = 1)) +
  scale_size_continuous(range = c(3.2, 7.2)) +
  labs(
    title = "Metformin perturbation against GSE166707-derived NiV signatures",
    subtitle = "GSE46263 normal-glucose HUVECs: metformin vs no metformin",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    color = NULL,
    shape = NULL,
    size = expression(-log[10](FDR))
  ) +
  theme_sci(11) +
  theme(
    legend.position = "right",
    strip.text.y = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10.5)
  )

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_GSEA_NES_dotplot_SCI.png"),
  p_dot,
  width = 10.5,
  height = 8.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_GSEA_NES_dotplot_SCI.pdf"),
  p_dot,
  width = 10.5,
  height = 8.2
)

# Score barplot
score_plot_df <- score_table %>%
  dplyr::mutate(
    analysis_label = factor(
      analysis_label_map[as.character(analysis_level)],
      levels = analysis_label_map[c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )]
    ),
    match_status = ifelse(
      both_direction_matched,
      "Both directions matched",
      "Partial or no match"
    )
  )

p_score <- ggplot(
  score_plot_df,
  aes(x = analysis_label, y = score, fill = match_status)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.2f", score)),
    vjust = ifelse(score_plot_df$score >= 0, -0.35, 1.25),
    size = 3.6
  ) +
  scale_fill_manual(
    values = c(
      "Both directions matched" = "#5AAE61",
      "Partial or no match" = "grey72"
    )
  ) +
  labs(
    title = "Metformin signature-level directional scores",
    subtitle = "Reversal scores for core/residual signatures; functional modulation score for virus-host factors",
    x = NULL,
    y = "Directional score",
    fill = NULL
  ) +
  theme_sci(11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10.5)
  )

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_directional_score_barplot_SCI.png"),
  p_score,
  width = 9.2,
  height = 5.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_directional_score_barplot_SCI.pdf"),
  p_score,
  width = 9.2,
  height = 5.8
)

# Heatmap-style score overview
p_heat <- ggplot(
  score_plot_df,
  aes(x = "Metformin\nnormal glucose", y = analysis_label, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 4, color = "black") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  labs(
    title = "Directional score overview",
    subtitle = "Positive scores indicate expected reversal/modulation direction",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_sci(11) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_score_heatmap_SCI.png"),
  p_heat,
  width = 6.8,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "07_plots", "Metformin_GSE166707_516_score_heatmap_SCI.pdf"),
  p_heat,
  width = 6.8,
  height = 4.5
)

############################################################
# 15. Optional leading-edge heatmap
############################################################

if (nrow(leading_edge_df) > 0) {
  
  leading_gene_count_df <- leading_edge_df %>%
    dplyr::count(gene_symbol, sort = TRUE)
  
  leading_gene_count <- min(top_n_heatmap, nrow(leading_gene_count_df))
  
  leading_genes_top <- leading_gene_count_df %>%
    dplyr::slice_head(n = leading_gene_count) %>%
    dplyr::pull(gene_symbol)
  
  leading_genes_top <- intersect(leading_genes_top, rownames(expr_gene))
  
  write_csv_safe(
    data.frame(gene_symbol = leading_genes_top),
    file.path(out_root, "06_leading_edge", "top_leading_edge_genes_for_heatmap.csv")
  )
  
  if (length(leading_genes_top) >= 2) {
    
    heat_mat <- expr_gene[leading_genes_top, , drop = FALSE]
    heat_mat <- impute_row_median(heat_mat)
    heat_z <- row_zscore(heat_mat)
    
    ann_col <- data.frame(Group = sample_info$comparison_group)
    rownames(ann_col) <- sample_info$sample_id
    
    annotation_colors <- list(
      Group = c(
        "NoMet" = "grey60",
        "Met" = "#6A9BCB"
      )
    )
    
    pheatmap::pheatmap(
      heat_z,
      annotation_col = ann_col,
      annotation_colors = annotation_colors,
      show_colnames = TRUE,
      show_rownames = length(leading_genes_top) <= 80,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      fontsize = 9,
      fontsize_row = 6,
      fontsize_col = 8,
      border_color = NA,
      main = "Leading-edge genes from metformin GSEA",
      filename = file.path(out_root, "07_plots", "Metformin_GSE166707_516_leading_edge_heatmap_SCI.png"),
      width = 8.5,
      height = 10
    )
    
    pdf(
      file.path(out_root, "07_plots", "Metformin_GSE166707_516_leading_edge_heatmap_SCI.pdf"),
      width = 8.5,
      height = 10
    )
    pheatmap::pheatmap(
      heat_z,
      annotation_col = ann_col,
      annotation_colors = annotation_colors,
      show_colnames = TRUE,
      show_rownames = length(leading_genes_top) <= 80,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      fontsize = 9,
      fontsize_row = 6,
      fontsize_col = 8,
      border_color = NA,
      main = "Leading-edge genes from metformin GSEA"
    )
    dev.off()
  }
}

############################################################
# 16. Summary report
############################################################

analysis_report <- c(
  "GSE46263 metformin validation of GSE166707 516-core signatures",
  "============================================================",
  paste0("Working directory: ", work_dir),
  paste0("Expression file: ", expr_file),
  paste0("Probe annotation file: ", probe_file),
  paste0("Signature file: ", signature_file),
  paste0("Output root: ", out_root),
  "",
  "Dataset:",
  "GSE46263 includes 66 primary HUVEC samples across 22 conditions, each with 3 biological replicates.",
  "GEO reports Partek GS 6.5 and GCRMA as the normalization method.",
  "",
  "Primary comparison used in this script:",
  "5.5 mM Glu + Met_Normal vs 5.5 mM Glu_No Met_Normal",
  "",
  "No Met normal-glucose samples:",
  paste(primary_no_met, collapse = ", "),
  "Met normal-glucose samples:",
  paste(primary_met, collapse = ", "),
  "",
  "limma model:",
  "Design: ~ 0 + comparison_group",
  "Contrast: Met - NoMet under normal glucose",
  "eBayes: trend = TRUE, robust = TRUE",
  "",
  "Three validation layers:",
  "",
  "1. Overall 516 expression-direction reversal",
  "   Core516_UP expected NES < 0",
  "   Core516_DOWN expected NES > 0",
  "   Score = (NES_Core516_DOWN - NES_Core516_UP) / 2",
  "",
  "2. Residual non-viral expression-direction reversal",
  "   Residual_UP expected NES < 0",
  "   Residual_DOWN expected NES > 0",
  "   Score = (NES_Residual_DOWN - NES_Residual_UP) / 2",
  "",
  "3. Virus-host factor functional direction validation",
  "   HRF_desired_up expected NES > 0",
  "   HDF_desired_down expected NES < 0",
  "   Score = (NES_HRF_desired_up - NES_HDF_desired_down) / 2",
  "",
  "Why HRF/HDF use functional directions:",
  "In the original CMap prediction, HRFs were placed into down genes and HDFs were placed into up genes.",
  "For negative connectivity drugs, this implies that the expected drug effect is to upregulate HRFs",
  "and downregulate HDFs. Therefore, validation of HRF/HDF follows this predefined functional direction",
  "rather than the original infection logFC direction.",
  "",
  "Main output tables:",
  "01_limma_metformin/DEG_all_Met_vs_NoMet_NormalGlu.csv",
  "02_ranked_gene_list/ranked_gene_list_t_Met_vs_NoMet_NormalGlu.csv",
  "03_signature_sets/GSE166707_516_signature_set_summary.csv",
  "04_GSEA_results/GSEA_all_layers_Met_vs_NoMet_NormalGlu_long.csv",
  "05_score_tables/Metformin_GSE166707_516_validation_score_table.csv",
  "06_leading_edge/GSE166707_516_metformin_GSEA_leading_edge_genes.csv",
  "",
  "Main figures:",
  "07_plots/Metformin_GSE166707_516_GSEA_NES_dotplot_SCI.pdf/png",
  "07_plots/Metformin_GSE166707_516_directional_score_barplot_SCI.pdf/png",
  "07_plots/Metformin_GSE166707_516_score_heatmap_SCI.pdf/png",
  "07_plots/Metformin_GSE166707_516_leading_edge_heatmap_SCI.pdf/png"
)

writeLines(
  analysis_report,
  file.path(out_root, "08_summary", "GSE46263_metformin_GSE166707_516_validation_report.txt")
)

write_csv_safe(
  qc_filter_stat,
  file.path(out_root, "08_summary", "QC_filter_transform_summary.csv")
)

write_csv_safe(
  signature_summary,
  file.path(out_root, "08_summary", "GSE166707_516_signature_set_summary.csv")
)

write_csv_safe(
  mapped_summary,
  file.path(out_root, "08_summary", "GSE166707_516_signature_mapping_summary.csv")
)

write_csv_safe(
  score_table,
  file.path(out_root, "08_summary", "Metformin_GSE166707_516_validation_score_table.csv")
)

cat("============================================================\n")
cat("GSE46263 metformin validation finished.\n")
cat("Output root:\n")
cat(out_root, "\n\n")
cat("Key outputs:\n")
cat("  00_preprocess_QC/sample_info_primary_normal_glucose.csv\n")
cat("  00_preprocess_QC/expr_gene_matrix_primary_normal_glucose.csv\n")
cat("  01_limma_metformin/DEG_all_Met_vs_NoMet_NormalGlu.csv\n")
cat("  02_ranked_gene_list/ranked_gene_list_t_Met_vs_NoMet_NormalGlu.csv\n")
cat("  03_signature_sets/GSE166707_516_signature_set_summary.csv\n")
cat("  04_GSEA_results/GSEA_all_layers_Met_vs_NoMet_NormalGlu_long.csv\n")
cat("  05_score_tables/Metformin_GSE166707_516_validation_score_table.csv\n")
cat("  07_plots/Metformin_GSE166707_516_GSEA_NES_dotplot_SCI.png/pdf\n")
cat("  07_plots/Metformin_GSE166707_516_directional_score_barplot_SCI.png/pdf\n")
cat("  07_plots/Metformin_GSE166707_516_score_heatmap_SCI.png/pdf\n")
cat("  08_summary/GSE46263_metformin_GSE166707_516_validation_report.txt\n")
cat("============================================================\n")












############################################################
# GSE15483_Fenofibrate_GSE166707_516_signature_validation.R
#
# Purpose:
#   Use external fenofibrate perturbation ranked lists from GSE15483
#   to validate GSE166707-derived 516 temporal core signatures.
#
# Three validation layers:
#
# 1. Overall 516 expression-direction reversal
#    Core516_UP   expected NES < 0
#    Core516_DOWN expected NES > 0
#
# 2. Residual non-viral expression-direction reversal
#    Residual_UP   expected NES < 0
#    Residual_DOWN expected NES > 0
#
# 3. Virus-host factor functional direction validation
#    HRF desired-up  expected NES > 0
#    HDF desired-down expected NES < 0
#
# Input:
#   D:/WeiRuan/NiV/非诺贝特GSE15483/GSE166707验证/516带表达方向带病毒标签.csv
#
# Ranked lists:
#   D:/WeiRuan/NiV/非诺贝特GSE15483/fenofibrate_reversal_validation_FINAL/03_ranked_gene_lists
#
# Output:
#   D:/WeiRuan/NiV/非诺贝特GSE15483/GSE166707验证
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 600)
set.seed(123)

############################################################
# 0. Paths and parameters
############################################################

work_dir <- "D:/WeiRuan/NiV/非诺贝特GSE15483/GSE166707验证"
setwd(work_dir)

signature_file <- "D:/WeiRuan/NiV/非诺贝特GSE15483/GSE166707验证/516带表达方向带病毒标签.csv"

ranked_dir <- "D:/WeiRuan/NiV/非诺贝特GSE15483/fenofibrate_reversal_validation_FINAL/03_ranked_gene_lists"

out_root <- "D:/WeiRuan/NiV/非诺贝特GSE15483/GSE166707验证"

min_gs_size <- 5
max_gs_size <- 5000

dir_list <- c(
  out_root,
  file.path(out_root, "00_signature_sets"),
  file.path(out_root, "01_GSEA_results"),
  file.path(out_root, "02_score_tables"),
  file.path(out_root, "03_plots"),
  file.path(out_root, "04_leading_edge"),
  file.path(out_root, "05_summary")
)

invisible(lapply(dir_list, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!file.exists(signature_file)) {
  stop("Cannot find signature file: ", signature_file)
}

if (!dir.exists(ranked_dir)) {
  stop("Cannot find ranked list directory: ", ranked_dir)
}

############################################################
# 1. Packages
############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "clusterProfiler",
  "pheatmap"
)

missing_packages <- required_packages[
  !sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  bioc_avail <- BiocManager::available()
  bioc_to_install <- intersect(missing_packages, bioc_avail)
  cran_to_install <- setdiff(missing_packages, bioc_to_install)
  
  if (length(bioc_to_install) > 0) {
    BiocManager::install(bioc_to_install, update = FALSE, ask = FALSE)
  }
  
  if (length(cran_to_install) > 0) {
    install.packages(cran_to_install)
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

############################################################
# 2. Helper functions
############################################################

write_csv_safe <- function(x, file) {
  write.csv(x, file, row.names = FALSE, quote = FALSE, na = "")
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_split(
    x,
    "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*",
    simplify = TRUE
  )[, 1]
  x <- stringr::str_trim(x)
  x[x == ""] <- NA_character_
  x[x == "NA"] <- NA_character_
  x[x == "NULL"] <- NA_character_
  x
}

theme_sci <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black")
    )
}

break_ties_deterministic <- function(x) {
  x + seq_along(x) * 1e-12
}

format_p_label <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

read_ranked_list <- function(file) {
  
  df <- read.csv(file, check.names = FALSE)
  
  if (!"gene_symbol" %in% colnames(df)) {
    stop("Ranked list missing gene_symbol column: ", file)
  }
  
  if ("t_ranked" %in% colnames(df)) {
    rank_col <- "t_ranked"
  } else if ("t" %in% colnames(df)) {
    rank_col <- "t"
  } else {
    stop("Ranked list missing t_ranked or t column: ", file)
  }
  
  out <- df %>%
    dplyr::mutate(
      gene_symbol = clean_symbol(gene_symbol),
      rank_value = as.numeric(.data[[rank_col]])
    ) %>%
    dplyr::filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      !is.na(rank_value),
      is.finite(rank_value)
    ) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::slice_max(order_by = abs(rank_value), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_value))
  
  geneList <- out$rank_value
  names(geneList) <- out$gene_symbol
  
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- break_ties_deterministic(geneList)
  
  list(
    df = out,
    geneList = geneList
  )
}

run_one_gsea <- function(geneList, term2gene, analysis_level, comparison) {
  
  term2gene_use <- term2gene %>%
    dplyr::filter(gene %in% names(geneList))
  
  if (nrow(term2gene_use) == 0) {
    return(NULL)
  }
  
  gsea_res <- tryCatch(
    clusterProfiler::GSEA(
      geneList = geneList,
      TERM2GENE = term2gene_use,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) {
      message("GSEA failed: ", analysis_level, " | ", comparison, " | ", e$message)
      NULL
    }
  )
  
  if (is.null(gsea_res)) {
    return(NULL)
  }
  
  gsea_df <- as.data.frame(gsea_res)
  
  if (nrow(gsea_df) == 0) {
    return(NULL)
  }
  
  gsea_df$analysis_level <- analysis_level
  gsea_df$comparison <- comparison
  
  gsea_df
}

get_expected_direction <- function(analysis_level, gene_set) {
  
  if (analysis_level %in% c("Overall_516_reversal", "Residual_nonviral_reversal")) {
    if (stringr::str_detect(gene_set, "_UP$")) {
      return("NES < 0")
    }
    if (stringr::str_detect(gene_set, "_DOWN$")) {
      return("NES > 0")
    }
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    if (gene_set == "HRF_desired_up") {
      return("NES > 0")
    }
    if (gene_set == "HDF_desired_down") {
      return("NES < 0")
    }
  }
  
  return(NA_character_)
}

is_direction_matched <- function(analysis_level, gene_set, NES) {
  
  if (is.na(NES)) return(FALSE)
  
  if (analysis_level %in% c("Overall_516_reversal", "Residual_nonviral_reversal")) {
    if (stringr::str_detect(gene_set, "_UP$")) {
      return(NES < 0)
    }
    if (stringr::str_detect(gene_set, "_DOWN$")) {
      return(NES > 0)
    }
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    if (gene_set == "HRF_desired_up") {
      return(NES > 0)
    }
    if (gene_set == "HDF_desired_down") {
      return(NES < 0)
    }
  }
  
  FALSE
}

score_one_level <- function(df_one) {
  
  analysis_level <- unique(df_one$analysis_level)
  comparison <- unique(df_one$comparison)
  
  if (length(analysis_level) != 1 || length(comparison) != 1) {
    stop("score_one_level requires one analysis_level and one comparison.")
  }
  
  get_nes <- function(term) {
    x <- df_one$NES[df_one$ID == term]
    if (length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }
  
  get_padj <- function(term) {
    x <- df_one$p.adjust[df_one$ID == term]
    if (length(x) == 0) return(NA_real_)
    as.numeric(x[1])
  }
  
  if (analysis_level == "Overall_516_reversal") {
    
    nes_up <- get_nes("Core516_UP")
    nes_down <- get_nes("Core516_DOWN")
    
    score <- (nes_down - nes_up) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "Core516_DOWN",
      term_negative_expected = "Core516_UP",
      NES_positive_expected = nes_down,
      NES_negative_expected = nes_up,
      FDR_positive_expected = get_padj("Core516_DOWN"),
      FDR_negative_expected = get_padj("Core516_UP"),
      score = score,
      both_direction_matched = (!is.na(nes_down) && nes_down > 0) &&
        (!is.na(nes_up) && nes_up < 0),
      score_name = "overall_reversal_score",
      interpretation = ifelse(
        (!is.na(nes_down) && nes_down > 0) && (!is.na(nes_up) && nes_up < 0),
        "Directional reversal supported",
        "Directional reversal not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  if (analysis_level == "Residual_nonviral_reversal") {
    
    nes_up <- get_nes("Residual_UP")
    nes_down <- get_nes("Residual_DOWN")
    
    score <- (nes_down - nes_up) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "Residual_DOWN",
      term_negative_expected = "Residual_UP",
      NES_positive_expected = nes_down,
      NES_negative_expected = nes_up,
      FDR_positive_expected = get_padj("Residual_DOWN"),
      FDR_negative_expected = get_padj("Residual_UP"),
      score = score,
      both_direction_matched = (!is.na(nes_down) && nes_down > 0) &&
        (!is.na(nes_up) && nes_up < 0),
      score_name = "residual_reversal_score",
      interpretation = ifelse(
        (!is.na(nes_down) && nes_down > 0) && (!is.na(nes_up) && nes_up < 0),
        "Residual signature reversal supported",
        "Residual signature reversal not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  if (analysis_level == "Virus_factor_functional_direction") {
    
    nes_hrf <- get_nes("HRF_desired_up")
    nes_hdf <- get_nes("HDF_desired_down")
    
    score <- (nes_hrf - nes_hdf) / 2
    
    return(data.frame(
      analysis_level = analysis_level,
      comparison = comparison,
      term_positive_expected = "HRF_desired_up",
      term_negative_expected = "HDF_desired_down",
      NES_positive_expected = nes_hrf,
      NES_negative_expected = nes_hdf,
      FDR_positive_expected = get_padj("HRF_desired_up"),
      FDR_negative_expected = get_padj("HDF_desired_down"),
      score = score,
      both_direction_matched = (!is.na(nes_hrf) && nes_hrf > 0) &&
        (!is.na(nes_hdf) && nes_hdf < 0),
      score_name = "virus_factor_modulation_score",
      interpretation = ifelse(
        (!is.na(nes_hrf) && nes_hrf > 0) && (!is.na(nes_hdf) && nes_hdf < 0),
        "Virus-host factor functional modulation supported",
        "Virus-host factor functional modulation not fully supported"
      ),
      stringsAsFactors = FALSE
    ))
  }
  
  stop("Unknown analysis_level: ", analysis_level)
}

############################################################
# 3. Load 516 signature table
############################################################

cat("============================================================\n")
cat("Loading GSE166707 516 signature table\n")
cat("File: ", signature_file, "\n", sep = "")
cat("============================================================\n")

sig_df <- read.csv(signature_file, check.names = FALSE)

required_cols <- c("Genes", "direction_type", "Host factor type")

missing_cols <- setdiff(required_cols, colnames(sig_df))

if (length(missing_cols) > 0) {
  stop("Signature file missing required columns: ", paste(missing_cols, collapse = ", "))
}

sig_df <- sig_df %>%
  dplyr::mutate(
    Genes = clean_symbol(Genes),
    direction_type = as.character(direction_type),
    direction_type = stringr::str_trim(direction_type),
    host_factor_type_raw = as.character(`Host factor type`),
    host_factor_type_raw = stringr::str_trim(host_factor_type_raw),
    host_factor_type = dplyr::case_when(
      toupper(host_factor_type_raw) == "HDF" ~ "HDF",
      toupper(host_factor_type_raw) == "HRF" ~ "HRF",
      TRUE ~ "Non_viral"
    )
  ) %>%
  dplyr::filter(!is.na(Genes), Genes != "")

sig_df <- sig_df %>%
  dplyr::distinct(Genes, .keep_all = TRUE)

############################################################
# 4. Build three layers of gene sets
############################################################

# Layer 1: Overall 516 expression-direction reversal
core_up <- sig_df %>%
  dplyr::filter(direction_type == "all_up") %>%
  dplyr::pull(Genes) %>%
  unique()

core_down <- sig_df %>%
  dplyr::filter(direction_type == "all_down") %>%
  dplyr::pull(Genes) %>%
  unique()

# Layer 2: Residual non-viral genes, expression-direction reversal
residual_up <- sig_df %>%
  dplyr::filter(host_factor_type == "Non_viral", direction_type == "all_up") %>%
  dplyr::pull(Genes) %>%
  unique()

residual_down <- sig_df %>%
  dplyr::filter(host_factor_type == "Non_viral", direction_type == "all_down") %>%
  dplyr::pull(Genes) %>%
  unique()

# Layer 3: Virus-host factor functional direction
# Your CMap design:
#   HRF was put into down genes -> negative CMap drug should upregulate HRF.
#   HDF was put into up genes   -> negative CMap drug should downregulate HDF.
hrf_desired_up <- sig_df %>%
  dplyr::filter(host_factor_type == "HRF") %>%
  dplyr::pull(Genes) %>%
  unique()

hdf_desired_down <- sig_df %>%
  dplyr::filter(host_factor_type == "HDF") %>%
  dplyr::pull(Genes) %>%
  unique()

signature_summary <- data.frame(
  analysis_level = c(
    "Overall_516_reversal",
    "Overall_516_reversal",
    "Residual_nonviral_reversal",
    "Residual_nonviral_reversal",
    "Virus_factor_functional_direction",
    "Virus_factor_functional_direction"
  ),
  gene_set = c(
    "Core516_UP",
    "Core516_DOWN",
    "Residual_UP",
    "Residual_DOWN",
    "HRF_desired_up",
    "HDF_desired_down"
  ),
  expected_direction = c(
    "NES < 0",
    "NES > 0",
    "NES < 0",
    "NES > 0",
    "NES > 0",
    "NES < 0"
  ),
  biological_meaning = c(
    "Drug suppresses GSE166707 core genes upregulated during NiV infection",
    "Drug restores GSE166707 core genes downregulated during NiV infection",
    "Drug suppresses non-viral-label residual upregulated host-response genes",
    "Drug restores non-viral-label residual downregulated host-response genes",
    "Drug upregulates host restriction / antiviral factors",
    "Drug downregulates host dependency / proviral factors"
  ),
  n_genes = c(
    length(core_up),
    length(core_down),
    length(residual_up),
    length(residual_down),
    length(hrf_desired_up),
    length(hdf_desired_down)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  signature_summary,
  file.path(out_root, "00_signature_sets", "GSE166707_516_signature_set_summary.csv")
)

write_csv_safe(data.frame(gene_symbol = core_up),
               file.path(out_root, "00_signature_sets", "Core516_UP.csv"))
write_csv_safe(data.frame(gene_symbol = core_down),
               file.path(out_root, "00_signature_sets", "Core516_DOWN.csv"))
write_csv_safe(data.frame(gene_symbol = residual_up),
               file.path(out_root, "00_signature_sets", "Residual_UP.csv"))
write_csv_safe(data.frame(gene_symbol = residual_down),
               file.path(out_root, "00_signature_sets", "Residual_DOWN.csv"))
write_csv_safe(data.frame(gene_symbol = hrf_desired_up),
               file.path(out_root, "00_signature_sets", "HRF_desired_up.csv"))
write_csv_safe(data.frame(gene_symbol = hdf_desired_down),
               file.path(out_root, "00_signature_sets", "HDF_desired_down.csv"))

print(signature_summary)

############################################################
# 5. Build TERM2GENE tables
############################################################

term2gene_overall <- data.frame(
  term = c(
    rep("Core516_UP", length(core_up)),
    rep("Core516_DOWN", length(core_down))
  ),
  gene = c(core_up, core_down),
  analysis_level = "Overall_516_reversal",
  stringsAsFactors = FALSE
)

term2gene_residual <- data.frame(
  term = c(
    rep("Residual_UP", length(residual_up)),
    rep("Residual_DOWN", length(residual_down))
  ),
  gene = c(residual_up, residual_down),
  analysis_level = "Residual_nonviral_reversal",
  stringsAsFactors = FALSE
)

term2gene_virus <- data.frame(
  term = c(
    rep("HRF_desired_up", length(hrf_desired_up)),
    rep("HDF_desired_down", length(hdf_desired_down))
  ),
  gene = c(hrf_desired_up, hdf_desired_down),
  analysis_level = "Virus_factor_functional_direction",
  stringsAsFactors = FALSE
)

term2gene_all <- dplyr::bind_rows(
  term2gene_overall,
  term2gene_residual,
  term2gene_virus
)

write_csv_safe(
  term2gene_all,
  file.path(out_root, "00_signature_sets", "TERM2GENE_all_signature_layers.csv")
)

############################################################
# 6. Load fenofibrate ranked lists
############################################################

cat("============================================================\n")
cat("Loading fenofibrate ranked lists\n")
cat("Directory: ", ranked_dir, "\n", sep = "")
cat("============================================================\n")

rank_files <- list.files(
  ranked_dir,
  pattern = "^ranked_gene_list_t_.*\\.csv$",
  full.names = TRUE
)

if (length(rank_files) == 0) {
  stop("No ranked_gene_list_t_*.csv found in: ", ranked_dir)
}

extract_comparison <- function(file) {
  x <- basename(file)
  x <- sub("^ranked_gene_list_t_", "", x)
  x <- sub("\\.csv$", "", x)
  x
}

rank_info <- data.frame(
  file = rank_files,
  comparison = sapply(rank_files, extract_comparison),
  stringsAsFactors = FALSE
)

expected_order <- c(
  "Feno_2h_vs_Untreated",
  "Feno_4h_vs_Untreated",
  "Feno_6h_vs_Untreated",
  "Feno_8h_vs_Untreated",
  "Feno_18h_vs_Untreated"
)

rank_info <- rank_info %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = expected_order)
  ) %>%
  dplyr::arrange(comparison)

write_csv_safe(
  rank_info,
  file.path(out_root, "05_summary", "ranked_list_files_used.csv")
)

print(rank_info)

############################################################
# 7. Run GSEA for all layers and all time points
############################################################

cat("============================================================\n")
cat("Running GSEA validation\n")
cat("============================================================\n")

gsea_all_list <- list()
mapped_summary_list <- list()

for (i in seq_len(nrow(rank_info))) {
  
  file <- rank_info$file[i]
  comparison <- as.character(rank_info$comparison[i])
  
  cat("Processing: ", comparison, "\n", sep = "")
  
  rank_obj <- read_ranked_list(file)
  geneList <- rank_obj$geneList
  rank_df <- rank_obj$df
  
  write_csv_safe(
    rank_df,
    file.path(out_root, "01_GSEA_results", paste0("ranked_list_used_", comparison, ".csv"))
  )
  
  # Mapping summary for each gene set in this ranked list
  mapped_summary <- term2gene_all %>%
    dplyr::group_by(analysis_level, term) %>%
    dplyr::summarise(
      n_original = dplyr::n_distinct(gene),
      n_mapped_to_ranked_list = length(intersect(unique(gene), names(geneList))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(comparison = comparison)
  
  mapped_summary_list[[comparison]] <- mapped_summary
  
  # Run three layers
  gsea_overall <- run_one_gsea(
    geneList = geneList,
    term2gene = term2gene_overall,
    analysis_level = "Overall_516_reversal",
    comparison = comparison
  )
  
  gsea_residual <- run_one_gsea(
    geneList = geneList,
    term2gene = term2gene_residual,
    analysis_level = "Residual_nonviral_reversal",
    comparison = comparison
  )
  
  gsea_virus <- run_one_gsea(
    geneList = geneList,
    term2gene = term2gene_virus,
    analysis_level = "Virus_factor_functional_direction",
    comparison = comparison
  )
  
  gsea_this <- dplyr::bind_rows(
    gsea_overall,
    gsea_residual,
    gsea_virus
  )
  
  if (!is.null(gsea_this) && nrow(gsea_this) > 0) {
    
    gsea_this <- gsea_this %>%
      dplyr::mutate(
        expected_direction = mapply(get_expected_direction, analysis_level, ID),
        direction_matched = mapply(is_direction_matched, analysis_level, ID, NES),
        direction_status = ifelse(direction_matched, "Matched", "Not matched")
      )
    
    gsea_all_list[[comparison]] <- gsea_this
    
    write_csv_safe(
      gsea_this,
      file.path(out_root, "01_GSEA_results", paste0("GSEA_all_layers_", comparison, ".csv"))
    )
  }
}

gsea_all <- dplyr::bind_rows(gsea_all_list)
mapped_summary_all <- dplyr::bind_rows(mapped_summary_list)

if (nrow(gsea_all) == 0) {
  stop("No GSEA results generated. Please check gene sets and ranked lists.")
}

write_csv_safe(
  gsea_all,
  file.path(out_root, "01_GSEA_results", "GSEA_all_layers_all_timepoints_long.csv")
)

write_csv_safe(
  mapped_summary_all,
  file.path(out_root, "00_signature_sets", "signature_mapping_to_ranked_lists_summary.csv")
)

############################################################
# 8. Compute scores
############################################################

score_table <- gsea_all %>%
  dplyr::group_by(analysis_level, comparison) %>%
  dplyr::group_split() %>%
  lapply(score_one_level) %>%
  dplyr::bind_rows()

score_table <- score_table %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = expected_order),
    analysis_level = factor(
      analysis_level,
      levels = c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )
    )
  ) %>%
  dplyr::arrange(analysis_level, comparison)

write_csv_safe(
  score_table,
  file.path(out_root, "02_score_tables", "Fenofibrate_GSE166707_516_validation_score_table.csv")
)

print(score_table)

############################################################
# 9. Leading-edge genes
############################################################

leading_edge_df <- gsea_all %>%
  dplyr::select(
    analysis_level,
    comparison,
    ID,
    NES,
    pvalue,
    p.adjust,
    expected_direction,
    direction_matched,
    core_enrichment
  ) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(gene_symbol = core_enrichment) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

write_csv_safe(
  leading_edge_df,
  file.path(out_root, "04_leading_edge", "GSE166707_516_fenofibrate_GSEA_leading_edge_genes.csv")
)

leading_edge_count <- leading_edge_df %>%
  dplyr::count(analysis_level, ID, gene_symbol, sort = TRUE)

write_csv_safe(
  leading_edge_count,
  file.path(out_root, "04_leading_edge", "GSE166707_516_fenofibrate_leading_edge_gene_frequency.csv")
)

############################################################
# 10. SCI-style plots
############################################################

cat("============================================================\n")
cat("Plotting SCI-style figures\n")
cat("============================================================\n")

comparison_label_map <- c(
  "Feno_2h_vs_Untreated" = "Fenofibrate 2 h",
  "Feno_4h_vs_Untreated" = "Fenofibrate 4 h",
  "Feno_6h_vs_Untreated" = "Fenofibrate 6 h",
  "Feno_8h_vs_Untreated" = "Fenofibrate 8 h",
  "Feno_18h_vs_Untreated" = "Fenofibrate 18 h"
)

analysis_label_map <- c(
  "Overall_516_reversal" = "Overall 516 reversal",
  "Residual_nonviral_reversal" = "Residual non-viral reversal",
  "Virus_factor_functional_direction" = "Virus-host factor direction"
)

gene_set_label_map <- c(
  "Core516_UP" = "Core516 UP",
  "Core516_DOWN" = "Core516 DOWN",
  "Residual_UP" = "Residual UP",
  "Residual_DOWN" = "Residual DOWN",
  "HRF_desired_up" = "HRF desired-up",
  "HDF_desired_down" = "HDF desired-down"
)

plot_df <- gsea_all %>%
  dplyr::mutate(
    comparison = factor(comparison, levels = expected_order),
    comparison_label = factor(
      comparison_label_map[as.character(comparison)],
      levels = comparison_label_map[expected_order]
    ),
    analysis_label = factor(
      analysis_label_map[analysis_level],
      levels = analysis_label_map[c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )]
    ),
    gene_set_label = gene_set_label_map[ID],
    neg_log10_fdr = -log10(pmax(p.adjust, 1e-300)),
    direction_status = ifelse(direction_matched, "Matched", "Not matched")
  )

p_dot <- ggplot(
  plot_df,
  aes(x = NES, y = comparison_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_point(
    aes(
      color = gene_set_label,
      shape = direction_status,
      size = neg_log10_fdr
    ),
    alpha = 0.92,
    stroke = 0.9
  ) +
  facet_grid(analysis_label ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c("Matched" = 16, "Not matched" = 1)) +
  scale_size_continuous(range = c(2.8, 7.2)) +
  labs(
    title = "Fenofibrate perturbation against GSE166707-derived NiV signatures",
    subtitle = "External GSE15483 ranked lists tested by GSEA",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    color = NULL,
    shape = NULL,
    size = expression(-log[10](FDR))
  ) +
  theme_sci(11) +
  theme(
    legend.position = "right",
    strip.text.y = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold")
  )

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_GSEA_NES_dotplot_SCI.png"),
  p_dot,
  width = 10.2,
  height = 8.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_GSEA_NES_dotplot_SCI.pdf"),
  p_dot,
  width = 10.2,
  height = 8.8
)

score_plot_df <- score_table %>%
  dplyr::mutate(
    comparison_label = factor(
      comparison_label_map[as.character(comparison)],
      levels = comparison_label_map[expected_order]
    ),
    analysis_label = factor(
      analysis_label_map[as.character(analysis_level)],
      levels = analysis_label_map[c(
        "Overall_516_reversal",
        "Residual_nonviral_reversal",
        "Virus_factor_functional_direction"
      )]
    ),
    match_status = ifelse(both_direction_matched, "Both directions matched", "Partial or no match")
  )

p_score <- ggplot(
  score_plot_df,
  aes(x = comparison_label, y = score, fill = match_status)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.2f", score)),
    vjust = ifelse(score_plot_df$score >= 0, -0.35, 1.25),
    size = 3.2
  ) +
  facet_wrap(~ analysis_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Both directions matched" = "#5AAE61",
      "Partial or no match" = "grey72"
    )
  ) +
  labs(
    title = "Fenofibrate signature-level directional scores",
    subtitle = "Reversal scores for core/residual signatures; functional modulation score for virus-host factors",
    x = NULL,
    y = "Directional score",
    fill = NULL
  ) +
  theme_sci(11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_directional_score_barplot_SCI.png"),
  p_score,
  width = 9.5,
  height = 9.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_directional_score_barplot_SCI.pdf"),
  p_score,
  width = 9.5,
  height = 9.5
)

# Heatmap-style score overview
p_heat <- ggplot(
  score_plot_df,
  aes(x = comparison_label, y = analysis_label, fill = score)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3.6, color = "black") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0
  ) +
  labs(
    title = "Directional score overview",
    subtitle = "Positive scores indicate expected reversal/modulation direction",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_sci(11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_score_heatmap_SCI.png"),
  p_heat,
  width = 9.2,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_root, "03_plots", "Fenofibrate_GSE166707_516_score_heatmap_SCI.pdf"),
  p_heat,
  width = 9.2,
  height = 4.8
)

############################################################
# 11. Summary report
############################################################

analysis_report <- c(
  "GSE166707 516-core signature validation using GSE15483 fenofibrate perturbation",
  "============================================================",
  paste0("Signature file: ", signature_file),
  paste0("Ranked list directory: ", ranked_dir),
  paste0("Output root: ", out_root),
  "",
  "Conceptual design:",
  "This analysis uses external GSE15483 fenofibrate-treated HUVEC ranked lists",
  "to test GSE166707-derived NiV temporal core signatures.",
  "",
  "Three validation layers:",
  "",
  "1. Overall 516 expression-direction reversal",
  "   Core516_UP expected NES < 0",
  "   Core516_DOWN expected NES > 0",
  "   Score = (NES_Core516_DOWN - NES_Core516_UP) / 2",
  "",
  "2. Residual non-viral expression-direction reversal",
  "   Residual_UP expected NES < 0",
  "   Residual_DOWN expected NES > 0",
  "   Score = (NES_Residual_DOWN - NES_Residual_UP) / 2",
  "",
  "3. Virus-host factor functional direction validation",
  "   HRF_desired_up expected NES > 0",
  "   HDF_desired_down expected NES < 0",
  "   Score = (NES_HRF_desired_up - NES_HDF_desired_down) / 2",
  "",
  "Why HRF/HDF use functional directions:",
  "In the original CMap prediction, HRFs were placed into down genes and HDFs were placed into up genes.",
  "For negative connectivity drugs, this implies that the expected drug effect is to upregulate HRFs",
  "and downregulate HDFs. Therefore, validation of HRF/HDF should follow this predefined functional direction",
  "rather than the original infection logFC direction.",
  "",
  "Main output tables:",
  "00_signature_sets/GSE166707_516_signature_set_summary.csv",
  "01_GSEA_results/GSEA_all_layers_all_timepoints_long.csv",
  "02_score_tables/Fenofibrate_GSE166707_516_validation_score_table.csv",
  "04_leading_edge/GSE166707_516_fenofibrate_GSEA_leading_edge_genes.csv",
  "",
  "Main figures:",
  "03_plots/Fenofibrate_GSE166707_516_GSEA_NES_dotplot_SCI.pdf/png",
  "03_plots/Fenofibrate_GSE166707_516_directional_score_barplot_SCI.pdf/png",
  "03_plots/Fenofibrate_GSE166707_516_score_heatmap_SCI.pdf/png"
)

writeLines(
  analysis_report,
  file.path(out_root, "05_summary", "GSE166707_516_fenofibrate_validation_report.txt")
)

write_csv_safe(
  signature_summary,
  file.path(out_root, "05_summary", "GSE166707_516_signature_set_summary.csv")
)

write_csv_safe(
  mapped_summary_all,
  file.path(out_root, "05_summary", "GSE166707_516_signature_mapping_summary.csv")
)

write_csv_safe(
  score_table,
  file.path(out_root, "05_summary", "GSE166707_516_fenofibrate_validation_score_table.csv")
)

cat("============================================================\n")
cat("GSE166707 516-core signature validation finished.\n")
cat("Output root:\n")
cat(out_root, "\n\n")
cat("Key outputs:\n")
cat("  00_signature_sets/GSE166707_516_signature_set_summary.csv\n")
cat("  01_GSEA_results/GSEA_all_layers_all_timepoints_long.csv\n")
cat("  02_score_tables/Fenofibrate_GSE166707_516_validation_score_table.csv\n")
cat("  03_plots/Fenofibrate_GSE166707_516_GSEA_NES_dotplot_SCI.png/pdf\n")
cat("  03_plots/Fenofibrate_GSE166707_516_directional_score_barplot_SCI.png/pdf\n")
cat("  03_plots/Fenofibrate_GSE166707_516_score_heatmap_SCI.png/pdf\n")
cat("  04_leading_edge/GSE166707_516_fenofibrate_GSEA_leading_edge_genes.csv\n")
cat("  05_summary/GSE166707_516_fenofibrate_validation_report.txt\n")
cat("============================================================\n")