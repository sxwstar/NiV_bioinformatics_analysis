############################################################
# GSE166707 – FINAL FULL PIPELINE (TIMEPOINT-BASED, MICROARRAY)
#
# Design:
#   4h/8h/12h/16h timepoint-specific limma
#   Infection vs Uninfected (5 vs 5 at each timepoint)
#
# DEG cutoff:
#   P.Value < 0.05
#   |logFC|  > 0.58
#
# Enrichment cutoff:
#   P.Value < 0.05
#
# Includes:
#   00_QC
#   01_PCA
#   02_DEG_timepoint
#   03_GO_KEGG_timepoint
#   04_GSEA_timepoint
#   Figure_1
#   Figure_2
#   Figure_3
#
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

############################
# 0) Parameters & paths
############################
deg_p_cutoff  <- 0.05
deg_fc_cutoff <- 0.58

enrich_p_cutoff <- 0.05  # ORA enrichment retained by raw P.Value < 0.05
gsea_p_cutoff   <- 0.05  # GSEA retained by raw P.Value < 0.05

expr_file  <- "D:/WeiRuan/NiV/GSE166707/counts.xlsx"
probe_file <- "D:/WeiRuan/NiV/GSE166707/probe.xlsx"
out_root   <- "D:/WeiRuan/NiV/GSE166707/results_FINAL_timepoint_FULL_SCI"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

subdirs <- c(
  "00_QC",
  "01_PCA",
  "02_DEG_4h","02_DEG_8h","02_DEG_12h","02_DEG_16h",
  "03_GO_KEGG_4h","03_GO_KEGG_8h","03_GO_KEGG_12h","03_GO_KEGG_16h",
  "04_GSEA_4h","04_GSEA_8h","04_GSEA_12h","04_GSEA_16h",
  "Figure_1",
  "Figure_2",
  "Figure_3"
)
invisible(lapply(file.path(out_root, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE))

cat("============================================================\n")
cat("GSE166707 MICROARRAY TIMEPOINT PIPELINE START\n")
cat("DEG cutoff: P.Value <", deg_p_cutoff, " & |logFC| >", deg_fc_cutoff, "\n")
cat("Enrichment cutoff: P.Value <", enrich_p_cutoff, "\n")
cat("Output root:", out_root, "\n")
cat("============================================================\n\n")

############################
# 1) Packages
############################
need_pkgs <- c(
  "readxl","stringr","dplyr","tibble","tidyr",
  "limma","ggplot2","pheatmap","RColorBrewer","patchwork",
  "clusterProfiler","org.Hs.eg.db","enrichplot",
  "grid","png"
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
  library(readxl)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(grid)
  library(png)
})


############################
# 2) Sample information
############################
inf_4h  <- c("GSM5079393","GSM5079394","GSM5079395","GSM5079396","GSM5079397")
inf_8h  <- c("GSM5079398","GSM5079399","GSM5079400","GSM5079401","GSM5079402")
inf_12h <- c("GSM5079403","GSM5079404","GSM5079405","GSM5079406","GSM5079407")
inf_16h <- c("GSM5079408","GSM5079409","GSM5079410","GSM5079411","GSM5079412")

ctl_4h  <- c("GSM5079413","GSM5079414","GSM5079415","GSM5079416","GSM5079417")
ctl_8h  <- c("GSM5079418","GSM5079419","GSM5079420","GSM5079421","GSM5079422")
ctl_12h <- c("GSM5079423","GSM5079424","GSM5079425","GSM5079426","GSM5079427")
ctl_16h <- c("GSM5079428","GSM5079429","GSM5079430","GSM5079431","GSM5079432")

time_list <- list(
  `4h`  = list(inf = inf_4h,  ctl = ctl_4h),
  `8h`  = list(inf = inf_8h,  ctl = ctl_8h),
  `12h` = list(inf = inf_12h, ctl = ctl_12h),
  `16h` = list(inf = inf_16h, ctl = ctl_16h)
)

all_samples <- c(inf_4h,inf_8h,inf_12h,inf_16h, ctl_4h,ctl_8h,ctl_12h,ctl_16h)
stopifnot(length(all_samples) == 40)

sample_order <- c(
  ctl_4h,  inf_4h,
  ctl_8h,  inf_8h,
  ctl_12h, inf_12h,
  ctl_16h, inf_16h
)

############################
# 3) Utility functions
############################
theme_sci <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

group_cols <- c("Uninfected" = "#3C5488FF", "Infection" = "#D73027")
time_cols  <- c("4h" = "#1B9E77", "8h" = "#D95F02", "12h" = "#7570B3", "16h" = "#E7298A")
deg_cols   <- c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")

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

break_ties <- function(x) x + seq_along(x) * 1e-12

safe_save_plot <- function(p, filename, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = dpi, bg = "white")
}

get_gaps_by_time <- function(samples, sample_info) {
  if (length(samples) <= 1) return(NULL)
  times <- as.character(sample_info[samples, "time"])
  r <- rle(times)
  gaps <- cumsum(r$lengths)
  gaps <- gaps[-length(gaps)]
  if (length(gaps) == 0) return(NULL)
  gaps
}

subset_enrich_obj <- function(obj, df_keep) {
  obj2 <- obj
  obj2@result <- df_keep
  obj2
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

safe_dotplot <- function(obj, file, title, showCategory = 15) {
  df0 <- as.data.frame(obj)
  if (nrow(df0) == 0) {
    writeLines("dotplot skipped: empty enrichment result.", paste0(file, "_DOTPLOT_EMPTY.txt"))
    return(invisible(FALSE))
  }
  show_n <- min(showCategory, nrow(df0))
  p <- tryCatch(
    enrichplot::dotplot(obj, showCategory = show_n) +
      ggtitle(title) +
      theme_sci(12),
    error = function(e) e
  )
  if (inherits(p, "error")) {
    writeLines(c("dotplot failed", p$message), paste0(file, "_DOTPLOT_ERROR.txt"))
    return(invisible(FALSE))
  }
  safe_save_plot(p, file, width = 10, height = 7, dpi = 320)
  invisible(TRUE)
}

safe_barplot <- function(obj, file, title, showCategory = 15) {
  df0 <- as.data.frame(obj)
  if (nrow(df0) == 0) {
    writeLines("barplot skipped: empty enrichment result.", paste0(file, "_BARPLOT_EMPTY.txt"))
    return(invisible(FALSE))
  }
  show_n <- min(showCategory, nrow(df0))
  p <- tryCatch(
    clusterProfiler::barplot(obj, showCategory = show_n) +
      ggtitle(title) +
      theme_sci(12),
    error = function(e) e
  )
  if (inherits(p, "error")) {
    writeLines(c("barplot failed", p$message), paste0(file, "_BARPLOT_ERROR.txt"))
    return(invisible(FALSE))
  }
  safe_save_plot(p, file, width = 10, height = 7, dpi = 320)
  invisible(TRUE)
}

safe_gseaplot2 <- function(gsea_obj, geneSetID, file, title = NULL, width = 8, height = 6) {
  p <- tryCatch({
    enrichplot::gseaplot2(
      gsea_obj,
      geneSetID = geneSetID,
      title = title,
      base_size = 12,
      pvalue_table = TRUE
    )
  }, error = function(e) e)
  
  if (inherits(p, "error")) {
    writeLines(c("gseaplot2 failed", p$message), paste0(file, "_GSEAPLOT2_ERROR.txt"))
    return(invisible(FALSE))
  }
  
  ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = 320, bg = "white")
  invisible(TRUE)
}

make_volcano <- function(deg_df, out_file, title_txt) {
  df <- deg_df
  df$gene <- rownames(df)
  df$group <- "NS"
  df$group[df$P.Value < deg_p_cutoff & df$logFC >  deg_fc_cutoff] <- "Up"
  df$group[df$P.Value < deg_p_cutoff & df$logFC < -deg_fc_cutoff] <- "Down"
  df$negLogP <- -log10(df$P.Value)
  
  p <- ggplot(df, aes(x = logFC, y = negLogP, color = group)) +
    geom_point(size = 1.6, alpha = 0.75) +
    geom_vline(xintercept = c(-deg_fc_cutoff, deg_fc_cutoff), linetype = 2, linewidth = 0.5) +
    geom_hline(yintercept = -log10(deg_p_cutoff), linetype = 2, linewidth = 0.5) +
    scale_color_manual(values = deg_cols) +
    labs(
      title = title_txt,
      x = "log2 Fold Change",
      y = expression(-log[10](italic(P)~value)),
      color = NULL
    ) +
    theme_sci(12)
  safe_save_plot(p, out_file, width = 6, height = 5, dpi = 320)
  p
}

make_heatmap_deg <- function(expr_mat, genes_use, sample_info, out_file, title_txt, cluster_cols = FALSE) {
  genes_use <- unique(genes_use)
  genes_use <- genes_use[genes_use %in% rownames(expr_mat)]
  if (length(genes_use) < 2) {
    writeLines("Too few genes for heatmap.", paste0(out_file, "_SKIPPED.txt"))
    return(invisible(FALSE))
  }
  
  if (length(genes_use) > 100) {
    genes_use <- genes_use[1:100]
  }
  
  sample_keep <- intersect(sample_order, colnames(expr_mat))
  hm <- expr_mat[genes_use, sample_keep, drop = FALSE]
  hm <- t(scale(t(hm)))
  hm[!is.finite(hm)] <- 0
  
  ann_col <- data.frame(
    Group = sample_info[sample_keep, "group"],
    Time  = sample_info[sample_keep, "time"]
  )
  rownames(ann_col) <- sample_keep
  
  ann_colors <- list(
    Group = group_cols,
    Time  = time_cols
  )
  
  gaps <- get_gaps_by_time(sample_keep, sample_info)
  
  pheatmap::pheatmap(
    hm,
    color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_cols = cluster_cols,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    border_color = NA,
    fontsize = 10,
    fontsize_row = 7,
    main = title_txt,
    gaps_col = gaps,
    filename = out_file,
    width = 9,
    height = 10
  )
  invisible(TRUE)
}

read_png_as_plot <- function(img_path) {
  if (!file.exists(img_path)) {
    return(ggplot() + theme_void() + labs(title = paste("Missing:", basename(img_path))))
  }
  img <- png::readPNG(img_path)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void()
}

convert_gene_ratio <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.character(x)) {
    return(sapply(x, function(y) {
      if (grepl("/", y, fixed = TRUE)) {
        eval(parse(text = y))
      } else {
        suppressWarnings(as.numeric(y))
      }
    }))
  }
  suppressWarnings(as.numeric(x))
}

make_enrichment_bubble <- function(df, out_file, title_txt, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) {
    writeLines("Empty enrichment table.", paste0(out_file, "_EMPTY.txt"))
    return(NULL)
  }
  
  df2 <- df[1:min(top_n, nrow(df)), , drop = FALSE]
  df2$GeneRatio_num <- convert_gene_ratio(df2$GeneRatio)
  df2$Description <- factor(df2$Description, levels = rev(df2$Description))
  
  p <- ggplot(df2, aes(x = GeneRatio_num, y = Description, size = Count, color = pvalue)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "reverse") +
    labs(
      title = title_txt,
      x = "Gene Ratio",
      y = NULL,
      size = "Count",
      color = "P.Value"
    ) +
    theme_sci(12)
  
  safe_save_plot(p, out_file, width = 8.5, height = 6.5, dpi = 320)
  p
}

make_gsea_nes_plot <- function(df, out_file, title_txt, top_n = 10) {
  if (is.null(df) || nrow(df) == 0) {
    writeLines("Empty GSEA table.", paste0(out_file, "_EMPTY.txt"))
    return(NULL)
  }
  
  df2 <- df[order(df$pvalue, -abs(df$NES)), , drop = FALSE]
  up_df   <- df2[df2$NES > 0, , drop = FALSE]
  down_df <- df2[df2$NES < 0, , drop = FALSE]
  
  up_df   <- head(up_df, top_n)
  down_df <- head(down_df, top_n)
  
  plot_df <- rbind(up_df, down_df)
  if (nrow(plot_df) == 0) {
    writeLines("No NES terms to plot.", paste0(out_file, "_EMPTY.txt"))
    return(NULL)
  }
  
  plot_df$Direction <- ifelse(plot_df$NES > 0, "Activated", "Suppressed")
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description[order(plot_df$NES)])
  
  p <- ggplot(plot_df, aes(x = NES, y = Description, fill = Direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Activated" = "#D73027", "Suppressed" = "#4575B4")) +
    labs(
      title = title_txt,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = NULL
    ) +
    theme_sci(12)
  
  safe_save_plot(p, out_file, width = 8.5, height = 6.5, dpi = 320)
  p
}

pick_latest_nonempty <- function(lst) {
  ord <- c("16h","12h","8h","4h")
  for (x in ord) {
    if (!is.null(lst[[x]]) && nrow(lst[[x]]) > 0) return(lst[[x]])
  }
  NULL
}

############################
# 4) Read expression & probe annotation
############################
cat("=== [00] Read counts.xlsx and probe.xlsx ===\n")

if (!file.exists(expr_file))  stop("Missing expr_file: ", expr_file)
if (!file.exists(probe_file)) stop("Missing probe_file: ", probe_file)

expr_df <- readxl::read_excel(expr_file)
if (!"ID" %in% colnames(expr_df)) stop("counts.xlsx must contain column 'ID'")

miss_cols <- setdiff(all_samples, colnames(expr_df))
if (length(miss_cols) > 0) stop("counts.xlsx missing sample columns: ", paste(miss_cols, collapse = ", "))

for (cc in all_samples) expr_df[[cc]] <- suppressWarnings(as.numeric(expr_df[[cc]]))

expr_mat0 <- as.matrix(expr_df[, all_samples])
rownames(expr_mat0) <- as.character(expr_df$ID)
storage.mode(expr_mat0) <- "numeric"

probe_df <- readxl::read_excel(probe_file)
if (!all(c("ID","Gene symbol") %in% colnames(probe_df))) {
  stop("probe.xlsx must contain columns: ID and 'Gene symbol'")
}

anno <- data.frame(
  ID = as.character(probe_df$ID),
  gene_symbol = as.character(probe_df[["Gene symbol"]]),
  stringsAsFactors = FALSE
)
anno$gene_symbol <- stringr::str_trim(anno$gene_symbol)
anno$gene_symbol <- stringr::str_split(
  anno$gene_symbol,
  "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*",
  simplify = TRUE
)[,1]
anno$gene_symbol <- stringr::str_trim(anno$gene_symbol)

expr_mat <- expr_mat0
expr_mat[expr_mat == -9999] <- NA

keep_probe_basic <- rowSums(!is.na(expr_mat)) >= 3
expr_mat <- expr_mat[keep_probe_basic, , drop = FALSE]

common_ids <- intersect(rownames(expr_mat), anno$ID)
expr_mat <- expr_mat[common_ids, , drop = FALSE]
anno2 <- anno[match(common_ids, anno$ID), , drop = FALSE]

ok_gene <- !is.na(anno2$gene_symbol) & anno2$gene_symbol != "" & anno2$gene_symbol != "NA"
expr_mat <- expr_mat[ok_gene, , drop = FALSE]
anno2 <- anno2[ok_gene, , drop = FALSE]

expr_mat <- expr_mat[, sample_order, drop = FALSE]

qc_note <- c(
  paste0("Raw probes: ", nrow(expr_mat0)),
  paste0("After filtering and annotation: ", nrow(expr_mat)),
  paste0("Samples: ", ncol(expr_mat)),
  "Microarray data were used directly for downstream analysis."
)
writeLines(qc_note, file.path(out_root, "00_QC", "QC_preprocess_notes.txt"))

df_long <- data.frame(
  value = as.vector(expr_mat),
  sample = rep(colnames(expr_mat), each = nrow(expr_mat))
)
df_long$sample <- factor(df_long$sample, levels = sample_order)

p_box <- ggplot(df_long, aes(x = sample, y = value, fill = sample)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.35) +
  scale_fill_manual(values = rep(c("#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF"), length.out = length(unique(df_long$sample)))) +
  labs(title = "Expression Distribution Across Samples", x = NULL, y = "Expression Intensity") +
  theme_sci(12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), legend.position = "none")
safe_save_plot(p_box, file.path(out_root, "00_QC", "Boxplot_expr.png"), width = 15, height = 6, dpi = 320)

############################
# 5) Collapse probe to gene
############################
cat("=== [00] Collapse probes to genes ===\n")

probe_mean <- rowMeans(expr_mat, na.rm = TRUE)

tmp <- data.frame(
  probe_id = rownames(expr_mat),
  gene_symbol = anno2$gene_symbol,
  mean_expr = probe_mean,
  stringsAsFactors = FALSE
)
tmp <- tmp[order(tmp$gene_symbol, -tmp$mean_expr), ]
tmp_uniq <- tmp[!duplicated(tmp$gene_symbol), ]

expr_gene <- expr_mat[tmp_uniq$probe_id, , drop = FALSE]
rownames(expr_gene) <- tmp_uniq$gene_symbol
storage.mode(expr_gene) <- "numeric"
expr_gene <- expr_gene[, sample_order, drop = FALSE]

write_table_dual(
  data.frame(gene_symbol = rownames(expr_gene), expr_gene, check.names = FALSE),
  file.path(out_root, "00_QC", "expr_gene_matrix.csv"),
  row.names = FALSE
)

############################
# 6) Sample metadata
############################
sample_info_all <- data.frame(
  sample = all_samples,
  group  = rep(NA_character_, length(all_samples)),
  time   = rep(NA_character_, length(all_samples)),
  stringsAsFactors = FALSE
)

for (tt in names(time_list)) {
  sample_info_all$group[sample_info_all$sample %in% time_list[[tt]]$inf] <- "Infection"
  sample_info_all$group[sample_info_all$sample %in% time_list[[tt]]$ctl] <- "Uninfected"
  sample_info_all$time[sample_info_all$sample %in% c(time_list[[tt]]$inf, time_list[[tt]]$ctl)] <- tt
}

sample_info_all$group <- factor(sample_info_all$group, levels = c("Uninfected","Infection"))
sample_info_all$time  <- factor(sample_info_all$time, levels = c("4h","8h","12h","16h"))
rownames(sample_info_all) <- sample_info_all$sample
sample_info_all <- sample_info_all[sample_order, , drop = FALSE]

write_table_dual(
  sample_info_all,
  file.path(out_root, "00_QC", "sample_info_all.csv"),
  row.names = FALSE
)

############################
# 7) PCA and sample correlation
############################
cat("=== [01] PCA ===\n")

expr_pca <- impute_row_median(expr_gene)
expr_pca <- expr_pca[, sample_order, drop = FALSE]

pca <- prcomp(t(expr_pca), scale. = FALSE)
var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = rownames(pca$x),
  group = sample_info_all[rownames(pca$x), "group"],
  time  = sample_info_all[rownames(pca$x), "time"],
  stringsAsFactors = FALSE
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = time)) +
  geom_point(size = 4, alpha = 0.9, stroke = 1.1) +
  stat_ellipse(aes(group = interaction(group, time)), linetype = 2, linewidth = 0.5, alpha = 0.4) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_sci(13) +
  theme(legend.position = "right")

p_pca_noellipse <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = time)) +
  geom_point(size = 4, alpha = 0.9, stroke = 1.1) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_sci(13) +
  theme(legend.position = "right")

safe_save_plot(p_pca, file.path(out_root, "01_PCA", "PCA_PC1_PC2_with_ellipse.png"), width = 8.5, height = 6.8, dpi = 320)
safe_save_plot(p_pca, file.path(out_root, "01_PCA", "PCA_PC1_PC2_with_ellipse.pdf"), width = 8.5, height = 6.8, dpi = 320)
safe_save_plot(p_pca_noellipse, file.path(out_root, "01_PCA", "PCA_PC1_PC2_no_ellipse.png"), width = 8.5, height = 6.8, dpi = 320)
safe_save_plot(p_pca_noellipse, file.path(out_root, "01_PCA", "PCA_PC1_PC2_no_ellipse.pdf"), width = 8.5, height = 6.8, dpi = 320)

cor_mat <- cor(expr_pca, use = "pairwise.complete.obs")
cor_mat <- cor_mat[sample_order, sample_order]

ann <- data.frame(
  Time  = sample_info_all[sample_order, "time"],
  Group = sample_info_all[sample_order, "group"]
)
rownames(ann) <- sample_order

ann_colors <- list(
  Group = group_cols,
  Time  = time_cols
)

gaps_cor <- get_gaps_by_time(sample_order, sample_info_all)

pheatmap::pheatmap(
  cor_mat,
  annotation_col = ann,
  annotation_row = ann,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_cor,
  gaps_col = gaps_cor,
  color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100),
  border_color = NA,
  fontsize = 9,
  main = "Sample-to-sample correlation based on global gene expression profiles",
  filename = file.path(out_root, "01_PCA", "Sample_correlation_heatmap.png"),
  width = 11,
  height = 10
)

############################
# 8) DEG + GO/KEGG + GSEA
############################
cat("=== [02-04] Timepoint-specific analyses ===\n")

deg_union <- character(0)
deg_all_list <- list()
deg_summary_list <- list()
go_top_list <- list()
kegg_top_list <- list()
gsea_go_top_list <- list()
gsea_kegg_top_list <- list()

run_timepoint_limma <- function(tt) {
  inf <- time_list[[tt]]$inf
  ctl <- time_list[[tt]]$ctl
  
  expr_sub <- expr_gene[, c(ctl, inf), drop = FALSE]
  
  meta <- data.frame(
    sample = c(ctl, inf),
    group = factor(c(rep("Uninfected", length(ctl)), rep("Infection", length(inf))),
                   levels = c("Uninfected","Infection")),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  design <- model.matrix(~ group, data = meta)
  fit <- limma::lmFit(expr_sub, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  
  tt_all <- limma::topTable(fit, coef = "groupInfection", number = Inf, sort.by = "P")
  tt_all <- tt_all[!is.na(tt_all$P.Value) & !is.na(tt_all$logFC), , drop = FALSE]
  deg_sig <- tt_all[tt_all$P.Value < deg_p_cutoff & abs(tt_all$logFC) > deg_fc_cutoff, , drop = FALSE]
  
  list(tt = tt_all, sig = deg_sig, fit = fit, expr_sub = expr_sub, meta = meta)
}

for (tt in names(time_list)) {
  cat("  -> ", tt, "\n", sep = "")
  
  res <- run_timepoint_limma(tt)
  deg_all_list[[tt]] <- res$tt
  deg_union <- union(deg_union, rownames(res$sig))
  
  out_deg_dir    <- file.path(out_root, paste0("02_DEG_", tt))
  out_enrich_dir <- file.path(out_root, paste0("03_GO_KEGG_", tt))
  out_gsea_dir   <- file.path(out_root, paste0("04_GSEA_", tt))
  
  write_table_dual(res$tt,  file.path(out_deg_dir, "DEG_all.csv"), row.names = TRUE)
  write_table_dual(res$sig, file.path(out_deg_dir, "DEG_sig_P0.05_logFC0.58.csv"), row.names = TRUE)
  
  deg_summary_list[[tt]] <- data.frame(
    Time = tt,
    Up = sum(res$sig$logFC > 0),
    Down = sum(res$sig$logFC < 0),
    Total = nrow(res$sig),
    stringsAsFactors = FALSE
  )
  
  writeLines(
    c(
      paste0("Timepoint: ", tt),
      "Comparison: Infection vs Uninfected",
      paste0("DEG_all n = ", nrow(res$tt)),
      paste0("DEG_sig n = ", nrow(res$sig)),
      paste0("Up = ", sum(res$sig$logFC > 0)),
      paste0("Down = ", sum(res$sig$logFC < 0))
    ),
    file.path(out_deg_dir, "summary.txt")
  )
  
  make_volcano(
    res$tt,
    file.path(out_deg_dir, paste0("Volcano_", tt, ".png")),
    paste0("Volcano Plot (", tt, " post infection)")
  )
  
  if (nrow(res$sig) >= 2) {
    top_genes <- rownames(res$sig)[order(res$sig$P.Value)][1:min(50, nrow(res$sig))]
    sample_info_sub <- sample_info_all[colnames(res$expr_sub), , drop = FALSE]
    make_heatmap_deg(
      res$expr_sub,
      top_genes,
      sample_info_sub,
      file.path(out_deg_dir, paste0("Heatmap_topDEG_", tt, ".png")),
      paste0("Top Differentially Expressed Genes (", tt, " post infection)")
    )
  }
  
  sig_symbols <- rownames(res$sig)
  
  if (length(sig_symbols) < 3) {
    writeLines("DEG_sig < 3, enrichment skipped.", file.path(out_enrich_dir, "TOO_FEW_DEG_SIG.txt"))
  } else {
    eg <- tryCatch(
      clusterProfiler::bitr(sig_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error = function(e) NULL
    )
    
    if (!is.null(eg) && nrow(eg) > 0) {
      write_table_dual(eg, file.path(out_enrich_dir, "DEGsig_SYMBOL_to_ENTREZ.csv"), row.names = FALSE)
      entrez <- unique(eg$ENTREZID)
      
      for (ont in c("BP","CC","MF")) {
        ego <- tryCatch(
          clusterProfiler::enrichGO(
            gene = entrez,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = ont,
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            readable = TRUE
          ),
          error = function(e) NULL
        )
        
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          go_df <- as.data.frame(ego)
          go_df <- go_df[go_df$pvalue < enrich_p_cutoff, , drop = FALSE]
          go_df <- go_df[order(go_df$pvalue, -go_df$Count), , drop = FALSE]
          
          if (nrow(go_df) > 0) {
            write_table_dual(go_df, file.path(out_enrich_dir, paste0("GO_", ont, ".csv")), row.names = FALSE)
            ego_plot <- subset_enrich_obj(ego, go_df)
            
            safe_dotplot(
              ego_plot,
              file.path(out_enrich_dir, paste0("GO_", ont, "_dotplot.png")),
              paste0("GO ", ont, " Enrichment (", tt, " post infection)"),
              15
            )
            
            safe_barplot(
              ego_plot,
              file.path(out_enrich_dir, paste0("GO_", ont, "_barplot.png")),
              paste0("GO ", ont, " Enrichment (", tt, " post infection)"),
              15
            )
            
            if (ont == "BP") {
              go_top_list[[tt]] <- head(go_df, 15)
            }
          }
        }
      }
      
      ekegg <- tryCatch(
        clusterProfiler::enrichKEGG(
          gene = as.character(entrez),
          organism = "hsa",
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1
        ),
        error = function(e) NULL
      )
      
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        kegg_df <- as.data.frame(ekegg)
        kegg_df <- kegg_df[kegg_df$pvalue < enrich_p_cutoff, , drop = FALSE]
        kegg_df <- kegg_df[order(kegg_df$pvalue, -kegg_df$Count), , drop = FALSE]
        
        if (nrow(kegg_df) > 0) {
          write_table_dual(kegg_df, file.path(out_enrich_dir, "KEGG.csv"), row.names = FALSE)
          ekegg_plot <- subset_enrich_obj(ekegg, kegg_df)
          
          safe_dotplot(
            ekegg_plot,
            file.path(out_enrich_dir, "KEGG_dotplot.png"),
            paste0("KEGG Pathway Enrichment (", tt, " post infection)"),
            15
          )
          
          safe_barplot(
            ekegg_plot,
            file.path(out_enrich_dir, "KEGG_barplot.png"),
            paste0("KEGG Pathway Enrichment (", tt, " post infection)"),
            15
          )
          
          kegg_top_list[[tt]] <- head(kegg_df, 15)
        }
      }
    }
  }
  
  tstat <- res$fit$t[, "groupInfection"]
  names(tstat) <- rownames(res$expr_sub)
  
  eg_all <- tryCatch(
    clusterProfiler::bitr(names(tstat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  
  if (!is.null(eg_all) && nrow(eg_all) > 0) {
    tdf <- data.frame(SYMBOL = names(tstat), t = as.numeric(tstat), stringsAsFactors = FALSE)
    tdf <- merge(tdf, eg_all, by = "SYMBOL")
    tdf <- tdf[order(abs(tdf$t), decreasing = TRUE), ]
    tdf <- tdf[!duplicated(tdf$ENTREZID), ]
    
    geneList <- tdf$t
    names(geneList) <- tdf$ENTREZID
    geneList <- break_ties(geneList)
    geneList <- sort(geneList, decreasing = TRUE)
    
    write_table_dual(tdf, file.path(out_gsea_dir, "GSEA_geneList_ENTREZ_t.csv"), row.names = FALSE)
    
    gsea_go <- tryCatch(
      clusterProfiler::gseGO(
        geneList = geneList,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        minGSSize = 3,
        maxGSSize = 500,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
      gsea_go_df <- as.data.frame(gsea_go)
      gsea_go_df <- gsea_go_df[order(gsea_go_df$pvalue, -abs(gsea_go_df$NES)), , drop = FALSE]
      write_table_dual(gsea_go_df, file.path(out_gsea_dir, "GSEA_GO_BP_all.csv"), row.names = FALSE)
      
      gsea_go_df_sig <- gsea_go_df[gsea_go_df$pvalue < gsea_p_cutoff, , drop = FALSE]
      write_table_dual(gsea_go_df_sig, file.path(out_gsea_dir, "GSEA_GO_BP_sig.csv"), row.names = FALSE)
      
      if (nrow(gsea_go_df_sig) > 0) {
        gsea_go_plot <- subset_enrich_obj(gsea_go, gsea_go_df_sig)
        
        safe_dotplot(
          gsea_go_plot,
          file.path(out_gsea_dir, "GSEA_GO_BP_dotplot.png"),
          paste0("GSEA GO Biological Process (", tt, " post infection)"),
          15
        )
        
        gsea_go_top_list[[tt]] <- head(gsea_go_df_sig, 20)
        
        top_ids <- head(gsea_go_df_sig$ID, 3)
        for (i in seq_along(top_ids)) {
          safe_gseaplot2(
            gsea_obj = gsea_go,
            geneSetID = top_ids[i],
            file = file.path(out_gsea_dir, paste0("GSEA_GO_BP_curve_top", i, ".png")),
            title = gsea_go_df_sig$Description[match(top_ids[i], gsea_go_df_sig$ID)]
          )
        }
      }
    }
    
    gsea_kegg <- tryCatch(
      clusterProfiler::gseKEGG(
        geneList = geneList,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        minGSSize = 3,
        maxGSSize = 500,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      gsea_kegg_df <- as.data.frame(gsea_kegg)
      gsea_kegg_df <- gsea_kegg_df[order(gsea_kegg_df$pvalue, -abs(gsea_kegg_df$NES)), , drop = FALSE]
      write_table_dual(gsea_kegg_df, file.path(out_gsea_dir, "GSEA_KEGG_all.csv"), row.names = FALSE)
      
      gsea_kegg_df_sig <- gsea_kegg_df[gsea_kegg_df$pvalue < gsea_p_cutoff, , drop = FALSE]
      write_table_dual(gsea_kegg_df_sig, file.path(out_gsea_dir, "GSEA_KEGG_sig.csv"), row.names = FALSE)
      
      if (nrow(gsea_kegg_df_sig) > 0) {
        gsea_kegg_plot <- subset_enrich_obj(gsea_kegg, gsea_kegg_df_sig)
        
        safe_dotplot(
          gsea_kegg_plot,
          file.path(out_gsea_dir, "GSEA_KEGG_dotplot.png"),
          paste0("GSEA KEGG Pathways (", tt, " post infection)"),
          15
        )
        
        gsea_kegg_top_list[[tt]] <- head(gsea_kegg_df_sig, 20)
        
        top_ids <- head(gsea_kegg_df_sig$ID, 3)
        for (i in seq_along(top_ids)) {
          safe_gseaplot2(
            gsea_obj = gsea_kegg,
            geneSetID = top_ids[i],
            file = file.path(out_gsea_dir, paste0("GSEA_KEGG_curve_top", i, ".png")),
            title = gsea_kegg_df_sig$Description[match(top_ids[i], gsea_kegg_df_sig$ID)]
          )
        }
      }
    }
  }
}

deg_union <- sort(unique(deg_union))
write_table_dual(
  data.frame(gene_symbol = deg_union),
  file.path(out_root, "Figure_3", "DEG_union.csv"),
  row.names = FALSE
)

cat("\nDEG_union size =", length(deg_union), "\n\n")

############################
# 9) Figure 1
############################
cat("=== [Figure 1] ===\n")

volcano_plots <- list()
for (tt in names(time_list)) {
  tab <- deg_all_list[[tt]]
  df <- tab
  df$gene <- rownames(df)
  df$group <- "NS"
  df$group[df$P.Value < deg_p_cutoff & df$logFC >  deg_fc_cutoff] <- "Up"
  df$group[df$P.Value < deg_p_cutoff & df$logFC < -deg_fc_cutoff] <- "Down"
  df$negLogP <- -log10(df$P.Value)
  
  volcano_plots[[tt]] <- ggplot(df, aes(x = logFC, y = negLogP, color = group)) +
    geom_point(size = 1.2, alpha = 0.75) +
    geom_vline(xintercept = c(-deg_fc_cutoff, deg_fc_cutoff), linetype = 2, linewidth = 0.4) +
    geom_hline(yintercept = -log10(deg_p_cutoff), linetype = 2, linewidth = 0.4) +
    scale_color_manual(values = deg_cols) +
    labs(
      title = paste0(tt, " post infection"),
      x = "log2FC",
      y = expression(-log[10](italic(P)))
    ) +
    theme_sci(11) +
    theme(legend.position = "none")
}

deg_summary_df <- dplyr::bind_rows(deg_summary_list)

deg_plot_df <- deg_summary_df %>%
  dplyr::select(Time, Up, Down) %>%
  tidyr::pivot_longer(
    cols = c("Up", "Down"),
    names_to = "Direction",
    values_to = "Count"
  )
deg_plot_df$Direction <- factor(deg_plot_df$Direction, levels = c("Up","Down"))

p_deg_bar <- ggplot(deg_plot_df, aes(x = Time, y = Count, fill = Direction)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  labs(title = "Differentially expressed genes across time points", x = NULL, y = "Gene count") +
  theme_sci(12)

p_pca_small_ellipse   <- p_pca + theme(legend.position = "bottom")
p_pca_small_noellipse <- p_pca_noellipse + theme(legend.position = "bottom")
p_heatmap_img <- read_png_as_plot(file.path(out_root, "01_PCA", "Sample_correlation_heatmap.png"))

fig1_volcano <- (volcano_plots[["4h"]] | volcano_plots[["8h"]]) /
  (volcano_plots[["12h"]] | volcano_plots[["16h"]])

fig1_top <- p_pca_small_ellipse | p_pca_small_noellipse | p_heatmap_img
fig1_bottom <- fig1_volcano | p_deg_bar

fig1 <- fig1_top / fig1_bottom +
  plot_annotation(title = "Figure 1. Global transcriptomic landscape and differential expression overview")

safe_save_plot(fig1, file.path(out_root, "Figure_1", "Figure_1.png"), width = 22, height = 14, dpi = 320)
safe_save_plot(fig1, file.path(out_root, "Figure_1", "Figure_1.pdf"), width = 22, height = 14, dpi = 320)

############################
# 10) Figure 2
############################
cat("=== [Figure 2] ===\n")

go_fig_df        <- pick_latest_nonempty(go_top_list)
kegg_fig_df      <- pick_latest_nonempty(kegg_top_list)
gsea_go_fig_df   <- pick_latest_nonempty(gsea_go_top_list)
gsea_kegg_fig_df <- pick_latest_nonempty(gsea_kegg_top_list)

if (!is.null(go_fig_df)) {
  p_go_bubble <- make_enrichment_bubble(
    go_fig_df,
    file.path(out_root, "Figure_2", "Figure_2_GO_bubble.png"),
    "GO biological process enrichment",
    15
  )
} else {
  p_go_bubble <- ggplot() + theme_void() + labs(title = "No GO enrichment available")
}

if (!is.null(kegg_fig_df)) {
  p_kegg_bubble <- make_enrichment_bubble(
    kegg_fig_df,
    file.path(out_root, "Figure_2", "Figure_2_KEGG_bubble.png"),
    "KEGG pathway enrichment",
    15
  )
} else {
  p_kegg_bubble <- ggplot() + theme_void() + labs(title = "No KEGG enrichment available")
}

if (!is.null(gsea_go_fig_df)) {
  p_gsea_go <- make_gsea_nes_plot(
    gsea_go_fig_df,
    file.path(out_root, "Figure_2", "Figure_2_GSEA_GO_NES.png"),
    "GSEA GO biological process: NES direction",
    10
  )
} else {
  p_gsea_go <- ggplot() + theme_void() + labs(title = "No GSEA GO result available")
}

if (!is.null(gsea_kegg_fig_df)) {
  p_gsea_kegg <- make_gsea_nes_plot(
    gsea_kegg_fig_df,
    file.path(out_root, "Figure_2", "Figure_2_GSEA_KEGG_NES.png"),
    "GSEA KEGG pathways: NES direction",
    10
  )
} else {
  p_gsea_kegg <- ggplot() + theme_void() + labs(title = "No GSEA KEGG result available")
}

fig2 <- (p_go_bubble | p_kegg_bubble) / (p_gsea_go | p_gsea_kegg) +
  plot_annotation(title = "Figure 2. Functional enrichment and pathway activation landscape")

safe_save_plot(fig2, file.path(out_root, "Figure_2", "Figure_2.png"), width = 16, height = 12, dpi = 320)
safe_save_plot(fig2, file.path(out_root, "Figure_2", "Figure_2.pdf"), width = 16, height = 12, dpi = 320)

############################
# 11) Figure 3
############################
cat("=== [Figure 3] ===\n")

if (length(deg_union) >= 2) {
  mean_expr <- function(sids) rowMeans(expr_gene[, sids, drop = FALSE], na.rm = TRUE)
  expr_union_mean <- cbind(
    `4h`  = mean_expr(inf_4h),
    `8h`  = mean_expr(inf_8h),
    `12h` = mean_expr(inf_12h),
    `16h` = mean_expr(inf_16h)
  )
  expr_union_mean <- expr_union_mean[deg_union, , drop = FALSE]
  
  top_dynamic_genes <- names(sort(apply(expr_union_mean, 1, var, na.rm = TRUE), decreasing = TRUE))[1:min(20, nrow(expr_union_mean))]
  trend_df <- data.frame(
    gene = rep(top_dynamic_genes, times = ncol(expr_union_mean)),
    time = rep(colnames(expr_union_mean), each = length(top_dynamic_genes)),
    value = as.vector(t(scale(t(expr_union_mean[top_dynamic_genes, , drop = FALSE])))),
    stringsAsFactors = FALSE
  )
  trend_df$time <- factor(trend_df$time, levels = c("4h","8h","12h","16h"))
  
  p_trend <- ggplot(trend_df, aes(x = time, y = value, group = gene)) +
    geom_line(alpha = 0.35, color = "grey40") +
    stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.4, color = "#D73027") +
    stat_summary(aes(group = 1), fun = mean, geom = "point", size = 3, color = "#D73027") +
    labs(
      title = "Temporal trend of highly dynamic DEGs in infected samples",
      x = NULL,
      y = "Scaled expression"
    ) +
    theme_sci(12)
  
  safe_save_plot(p_trend, file.path(out_root, "Figure_3", "Figure_3_DEG_trend_plot.png"), width = 8, height = 6, dpi = 320)
  safe_save_plot(p_trend, file.path(out_root, "Figure_3", "Figure_3_DEG_trend_plot.pdf"), width = 8, height = 6, dpi = 320)
} else {
  p_trend <- ggplot() + theme_void() + labs(title = "No DEG trend available")
}

fig3 <- p_trend +
  plot_annotation(title = "Figure 3. Temporal trend of highly dynamic DEGs")

safe_save_plot(fig3, file.path(out_root, "Figure_3", "Figure_3.png"), width = 8, height = 6.5, dpi = 320)
safe_save_plot(fig3, file.path(out_root, "Figure_3", "Figure_3.pdf"), width = 8, height = 6.5, dpi = 320)

cat("\n============================================================\n")
cat("PIPELINE FINISHED.\n")
cat("Output:", out_root, "\n")
cat("Generated:\n")
cat("  PCA with ellipse + no ellipse\n")
cat("  GSEA curve plots (top 3 per timepoint when available)\n")
cat("  Figure 1: PCA + sample heatmap + volcano + DEG barplot\n")
cat("  Figure 2: enrichment bubble + NES direction\n")
cat("  Figure 3: DEG temporal trend only\n")
cat("  All CSV tables also exported as same-name TXT files\n")
cat("============================================================\n")