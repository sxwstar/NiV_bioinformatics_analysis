# ============================================================
# Drug reversal analysis using ONLY matched_drugs.csv
# Publication-style heatmap version




# ============================================================

setwd("/users/hzhang1/project/Cmap")

suppressPackageStartupMessages({
  library(cmapR)
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ==================== 参数 ====================
gctx_file      <- "level5_beta_trt_cp_n720216x12328.gctx"
matched_file   <- "matched_drugs.csv"
up_file        <- "up_genes.txt"
down_file      <- "down_genes.txt"
geneinfo_file  <- "geneinfo_beta.txt"

output_dir     <- "drug_reversal_matched_pubheatmap_virus_genes"

# 热图颜色范围截断分位数
# 如果你觉得白色仍然偏多，可改成 0.90
heatmap_quantile <- 0.95

# 是否导出每个signature的明细CSV
export_per_signature_csv <- TRUE

# 全基因热图是否显示列名
# 基因非常多时建议 FALSE
show_colnames_all <- FALSE

# ==================== 创建输出目录 ====================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "per_signature_csv"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "matrices"), showWarnings = FALSE)

# ==================== 工具函数 ====================
safe_filename <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  x <- gsub("[/:*?\"<>|\\\\]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

make_unique_labels <- function(x) {
  make.unique(as.character(x), sep = "_")
}

row_dist_fun <- function(m) {
  m <- as.matrix(m)
  as.dist(1 - cor(t(m), method = "spearman", use = "pairwise.complete.obs"))
}

col_dist_fun <- function(m) {
  m <- as.matrix(m)
  as.dist(1 - cor(m, method = "spearman", use = "pairwise.complete.obs"))
}

get_clip_value <- function(mat, q = 0.95) {
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (length(vals) == 0) return(2)
  x <- as.numeric(quantile(abs(vals), probs = q, na.rm = TRUE))
  if (!is.finite(x) || is.na(x) || x <= 0) {
    x <- max(abs(vals), na.rm = TRUE)
  }
  if (!is.finite(x) || is.na(x) || x <= 0) {
    x <- 2
  }
  x
}

clip_matrix <- function(mat, clip_val) {
  m <- mat
  m[m >  clip_val] <-  clip_val
  m[m < -clip_val] <- -clip_val
  m
}

# ==================== 顶刊风格热图函数 ====================
draw_pub_heatmap <- function(
    mat,
    row_anno_df,
    col_gene_set,
    out_prefix,
    title_text,
    show_colnames = FALSE,
    pdf_width = 12,
    pdf_height = 8,
    png_width = 12,
    png_height = 8,
    fontsize_row = 10,
    fontsize_col = 8
) {
  mat <- as.matrix(mat)
  
  keep_cols <- colSums(!is.na(mat)) > 0
  mat <- mat[, keep_cols, drop = FALSE]
  col_gene_set <- col_gene_set[colnames(mat), , drop = FALSE]
  
  clip_val <- get_clip_value(mat, q = heatmap_quantile)
  mat_plot <- clip_matrix(mat, clip_val)
  
  col_fun_main <- colorRamp2(
    c(-clip_val, 0, clip_val),
    c("#3B4CC0", "white", "#B40426")
  )
  
  gene_set_cols <- c(
    "up" = "#F8766D",
    "down" = "#00BFC4",
    "both" = "#7E62A3",
    "unknown" = "#BDBDBD"
  )
  
  score_rng <- range(row_anno_df$reversal_score, na.rm = TRUE)
  prop_rng  <- range(row_anno_df$prop_reversed, na.rm = TRUE)
  
  if (!all(is.finite(score_rng))) score_rng <- c(0, 1)
  if (!all(is.finite(prop_rng)))  prop_rng  <- c(0, 1)
  
  if (diff(score_rng) == 0) score_rng <- score_rng + c(-0.001, 0.001)
  if (diff(prop_rng)  == 0) prop_rng  <- prop_rng  + c(-0.001, 0.001)
  
  col_fun_score <- colorRamp2(
    c(score_rng[1], mean(score_rng), score_rng[2]),
    c("#EAF2F8", "#7FB3D5", "#1F618D")
  )
  col_fun_prop <- colorRamp2(
    c(prop_rng[1], mean(prop_rng), prop_rng[2]),
    c("#F4ECF7", "#BB8FCE", "#6C3483")
  )
  
  ha_top <- HeatmapAnnotation(
    gene_set = col_gene_set$gene_set,
    col = list(gene_set = gene_set_cols),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    gp = gpar(col = NA),
    simple_anno_size = unit(4, "mm")
  )
  
  ha_left <- rowAnnotation(
    reversal_score = row_anno_df$reversal_score,
    prop_reversed = row_anno_df$prop_reversed,
    col = list(
      reversal_score = col_fun_score,
      prop_reversed = col_fun_prop
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      reversal_score = list(title = "Reversal score"),
      prop_reversed = list(title = "Prop. reversed")
    ),
    gap = unit(1.5, "mm"),
    width = unit(c(5, 5), "mm")
  )
  
  column_split <- factor(
    col_gene_set$gene_set,
    levels = c("up", "down", "both", "unknown")
  )
  column_split <- droplevels(column_split)
  
  ht <- Heatmap(
    mat_plot,
    name = "Reversal\nvalue",
    col = col_fun_main,
    na_col = "#F2F2F2",
    top_annotation = ha_top,
    left_annotation = ha_left,
    
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = row_dist_fun,
    clustering_distance_columns = col_dist_fun,
    clustering_method_rows = "average",
    clustering_method_columns = "average",
    
    row_names_side = "right",
    row_names_gp = gpar(fontsize = fontsize_row),
    row_names_max_width = unit(8, "cm"),
    
    show_column_names = show_colnames,
    column_names_gp = gpar(fontsize = fontsize_col),
    column_names_rot = 45,
    
    row_title = NULL,
    column_title = title_text,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    border = FALSE,
    rect_gp = gpar(col = NA),
    
    column_split = column_split,
    column_gap = unit(1.5, "mm"),
    
    heatmap_legend_param = list(
      title = "Reversal value",
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(45, "mm"),
      at = c(-clip_val, 0, clip_val),
      labels = c(
        sprintf("%.2f", -clip_val),
        "0",
        sprintf("%.2f", clip_val)
      )
    ),
    
    use_raster = TRUE,
    raster_quality = 3
  )
  
  pdf(paste0(out_prefix, ".pdf"), width = pdf_width, height = pdf_height)
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(6, 6, 6, 6), "mm")
  )
  dev.off()
  
  png(
    filename = paste0(out_prefix, ".png"),
    width = png_width, height = png_height,
    units = "in", res = 300
  )
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = unit(c(6, 6, 6, 6), "mm")
  )
  dev.off()
  
  invisible(list(
    heatmap = ht,
    clip_value = clip_val,
    matrix_used = mat_plot
  ))
}

# ==================== 读取 gene info ====================
gene_info <- fread(geneinfo_file)

required_geneinfo_cols <- c("gene_id", "gene_symbol")
missing_geneinfo_cols <- setdiff(required_geneinfo_cols, colnames(gene_info))
if (length(missing_geneinfo_cols) > 0) {
  stop("geneinfo_beta.txt 缺少以下列: ",
       paste(missing_geneinfo_cols, collapse = ", "))
}

gene_info <- gene_info %>%
  mutate(
    gene_id = as.character(gene_id),
    gene_symbol = as.character(gene_symbol),
    gene_symbol_upper = str_to_upper(gene_symbol)
  )

symbol_to_entrez <- setNames(gene_info$gene_id, gene_info$gene_symbol_upper)
entrez_to_symbol <- setNames(gene_info$gene_symbol, gene_info$gene_id)

# ==================== 读取上下调基因 ====================
up_genes <- read_lines(up_file) %>%
  str_trim() %>%
  discard(~ .x == "") %>%
  str_to_upper() %>%
  unique()

down_genes <- read_lines(down_file) %>%
  str_trim() %>%
  discard(~ .x == "") %>%
  str_to_upper() %>%
  unique()

up_entrez <- unname(symbol_to_entrez[up_genes]) %>%
  na.omit() %>%
  as.character() %>%
  unique()

down_entrez <- unname(symbol_to_entrez[down_genes]) %>%
  na.omit() %>%
  as.character() %>%
  unique()

target_entrez <- unique(c(up_entrez, down_entrez))

overlap_entrez <- intersect(up_entrez, down_entrez)
if (length(overlap_entrez) > 0) {
  warning("up/down 基因有重叠，共 ", length(overlap_entrez),
          " 个。这些基因将标记为 both，且不参与 reversal_score 计算。")
}

gene_set_map <- data.frame(
  gene_id = target_entrez,
  stringsAsFactors = FALSE
) %>%
  mutate(
    gene_id = as.character(gene_id),
    in_up = gene_id %in% up_entrez,
    in_down = gene_id %in% down_entrez,
    gene_set = case_when(
      in_up & in_down ~ "both",
      in_up ~ "up",
      in_down ~ "down",
      TRUE ~ "unknown"
    ),
    input_symbol = unname(entrez_to_symbol[gene_id])
  ) %>%
  select(gene_id, gene_set, input_symbol)

cat("目标基因总数：", length(target_entrez), "\n")
cat("up基因数：", length(up_entrez), "\n")
cat("down基因数：", length(down_entrez), "\n")

# ==================== 读取 matched_drugs ====================
matched <- fread(matched_file)

required_cols <- c("id", "pert_iname", "cell_iname", "pert_idose", "pert_itime")
missing_cols <- setdiff(required_cols, colnames(matched))
if (length(missing_cols) > 0) {
  stop("matched_drugs.csv 缺少以下列: ", paste(missing_cols, collapse = ", "))
}

matched <- matched %>%
  mutate(
    id = as.character(id),
    pert_iname = as.character(pert_iname),
    cell_iname = as.character(cell_iname),
    pert_idose = as.character(pert_idose),
    pert_itime = as.character(pert_itime)
  )

sig_ids <- unique(matched$id)
if (length(sig_ids) == 0) {
  stop("matched_drugs.csv 中没有可用的 signature id。")
}

cat("matched_drugs 中 signature 数：", length(sig_ids), "\n")

# ==================== 核心函数：计算单个 signature ====================
compute_one_signature <- function(sig_id, meta_row, gctx_file, gene_set_map, entrez_to_symbol) {
  gct <- parse_gctx(gctx_file, cid = sig_id)
  
  if (ncol(gct@mat) != 1) {
    stop("sig_id=", sig_id, " 读取后不是单列矩阵。")
  }
  
  zscores <- gct@mat[, 1]
  gene_ids_in_gctx <- rownames(gct@mat)
  
  if (is.null(gene_ids_in_gctx)) {
    stop("sig_id=", sig_id, " 的表达矩阵没有行名。")
  }
  
  zscores <- as.numeric(zscores)
  names(zscores) <- as.character(gene_ids_in_gctx)
  
  ranked <- sort(zscores, decreasing = TRUE)
  
  ranked_df <- data.frame(
    gene_id = names(ranked),
    zscore = as.numeric(ranked),
    rank = seq_along(ranked),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene_id = as.character(gene_id),
      rank_percent = rank / n(),
      regulation_in_drug = case_when(
        zscore > 0 ~ "up",
        zscore < 0 ~ "down",
        TRUE ~ "neutral"
      ),
      abs_zscore = abs(zscore)
    )
  
  drug_gene_effect <- ranked_df %>%
    inner_join(gene_set_map, by = "gene_id") %>%
    mutate(
      gene_symbol = unname(entrez_to_symbol[gene_id]),
      gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", input_symbol, gene_symbol),
      
      reversal_value = case_when(
        gene_set == "up"   ~ -zscore,
        gene_set == "down" ~  zscore,
        gene_set == "both" ~ NA_real_,
        TRUE ~ NA_real_
      ),
      
      reversed = case_when(
        is.na(reversal_value) ~ NA,
        reversal_value > 0 ~ TRUE,
        reversal_value <= 0 ~ FALSE
      ),
      
      reversal_strength = case_when(
        is.na(reversal_value) ~ NA_real_,
        reversal_value > 0 ~ reversal_value,
        reversal_value <= 0 ~ 0
      ),
      
      same_direction = case_when(
        is.na(reversal_value) ~ NA,
        reversal_value < 0 ~ TRUE,
        reversal_value >= 0 ~ FALSE
      ),
      
      drug_name = meta_row$pert_iname,
      sig_id = sig_id,
      cell_id = meta_row$cell_iname,
      pert_dose = meta_row$pert_idose,
      pert_time = meta_row$pert_itime
    ) %>%
    select(
      drug_name,
      sig_id,
      cell_id,
      pert_dose,
      pert_time,
      gene_symbol,
      gene_id,
      gene_set,
      zscore,
      abs_zscore,
      regulation_in_drug,
      reversal_value,
      reversed,
      reversal_strength,
      same_direction,
      rank,
      rank_percent
    ) %>%
    arrange(desc(reversal_strength), desc(abs_zscore), rank)
  
  valid_df <- drug_gene_effect %>%
    filter(!is.na(reversal_value))
  
  n_total_genes <- nrow(valid_df)
  n_reversed <- sum(valid_df$reversed %in% TRUE, na.rm = TRUE)
  
  summary_row <- data.frame(
    drug_name = meta_row$pert_iname,
    sig_id = sig_id,
    cell_id = meta_row$cell_iname,
    pert_dose = meta_row$pert_idose,
    pert_time = meta_row$pert_itime,
    
    n_target_genes_found = nrow(drug_gene_effect),
    n_up_genes_found = sum(drug_gene_effect$gene_set %in% c("up", "both"), na.rm = TRUE),
    n_down_genes_found = sum(drug_gene_effect$gene_set %in% c("down", "both"), na.rm = TRUE),
    
    mean_abs_zscore = ifelse(nrow(drug_gene_effect) > 0,
                             mean(drug_gene_effect$abs_zscore, na.rm = TRUE),
                             NA_real_),
    mean_zscore = ifelse(nrow(drug_gene_effect) > 0,
                         mean(drug_gene_effect$zscore, na.rm = TRUE),
                         NA_real_),
    
    n_total_genes_for_reversal = n_total_genes,
    n_reversed = n_reversed,
    prop_reversed = ifelse(n_total_genes > 0, n_reversed / n_total_genes, NA_real_),
    mean_reversal_strength = ifelse(n_total_genes > 0,
                                    mean(valid_df$reversal_strength, na.rm = TRUE),
                                    NA_real_),
    reversal_score = ifelse(n_total_genes > 0,
                            sum(valid_df$reversal_strength, na.rm = TRUE) / n_total_genes,
                            NA_real_),
    
    stringsAsFactors = FALSE
  )
  
  list(
    gene_table = drug_gene_effect,
    summary_row = summary_row
  )
}

# ==================== 主循环 ====================
summary_list <- vector("list", length(sig_ids))
all_gene_effect_list <- vector("list", length(sig_ids))

pb <- txtProgressBar(min = 0, max = length(sig_ids), style = 3)

for (i in seq_along(sig_ids)) {
  sig_id <- sig_ids[i]
  
  tryCatch({
    meta_row <- matched %>% filter(id == sig_id) %>% slice(1)
    
    res <- compute_one_signature(
      sig_id = sig_id,
      meta_row = meta_row,
      gctx_file = gctx_file,
      gene_set_map = gene_set_map,
      entrez_to_symbol = entrez_to_symbol
    )
    
    summary_list[[i]] <- res$summary_row
    all_gene_effect_list[[i]] <- res$gene_table
    
    if (export_per_signature_csv) {
      out_name <- paste0(
        sprintf("%03d", i), "_",
        safe_filename(meta_row$pert_iname), "_",
        safe_filename(meta_row$cell_iname), "_",
        safe_filename(meta_row$pert_idose), "_",
        safe_filename(meta_row$pert_itime), "_",
        safe_filename(sig_id),
        ".csv"
      )
      
      write_csv(
        res$gene_table,
        file.path(output_dir, "per_signature_csv", out_name)
      )
    }
    
  }, error = function(e) {
    message("Error in sig_id ", sig_id, ": ", e$message)
    
    meta_row <- matched %>% filter(id == sig_id) %>% slice(1)
    
    summary_list[[i]] <<- data.frame(
      drug_name = ifelse(nrow(meta_row) > 0, meta_row$pert_iname, NA),
      sig_id = sig_id,
      cell_id = ifelse(nrow(meta_row) > 0, meta_row$cell_iname, NA),
      pert_dose = ifelse(nrow(meta_row) > 0, meta_row$pert_idose, NA),
      pert_time = ifelse(nrow(meta_row) > 0, meta_row$pert_itime, NA),
      
      n_target_genes_found = NA,
      n_up_genes_found = NA,
      n_down_genes_found = NA,
      mean_abs_zscore = NA,
      mean_zscore = NA,
      
      n_total_genes_for_reversal = NA,
      n_reversed = NA,
      prop_reversed = NA,
      mean_reversal_strength = NA,
      reversal_score = NA,
      stringsAsFactors = FALSE
    )
    
    all_gene_effect_list[[i]] <<- NULL
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

# ==================== 汇总输出 ====================
drug_gene_effect_summary <- bind_rows(summary_list) %>%
  arrange(desc(reversal_score), desc(prop_reversed), desc(mean_reversal_strength))

all_gene_effect_df <- bind_rows(all_gene_effect_list)

write_csv(drug_gene_effect_summary, file.path(output_dir, "drug_gene_effect_summary.csv"))
write_csv(all_gene_effect_df, file.path(output_dir, "all_drug_gene_effect_long_table.csv"))

# ==================== 构建热图矩阵 ====================
heatmap_df <- all_gene_effect_df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(drug_name, gene_symbol) %>%
  summarise(
    reversal_value = mean(reversal_value, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_mat <- heatmap_df %>%
  pivot_wider(names_from = gene_symbol, values_from = reversal_value) %>%
  as.data.frame()

rownames(heatmap_mat) <- make_unique_labels(heatmap_mat$drug_name)
drug_labels <- heatmap_mat$drug_name
heatmap_mat$drug_name <- NULL
heatmap_mat <- as.matrix(heatmap_mat)

ordered_drugs <- drug_gene_effect_summary$drug_name
row_map <- data.frame(
  row_id = rownames(heatmap_mat),
  drug_name = drug_labels,
  stringsAsFactors = FALSE
)

ordered_row_ids <- row_map %>%
  mutate(order_idx = match(drug_name, ordered_drugs)) %>%
  arrange(order_idx) %>%
  pull(row_id)

ordered_row_ids <- ordered_row_ids[ordered_row_ids %in% rownames(heatmap_mat)]
heatmap_mat <- heatmap_mat[ordered_row_ids, , drop = FALSE]

annotation_row_df <- row_map %>%
  mutate(order_idx = match(drug_name, ordered_drugs)) %>%
  arrange(order_idx) %>%
  left_join(
    drug_gene_effect_summary %>%
      select(drug_name, reversal_score, prop_reversed),
    by = "drug_name"
  )

rownames(annotation_row_df) <- annotation_row_df$row_id
annotation_row_df <- annotation_row_df[rownames(heatmap_mat), , drop = FALSE]

row_labels_final <- annotation_row_df$drug_name
rownames(annotation_row_df) <- row_labels_final
rownames(heatmap_mat) <- row_labels_final

annotation_row_df <- annotation_row_df %>%
  select(reversal_score, prop_reversed)

gene_anno <- all_gene_effect_df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, gene_set)

up_cols <- gene_anno %>%
  filter(gene_set == "up") %>%
  pull(gene_symbol) %>%
  unique() %>%
  sort()

down_cols <- gene_anno %>%
  filter(gene_set == "down") %>%
  pull(gene_symbol) %>%
  unique() %>%
  sort()

both_cols <- gene_anno %>%
  filter(gene_set == "both") %>%
  pull(gene_symbol) %>%
  unique() %>%
  sort()

ordered_cols <- c(up_cols, down_cols, both_cols)
ordered_cols <- intersect(ordered_cols, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, ordered_cols, drop = FALSE]

annotation_col_df <- data.frame(
  gene_set = case_when(
    colnames(heatmap_mat) %in% up_cols ~ "up",
    colnames(heatmap_mat) %in% down_cols ~ "down",
    colnames(heatmap_mat) %in% both_cols ~ "both",
    TRUE ~ "unknown"
  ),
  stringsAsFactors = FALSE
)
rownames(annotation_col_df) <- colnames(heatmap_mat)

# 输出矩阵和注释
write.csv(
  heatmap_mat,
  file.path(output_dir, "matrices", "heatmap_reversal_value_matrix.csv"),
  row.names = TRUE
)
write.csv(
  annotation_row_df,
  file.path(output_dir, "matrices", "heatmap_drug_annotation.csv"),
  row.names = TRUE
)
write.csv(
  annotation_col_df,
  file.path(output_dir, "matrices", "heatmap_gene_annotation.csv"),
  row.names = TRUE
)

# ==================== 全基因热图 ====================
all_heatmap_res <- draw_pub_heatmap(
  mat = heatmap_mat,
  row_anno_df = annotation_row_df,
  col_gene_set = annotation_col_df,
  out_prefix = file.path(output_dir, "plots", "heatmap_all_genes_publication_style"),
  title_text = paste0("Drug-Gene Reversal Heatmap (All genes, q=", heatmap_quantile, ")"),
  show_colnames = show_colnames_all,
  pdf_width = 13.5,
  pdf_height = 9,
  png_width = 13.5,
  png_height = 9,
  fontsize_row = 10,
  fontsize_col = 8
)

write.csv(
  all_heatmap_res$matrix_used,
  file.path(output_dir, "matrices", "heatmap_all_genes_matrix_clipped_for_plot.csv"),
  row.names = TRUE
)

# ==================== reversal_score条形图 ====================
plot_df1 <- drug_gene_effect_summary %>%
  filter(!is.na(reversal_score)) %>%
  arrange(reversal_score) %>%
  mutate(drug_name = factor(drug_name, levels = drug_name))

p1 <- ggplot(plot_df1, aes(x = drug_name, y = reversal_score)) +
  geom_col(fill = "#B40426", width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Drug reversal score",
    x = "Drug",
    y = "Reversal score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "plots", "barplot_reversal_score.png"),
  plot = p1, width = 9, height = 7, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "plots", "barplot_reversal_score.pdf"),
  plot = p1, width = 9, height = 7
)

# ==================== 散点图 ====================
plot_df2 <- drug_gene_effect_summary %>%
  filter(!is.na(prop_reversed), !is.na(reversal_score))

p2 <- ggplot(plot_df2, aes(x = prop_reversed, y = reversal_score, label = drug_name)) +
  geom_point(size = 3.2, color = "#B40426") +
  geom_text_repel(size = 3.5, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Drug reversal overview",
    x = "Proportion reversed",
    y = "Reversal score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "plots", "scatter_reversal_score_vs_prop.png"),
  plot = p2, width = 8.5, height = 7, dpi = 300
)
ggsave(
  filename = file.path(output_dir, "plots", "scatter_reversal_score_vs_prop.pdf"),
  plot = p2, width = 8.5, height = 7
)

# ==================== 输出说明 ====================
readme_txt <- c(
  "结果说明",
  "",
  "一、分析逻辑",
  "1. 仅使用 matched_drugs.csv 中提供的 signature",
  "2. 每个signature单独计算 reversal 指标",
  "3. 使用 reversal_value 绘制热图，而不再直接使用原始 zscore",
  "4. 不再绘制 top50 热图，只保留一个全基因热图",
  "",
  "二、关键字段定义",
  "reversal_value:",
  "  up gene   -> -zscore",
  "  down gene ->  zscore",
  "  >0 表示逆转，<0 表示同向",
  "",
  "reversed:",
  "  reversal_value > 0 为 TRUE，否则 FALSE",
  "",
  "reversal_strength:",
  "  若 reversal_value > 0，则等于 reversal_value",
  "  否则为 0",
  "",
  "reversal_score:",
  "  sum(reversal_strength) / n_total_genes_for_reversal",
  "",
  "三、热图风格说明",
  "1. 使用 ComplexHeatmap 绘制顶刊风格热图",
  "2. 左侧仅保留 reversal_score 与 prop_reversed 两个行注释",
  "3. 主热图颜色为蓝-白-红：负值=同向，正值=逆转",
  paste0("4. 颜色范围使用 abs(reversal_value) 的 ", heatmap_quantile * 100, "% 分位数截断"),
  "5. 全基因热图默认不显示列名，以避免拥挤",
  "",
  "四、主要输出文件",
  "1. drug_gene_effect_summary.csv",
  "2. all_drug_gene_effect_long_table.csv",
  "3. per_signature_csv/",
  "4. plots/heatmap_all_genes_publication_style.pdf/png",
  "5. matrices/"
)

writeLines(readme_txt, con = file.path(output_dir, "README_results.txt"))

cat("\n========================================\n")
cat("完成！结果输出目录：", output_dir, "\n")
cat("summary：drug_gene_effect_summary.csv\n")
cat("全基因热图：heatmap_all_genes_publication_style.pdf/png\n")
cat("========================================\n")





############################################################
## drug_reversal_virus_genes.R
##
## Purpose:
##   Drug reversal analysis for virus-tagged host factors.
##
## Working directory:
##   /users/hzhang1/project/Cmap
##
## Required input files:
##   level5_beta_trt_cp_n720216x12328.gctx
##   geneinfo_beta.txt
##   drug_list.txt
##   virus_matched_drugs.csv
##   virus_up_genes.txt
##   virus_down_genes.txt
##
## Output directory:
##   /users/hzhang1/project/Cmap/drug_reversal_virus_genes
##
## Biological convention:
##   virus_up_genes.txt   -> HDFs
##   virus_down_genes.txt -> HRFs
##
## Reversal definition:
##   HDFs / virus_up_genes.txt:
##     reversal_value = -zscore
##
##   HRFs / virus_down_genes.txt:
##     reversal_value =  zscore
##
##   reversal_value > 0:
##     drug-induced reversal of virus-associated host-factor dysregulation
##
##   reversal_value < 0:
##     same-direction perturbation
##
## Main SCI figure:
##   Virus-tagged host-factor reversal profiles of candidate drugs
##
## Final figure style:
##   Same plotting style as residual_gene_reversal_integrated_heatmap_FINAL_SCI_v12.R
##
## Figure components:
##   1. left narrow Reversal score strip
##   2. left narrow Reversal proportion strip
##   3. gene-level reversal-value heatmap
##   4. top HDFs / HRFs annotation
##   5. fixed right legend panel
##
##
############################################################


############################################################
## 0. Clean environment
############################################################
rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
set.seed(123)

start_time <- Sys.time()


############################################################
## 1. Load packages
############################################################
required_pkgs <- c(
  "cmapR",
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2",
  "ComplexHeatmap",
  "circlize",
  "grid"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "缺少 R 包：",
    paste(missing_pkgs, collapse = ", "),
    "\n请先在 r45 环境中安装这些包。"
  )
}

suppressPackageStartupMessages({
  library(cmapR)
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

ht_opt$message <- FALSE


############################################################
## 2. Working directory and paths
##    注意：输入与输出地址保持病毒基因分析原路径
############################################################
setwd("/users/hzhang1/project/Cmap")

gctx_file      <- "level5_beta_trt_cp_n720216x12328.gctx"
geneinfo_file  <- "geneinfo_beta.txt"
drug_list_file <- "drug_list.txt"

matched_file   <- "virus_matched_drugs.csv"
up_file        <- "virus_up_genes.txt"
down_file      <- "virus_down_genes.txt"

output_dir   <- "/users/hzhang1/project/Cmap/drug_reversal_virus_genes"
plots_dir    <- file.path(output_dir, "plots")
matrices_dir <- file.path(output_dir, "matrices")
per_sig_dir  <- file.path(output_dir, "per_signature_csv")
checks_dir   <- file.path(output_dir, "input_checks")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrices_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_sig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checks_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  gctx_file,
  geneinfo_file,
  drug_list_file,
  matched_file,
  up_file,
  down_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "以下输入文件不存在，请检查路径：\n",
    paste(missing_files, collapse = "\n")
  )
}

## 如果文件系统仍是只读，这里直接给出清晰报错
if (file.access(output_dir, 2) != 0) {
  stop(
    "输出目录当前不可写：\n",
    output_dir,
    "\n请先恢复该路径写权限后再运行。"
  )
}


############################################################
## 3. Global settings
##    与 residual v12 主图保持同一画图风格
############################################################
global_font_family <- "sans"

main_title <- "Virus-tagged host-factor reversal profiles of candidate drugs"

heatmap_clip_value <- 1.5

pdf_width  <- 15.8
pdf_height <- 8.2
plot_dpi   <- 600

export_per_signature_csv <- TRUE


############################################################
## 4. Color system
##    与 residual v12 保持同一配色逻辑
############################################################
## Core palette
col_text <- "#222222"
col_bg0  <- "#F2ECE4"   # not pure white, softer center

## Top annotation colors
## 这里用 residual up/down 的视觉风格映射 HDFs / HRFs
col_hdf <- "#FF4D2D"    # same bright red-orange style
col_hrf <- "#2C7FB8"    # same navy blue style
col_both <- "#7A7A7A"
col_unknown <- "#BDBDBD"

host_factor_cols <- c(
  "HDFs" = col_hdf,
  "HRFs" = col_hrf,
  "Both" = col_both,
  "Unknown" = col_unknown
)

## Heatmap main colors: wine red <-> soft neutral <-> navy blue
col_wine_dark  <- "#8E1B1B"
col_wine_mid   <- "#C13F2B"
col_wine_light <- "#F2B9A3"

col_navy_dark  <- "#083D77"
col_navy_mid   <- "#2C7FB8"
col_navy_light <- "#BFD7EA"

## Reversal score strip: blue scale
score_col_fun <- circlize::colorRamp2(
  c(0.20, 0.35, 0.50, 0.65),
  c(
    "#A7C6DF",
    "#6FA6CC",
    "#2C7FB8",
    "#12508D"
  )
)

## Reversal proportion strip: red-orange scale
prop_col_fun <- circlize::colorRamp2(
  c(0.40, 0.55, 0.70, 0.85),
  c(
    "#F4AE8C",
    "#EA7C58",
    "#D24A34",
    "#A91D1D"
  )
)

## Main gene-level reversal-value heatmap
## Use soft neutral instead of pure white at 0
reversal_col_fun <- circlize::colorRamp2(
  c(-1.5, -0.75, -0.25, 0, 0.25, 0.75, 1.5),
  c(
    col_navy_dark,
    col_navy_mid,
    col_navy_light,
    col_bg0,
    col_wine_light,
    col_wine_mid,
    col_wine_dark
  )
)

## Metric summary plot colors
col_metric_score <- "#12508D"
col_metric_prop  <- "#A91D1D"
col_border <- "#D9D9D9"
col_grid <- "#EFEFEF"
col_strip_fill <- "#F5F5F5"


############################################################
## 5. Helper functions
############################################################
open_pdf_device <- function(file, width, height) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF(
      file = file,
      width = width,
      height = height,
      family = global_font_family
    )
  } else if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = file,
      width = width,
      height = height,
      family = global_font_family
    )
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      family = global_font_family
    )
  }
}

open_tiff_device <- function(file, width, height, dpi = 600) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoTIFF(
      filename = file,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
  } else {
    grDevices::tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      type = if (capabilities("cairo")) "cairo" else "Xlib",
      bg = "white"
    )
  }
}

open_png_device <- function(file, width, height, dpi = 600) {
  grDevices::png(
    filename = file,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    type = "cairo-png",
    bg = "white"
  )
}

clean_text_vector <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

make_key <- function(x) {
  stringr::str_to_lower(stringr::str_trim(as.character(x)))
}

safe_filename <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x <- gsub("[/:*?\"<>|\\\\]", "_", x)
  x <- gsub("\\s+", "_", x)
  x
}

safe_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}

clip_matrix <- function(mat, clip_val = 1.5) {
  m <- mat
  m[m > clip_val] <- clip_val
  m[m < -clip_val] <- -clip_val
  m
}

cluster_cols_within_group <- function(mat, cols) {
  cols <- intersect(cols, colnames(mat))
  if (length(cols) <= 1) return(cols)
  
  submat <- mat[, cols, drop = FALSE]
  
  cor_mat <- suppressWarnings(
    stats::cor(
      submat,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  if (all(is.na(cor_mat))) return(cols)
  
  cor_mat[is.na(cor_mat)] <- 0
  d <- as.dist(1 - cor_mat)
  hc <- hclust(d, method = "average")
  cols[hc$order]
}

theme_pub_cmap <- function(base_size = 13, base_family = global_font_family) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = col_border,
        linewidth = 0.5
      ),
      panel.grid.major = element_line(
        color = col_grid,
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        color = col_text,
        margin = margin(b = 12)
      ),
      axis.text = element_text(color = col_text),
      axis.title = element_text(
        color = col_text,
        face = "bold"
      ),
      axis.ticks = element_line(color = col_border),
      legend.title = element_text(
        color = col_text,
        face = "bold"
      ),
      legend.text = element_text(color = col_text),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = unit(c(5, 8, 5, 5), "mm")
    )
}

save_ggplot_sci <- function(plot_obj, filename_prefix, width, height, dpi = 600) {
  open_pdf_device(
    file = paste0(filename_prefix, ".pdf"),
    width = width,
    height = height
  )
  print(plot_obj)
  dev.off()
  
  open_png_device(
    file = paste0(filename_prefix, ".png"),
    width = width,
    height = height,
    dpi = dpi
  )
  print(plot_obj)
  dev.off()
  
  open_tiff_device(
    file = paste0(filename_prefix, ".tiff"),
    width = width,
    height = height,
    dpi = dpi
  )
  print(plot_obj)
  dev.off()
  
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(
      file = paste0(filename_prefix, ".svg"),
      width = width,
      height = height,
      bg = "white"
    )
    print(plot_obj)
    dev.off()
  }
}


############################################################
## 6. Read drug list and filter matched signatures
############################################################
drug_list <- readr::read_lines(drug_list_file) %>%
  clean_text_vector()

if (length(drug_list) == 0) {
  stop("drug_list.txt 为空。")
}

drug_list_df <- data.frame(
  drug_name = drug_list,
  drug_key = make_key(drug_list),
  stringsAsFactors = FALSE
) %>%
  distinct(drug_key, .keep_all = TRUE)

write.csv(
  drug_list_df,
  file.path(checks_dir, "drug_list_used.csv"),
  row.names = FALSE,
  quote = FALSE
)

matched <- data.table::fread(matched_file)

required_matched_cols <- c(
  "id",
  "pert_iname",
  "cell_iname",
  "pert_idose",
  "pert_itime"
)

missing_matched_cols <- setdiff(required_matched_cols, colnames(matched))

if (length(missing_matched_cols) > 0) {
  stop(
    matched_file,
    " 缺少以下列：",
    paste(missing_matched_cols, collapse = ", ")
  )
}

matched <- matched %>%
  mutate(
    id = as.character(id),
    pert_iname = as.character(pert_iname),
    cell_iname = as.character(cell_iname),
    pert_idose = as.character(pert_idose),
    pert_itime = as.character(pert_itime),
    drug_key = make_key(pert_iname)
  )

matched_not_in_list <- matched %>%
  filter(!drug_key %in% drug_list_df$drug_key)

write.csv(
  matched_not_in_list,
  file.path(checks_dir, "matched_drugs_not_in_drug_list.csv"),
  row.names = FALSE,
  quote = FALSE
)

drug_missing_in_matched <- drug_list_df %>%
  filter(!drug_key %in% matched$drug_key)

write.csv(
  drug_missing_in_matched,
  file.path(checks_dir, "drug_list_missing_in_matched_drugs.csv"),
  row.names = FALSE,
  quote = FALSE
)

## 强制过滤：drug_list.txt 是最终药物过滤依据
matched <- matched %>%
  filter(drug_key %in% drug_list_df$drug_key)

if (nrow(matched) == 0) {
  stop(
    "使用 drug_list.txt 过滤后，",
    matched_file,
    " 中没有可分析的 signature。请检查药物名是否一致。"
  )
}

write.csv(
  matched,
  file.path(checks_dir, "virus_matched_drugs_after_drug_list_filter.csv"),
  row.names = FALSE,
  quote = FALSE
)

duplicated_signature_ids <- matched %>%
  count(id, name = "n") %>%
  filter(n > 1)

write.csv(
  duplicated_signature_ids,
  file.path(checks_dir, "duplicated_signature_ids.csv"),
  row.names = FALSE,
  quote = FALSE
)

matched <- matched %>%
  distinct(id, .keep_all = TRUE)

duplicated_drug_names <- matched %>%
  count(pert_iname, name = "n") %>%
  filter(n > 1)

write.csv(
  duplicated_drug_names,
  file.path(checks_dir, "duplicated_drug_names_in_matched.csv"),
  row.names = FALSE,
  quote = FALSE
)

matched_condition_counts <- matched %>%
  group_by(pert_iname, cell_iname, pert_idose, pert_itime) %>%
  summarise(
    n_signature = n(),
    sig_ids = paste(id, collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_signature), pert_iname)

write.csv(
  matched_condition_counts,
  file.path(checks_dir, "matched_condition_signature_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

sig_ids <- unique(matched$id)

if (length(sig_ids) == 0) {
  stop("matched 文件中没有可分析的 signature id。")
}

cat("病毒标签基因 matched signature 数：", length(sig_ids), "\n")


############################################################
## 7. Read gene annotation and build HDF/HRF map
############################################################
gene_info <- data.table::fread(geneinfo_file)

required_geneinfo_cols <- c("gene_id", "gene_symbol")
missing_geneinfo_cols <- setdiff(required_geneinfo_cols, colnames(gene_info))

if (length(missing_geneinfo_cols) > 0) {
  stop(
    geneinfo_file,
    " 缺少以下列：",
    paste(missing_geneinfo_cols, collapse = ", ")
  )
}

gene_info <- gene_info %>%
  mutate(
    gene_id = as.character(gene_id),
    gene_symbol = as.character(gene_symbol),
    gene_symbol_upper = str_to_upper(gene_symbol)
  )

gene_info_by_symbol <- gene_info %>%
  filter(
    !is.na(gene_symbol_upper),
    gene_symbol_upper != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct(gene_symbol_upper, .keep_all = TRUE)

gene_info_by_id <- gene_info %>%
  filter(
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

symbol_to_entrez <- setNames(
  gene_info_by_symbol$gene_id,
  gene_info_by_symbol$gene_symbol_upper
)

entrez_to_symbol <- setNames(
  gene_info_by_id$gene_symbol,
  gene_info_by_id$gene_id
)

up_genes <- readr::read_lines(up_file) %>%
  clean_text_vector() %>%
  str_to_upper()

down_genes <- readr::read_lines(down_file) %>%
  clean_text_vector() %>%
  str_to_upper()

if (length(up_genes) == 0) {
  stop("virus_up_genes.txt 为空。")
}

if (length(down_genes) == 0) {
  stop("virus_down_genes.txt 为空。")
}

overlap_symbols <- intersect(up_genes, down_genes)

if (length(overlap_symbols) > 0) {
  stop(
    "以下基因同时出现在 virus_up_genes.txt 和 virus_down_genes.txt 中，请先清理输入文件：\n",
    paste(overlap_symbols, collapse = ", ")
  )
}

up_map <- data.frame(
  input_symbol = up_genes,
  gene_id = unname(symbol_to_entrez[up_genes]),
  input_direction = "up",
  host_factor_type = "HDFs",
  column_group = 1,
  input_rank = seq_along(up_genes),
  stringsAsFactors = FALSE
)

down_map <- data.frame(
  input_symbol = down_genes,
  gene_id = unname(symbol_to_entrez[down_genes]),
  input_direction = "down",
  host_factor_type = "HRFs",
  column_group = 2,
  input_rank = seq_along(down_genes),
  stringsAsFactors = FALSE
)

gene_mapping_table <- bind_rows(up_map, down_map)

write.csv(
  gene_mapping_table,
  file.path(checks_dir, "gene_symbol_to_entrez_mapping_all.csv"),
  row.names = FALSE,
  quote = FALSE
)

unmatched_up_genes <- up_map %>%
  filter(is.na(gene_id) | gene_id == "") %>%
  select(input_symbol)

unmatched_down_genes <- down_map %>%
  filter(is.na(gene_id) | gene_id == "") %>%
  select(input_symbol)

write.csv(
  unmatched_up_genes,
  file.path(checks_dir, "unmatched_up_genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  unmatched_down_genes,
  file.path(checks_dir, "unmatched_down_genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

up_mapped <- up_map %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  mutate(gene_id = as.character(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

down_mapped <- down_map %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  mutate(gene_id = as.character(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)

if (nrow(up_mapped) == 0) {
  stop("virus_up_genes.txt 中的基因没有匹配到 geneinfo_beta.txt。")
}

if (nrow(down_mapped) == 0) {
  stop("virus_down_genes.txt 中的基因没有匹配到 geneinfo_beta.txt。")
}

overlap_entrez <- intersect(up_mapped$gene_id, down_mapped$gene_id)

if (length(overlap_entrez) > 0) {
  overlap_symbol <- unique(unname(entrez_to_symbol[overlap_entrez]))
  stop(
    "up/down 基因转换为 Entrez ID 后存在重叠，请检查 geneinfo 映射或输入基因：\n",
    paste(overlap_symbol, collapse = ", ")
  )
}

input_direction_map <- bind_rows(up_mapped, down_mapped) %>%
  select(
    gene_id,
    input_symbol,
    input_direction,
    host_factor_type,
    column_group,
    input_rank
  ) %>%
  mutate(
    gene_id = as.character(gene_id)
  )

direction_check <- input_direction_map %>%
  count(gene_id, name = "n_direction") %>%
  filter(n_direction > 1)

if (nrow(direction_check) > 0) {
  stop("部分 Entrez ID 同时对应 HDFs 和 HRFs，请检查输入文件。")
}

gene_mapping_summary <- data.frame(
  item = c(
    "input_HDF_symbols_from_virus_up",
    "input_HRF_symbols_from_virus_down",
    "mapped_HDF_entrez",
    "mapped_HRF_entrez",
    "total_mapped_entrez"
  ),
  n = c(
    length(up_genes),
    length(down_genes),
    nrow(up_mapped),
    nrow(down_mapped),
    nrow(input_direction_map)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  gene_mapping_summary,
  file.path(checks_dir, "gene_mapping_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("HDFs mapped from virus_up_genes.txt：", nrow(up_mapped), "\n")
cat("HRFs mapped from virus_down_genes.txt：", nrow(down_mapped), "\n")
cat("病毒标签宿主因子总数：", nrow(input_direction_map), "\n")


############################################################
## 8. Core function for one signature
############################################################
compute_one_signature <- function(
    sig_id,
    meta_row,
    gctx_file,
    input_direction_map,
    entrez_to_symbol
) {
  gct <- cmapR::parse_gctx(gctx_file, cid = sig_id)
  
  if (ncol(gct@mat) != 1) {
    stop("sig_id = ", sig_id, " 读取后不是单列矩阵。")
  }
  
  zscores <- gct@mat[, 1]
  gene_ids_in_gctx <- rownames(gct@mat)
  
  if (is.null(gene_ids_in_gctx)) {
    stop("sig_id = ", sig_id, " 的表达矩阵没有行名。")
  }
  
  zscores <- as.numeric(zscores)
  names(zscores) <- as.character(gene_ids_in_gctx)
  
  ranked <- sort(zscores, decreasing = TRUE)
  
  ranked_df <- data.frame(
    gene_id = names(ranked),
    zscore = as.numeric(ranked),
    rank = seq_along(ranked),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      gene_id = as.character(gene_id),
      rank_percent = rank / n(),
      regulation_in_drug = case_when(
        zscore > 0 ~ "up",
        zscore < 0 ~ "down",
        TRUE ~ "neutral"
      ),
      abs_zscore = abs(zscore)
    )
  
  drug_gene_effect <- ranked_df %>%
    inner_join(input_direction_map, by = "gene_id") %>%
    mutate(
      gene_symbol = unname(entrez_to_symbol[gene_id]),
      gene_symbol = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        input_symbol,
        gene_symbol
      ),
      
      ## Core reversal definition
      ## HDFs: drug down-regulation is reversal
      ## HRFs: drug up-regulation is reversal
      reversal_value = case_when(
        input_direction == "up" ~ -zscore,
        input_direction == "down" ~ zscore,
        TRUE ~ NA_real_
      ),
      
      reversed = reversal_value > 0,
      reversal_strength = ifelse(reversal_value > 0, reversal_value, 0),
      same_direction = reversal_value < 0,
      
      drug_name = meta_row$pert_iname,
      sig_id = sig_id,
      cell_id = meta_row$cell_iname,
      pert_dose = meta_row$pert_idose,
      pert_time = meta_row$pert_itime
    ) %>%
    select(
      drug_name,
      sig_id,
      cell_id,
      pert_dose,
      pert_time,
      gene_symbol,
      gene_id,
      input_symbol,
      input_direction,
      host_factor_type,
      zscore,
      abs_zscore,
      regulation_in_drug,
      reversal_value,
      reversed,
      reversal_strength,
      same_direction,
      rank,
      rank_percent,
      column_group,
      input_rank
    ) %>%
    arrange(
      column_group,
      input_rank,
      desc(reversal_strength),
      desc(abs_zscore),
      rank
    )
  
  valid_df <- drug_gene_effect %>%
    filter(!is.na(reversal_value))
  
  n_total_genes <- nrow(valid_df)
  n_reversed <- sum(valid_df$reversed %in% TRUE, na.rm = TRUE)
  
  summary_row <- data.frame(
    drug_name = meta_row$pert_iname,
    sig_id = sig_id,
    cell_id = meta_row$cell_iname,
    pert_dose = meta_row$pert_idose,
    pert_time = meta_row$pert_itime,
    
    n_target_genes_found = nrow(drug_gene_effect),
    
    n_HDF_genes_found = sum(
      drug_gene_effect$host_factor_type == "HDFs",
      na.rm = TRUE
    ),
    
    n_HRF_genes_found = sum(
      drug_gene_effect$host_factor_type == "HRFs",
      na.rm = TRUE
    ),
    
    mean_abs_zscore = ifelse(
      nrow(drug_gene_effect) > 0,
      mean(drug_gene_effect$abs_zscore, na.rm = TRUE),
      NA_real_
    ),
    
    mean_zscore = ifelse(
      nrow(drug_gene_effect) > 0,
      mean(drug_gene_effect$zscore, na.rm = TRUE),
      NA_real_
    ),
    
    n_total_genes_for_reversal = n_total_genes,
    n_reversed = n_reversed,
    
    prop_reversed = ifelse(
      n_total_genes > 0,
      n_reversed / n_total_genes,
      NA_real_
    ),
    
    mean_reversal_strength = ifelse(
      n_total_genes > 0,
      mean(valid_df$reversal_strength, na.rm = TRUE),
      NA_real_
    ),
    
    reversal_score = ifelse(
      n_total_genes > 0,
      sum(valid_df$reversal_strength, na.rm = TRUE) / n_total_genes,
      NA_real_
    ),
    
    stringsAsFactors = FALSE
  )
  
  list(
    gene_table = drug_gene_effect,
    summary_row = summary_row
  )
}


############################################################
## 9. Main loop
############################################################
summary_list <- vector("list", length(sig_ids))
all_gene_effect_list <- vector("list", length(sig_ids))

pb <- txtProgressBar(min = 0, max = length(sig_ids), style = 3)

for (i in seq_along(sig_ids)) {
  sig_id <- sig_ids[i]
  
  tryCatch({
    meta_row <- matched %>%
      filter(id == sig_id) %>%
      slice(1)
    
    res <- compute_one_signature(
      sig_id = sig_id,
      meta_row = meta_row,
      gctx_file = gctx_file,
      input_direction_map = input_direction_map,
      entrez_to_symbol = entrez_to_symbol
    )
    
    summary_list[[i]] <- res$summary_row
    all_gene_effect_list[[i]] <- res$gene_table
    
    if (export_per_signature_csv) {
      out_name <- paste0(
        sprintf("%03d", i),
        "_",
        safe_filename(meta_row$pert_iname),
        "_",
        safe_filename(meta_row$cell_iname),
        "_",
        safe_filename(meta_row$pert_idose),
        "_",
        safe_filename(meta_row$pert_itime),
        "_",
        safe_filename(sig_id),
        ".csv"
      )
      
      readr::write_csv(
        res$gene_table,
        file.path(per_sig_dir, out_name)
      )
    }
  }, error = function(e) {
    message("Error in sig_id ", sig_id, ": ", e$message)
    
    meta_row <- matched %>%
      filter(id == sig_id) %>%
      slice(1)
    
    summary_list[[i]] <<- data.frame(
      drug_name = ifelse(nrow(meta_row) > 0, meta_row$pert_iname, NA),
      sig_id = sig_id,
      cell_id = ifelse(nrow(meta_row) > 0, meta_row$cell_iname, NA),
      pert_dose = ifelse(nrow(meta_row) > 0, meta_row$pert_idose, NA),
      pert_time = ifelse(nrow(meta_row) > 0, meta_row$pert_itime, NA),
      n_target_genes_found = NA,
      n_HDF_genes_found = NA,
      n_HRF_genes_found = NA,
      mean_abs_zscore = NA,
      mean_zscore = NA,
      n_total_genes_for_reversal = NA,
      n_reversed = NA,
      prop_reversed = NA,
      mean_reversal_strength = NA,
      reversal_score = NA,
      stringsAsFactors = FALSE
    )
    
    all_gene_effect_list[[i]] <<- NULL
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)


############################################################
## 10. Summary output
############################################################
drug_gene_effect_summary <- bind_rows(summary_list) %>%
  arrange(
    desc(reversal_score),
    desc(prop_reversed),
    desc(mean_reversal_strength)
  )

if (nrow(drug_gene_effect_summary) == 0) {
  stop("没有生成任何 signature 结果，请检查输入数据。")
}

drug_gene_effect_summary <- drug_gene_effect_summary %>%
  mutate(
    heatmap_label = make.unique(as.character(drug_name), sep = "_")
  )

all_gene_effect_df <- bind_rows(all_gene_effect_list) %>%
  left_join(
    drug_gene_effect_summary %>%
      select(sig_id, heatmap_label),
    by = "sig_id"
  )

if (nrow(all_gene_effect_df) == 0) {
  stop("all_gene_effect_df 为空。")
}

readr::write_csv(
  drug_gene_effect_summary,
  file.path(output_dir, "drug_gene_effect_summary.csv")
)

readr::write_csv(
  all_gene_effect_df,
  file.path(output_dir, "all_drug_gene_effect_long_table.csv")
)

write.csv(
  drug_gene_effect_summary %>%
    select(
      heatmap_label,
      drug_name,
      sig_id,
      cell_id,
      pert_dose,
      pert_time,
      reversal_score,
      prop_reversed,
      n_reversed,
      n_total_genes_for_reversal,
      n_HDF_genes_found,
      n_HRF_genes_found
    ),
  file.path(matrices_dir, "heatmap_row_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)


############################################################
## 11. Build heatmap matrix
############################################################
heatmap_df <- all_gene_effect_df %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(heatmap_label)
  ) %>%
  group_by(heatmap_label, gene_symbol) %>%
  summarise(
    reversal_value = mean(reversal_value, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_df$reversal_value[is.nan(heatmap_df$reversal_value)] <- NA_real_

heatmap_mat <- heatmap_df %>%
  tidyr::pivot_wider(
    names_from = gene_symbol,
    values_from = reversal_value
  ) %>%
  as.data.frame()

if (!"heatmap_label" %in% colnames(heatmap_mat)) {
  stop("无法构建热图矩阵：heatmap_label 列不存在。")
}

rownames(heatmap_mat) <- heatmap_mat$heatmap_label
heatmap_mat$heatmap_label <- NULL
heatmap_mat <- as.matrix(heatmap_mat)

ordered_rows <- drug_gene_effect_summary$heatmap_label
ordered_rows <- ordered_rows[ordered_rows %in% rownames(heatmap_mat)]
heatmap_mat <- heatmap_mat[ordered_rows, , drop = FALSE]

gene_anno_raw <- all_gene_effect_df %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  select(
    gene_symbol,
    input_symbol,
    input_direction,
    host_factor_type,
    column_group,
    input_rank
  ) %>%
  distinct()

gene_anno_check <- gene_anno_raw %>%
  count(gene_symbol, name = "n_type") %>%
  filter(n_type > 1)

if (nrow(gene_anno_check) > 0) {
  stop(
    "以下 gene_symbol 对应多个输入方向或宿主因子类型，请检查输入文件或映射：\n",
    paste(gene_anno_check$gene_symbol, collapse = ", ")
  )
}

gene_anno <- gene_anno_raw %>%
  arrange(column_group, input_rank)

hdf_cols <- gene_anno %>%
  filter(host_factor_type == "HDFs") %>%
  pull(gene_symbol) %>%
  unique()

hrf_cols <- gene_anno %>%
  filter(host_factor_type == "HRFs") %>%
  pull(gene_symbol) %>%
  unique()

ordered_cols <- c(hdf_cols, hrf_cols)
ordered_cols <- intersect(ordered_cols, colnames(heatmap_mat))

if (length(ordered_cols) == 0) {
  stop("热图没有可用基因列。")
}

heatmap_mat <- heatmap_mat[, ordered_cols, drop = FALSE]

annotation_col_df <- gene_anno %>%
  filter(gene_symbol %in% colnames(heatmap_mat)) %>%
  select(
    gene_symbol,
    input_direction,
    host_factor_type,
    column_group,
    input_rank
  ) %>%
  distinct()

annotation_col_df <- annotation_col_df[
  match(colnames(heatmap_mat), annotation_col_df$gene_symbol),
  ,
  drop = FALSE
]

rownames(annotation_col_df) <- annotation_col_df$gene_symbol

heatmap_mat_plot <- clip_matrix(
  heatmap_mat,
  clip_val = heatmap_clip_value
)

write.csv(
  heatmap_mat,
  file.path(matrices_dir, "heatmap_reversal_value_matrix_raw.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  heatmap_mat_plot,
  file.path(matrices_dir, "heatmap_reversal_value_matrix_clipped.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  annotation_col_df,
  file.path(matrices_dir, "heatmap_gene_annotation.csv"),
  row.names = TRUE,
  quote = FALSE
)


############################################################
## 12. Narrow metric matrices
##     Use heatmap_label as row names, matching heatmap_mat_plot
############################################################
score_mat <- matrix(
  drug_gene_effect_summary$reversal_score,
  ncol = 1
)

rownames(score_mat) <- drug_gene_effect_summary$heatmap_label
colnames(score_mat) <- "Reversal score"

prop_mat <- matrix(
  drug_gene_effect_summary$prop_reversed,
  ncol = 1
)

rownames(prop_mat) <- drug_gene_effect_summary$heatmap_label
colnames(prop_mat) <- "Reversal proportion"

score_mat <- score_mat[rownames(heatmap_mat_plot), , drop = FALSE]
prop_mat  <- prop_mat[rownames(heatmap_mat_plot), , drop = FALSE]


############################################################
## 13. Top annotation
##     Same style as residual v12, but labels are HDFs / HRFs
############################################################
annotation_col_df$host_factor_type <- factor(
  annotation_col_df$host_factor_type,
  levels = c("HDFs", "HRFs")
)

top_ha <- HeatmapAnnotation(
  `Host factor type` = annotation_col_df$host_factor_type,
  col = list(
    `Host factor type` = c(
      "HDFs" = col_hdf,
      "HRFs" = col_hrf
    )
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(5.0, "mm"),
  gp = gpar(col = NA),
  show_legend = FALSE
)


############################################################
## 14. Legends
##     Same right-side legend style as residual v12
############################################################
lgd_host <- Legend(
  title = "Host factor type",
  labels = c("HDFs", "HRFs"),
  legend_gp = gpar(fill = c(col_hdf, col_hrf), col = NA),
  title_gp = gpar(
    fontsize = 12,
    fontface = "bold",
    family = global_font_family
  ),
  labels_gp = gpar(
    fontsize = 11,
    family = global_font_family
  ),
  grid_width = unit(5.0, "mm"),
  grid_height = unit(5.0, "mm")
)

lgd_value <- Legend(
  title = "Reversal value",
  col_fun = reversal_col_fun,
  at = c(-1.5, 0, 1.5),
  labels = c("-1.5", "0", "1.5"),
  title_gp = gpar(
    fontsize = 12,
    fontface = "bold",
    family = global_font_family
  ),
  labels_gp = gpar(
    fontsize = 11,
    family = global_font_family
  ),
  legend_height = unit(45, "mm")
)

legend_pack <- packLegend(
  lgd_host,
  lgd_value,
  direction = "vertical",
  gap = unit(10, "mm")
)


############################################################
## 15. Build integrated heatmap object
##     Style and layout copied from residual v12
############################################################
n_drug <- nrow(heatmap_mat_plot)

row_font_size <- ifelse(
  n_drug <= 8, 14,
  ifelse(n_drug <= 20, 10.5, 8.5)
)

body_height <- if (n_drug <= 8) {
  unit(86, "mm")
} else if (n_drug <= 20) {
  unit(max(95, n_drug * 5.6), "mm")
} else {
  unit(max(115, n_drug * 4.0), "mm")
}

## Same metric strip width as residual v12
metric_strip_width <- unit(5.8, "mm")

ht_score <- Heatmap(
  score_mat,
  name = "Reversal score",
  col = score_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = row_font_size,
    family = global_font_family
  ),
  row_names_max_width = unit(46, "mm"),
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 90,
  column_names_centered = TRUE,
  column_names_gp = gpar(
    fontsize = 10.5,
    fontface = "bold",
    family = global_font_family
  ),
  column_names_max_height = unit(40, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  border = FALSE,
  width = metric_strip_width,
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_prop <- Heatmap(
  prop_mat,
  name = "Reversal proportion",
  col = prop_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 90,
  column_names_centered = TRUE,
  column_names_gp = gpar(
    fontsize = 10.5,
    fontface = "bold",
    family = global_font_family
  ),
  column_names_max_height = unit(48, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  border = FALSE,
  width = metric_strip_width,
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_gene <- Heatmap(
  heatmap_mat_plot,
  name = "Reversal value",
  col = reversal_col_fun,
  top_annotation = top_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  rect_gp = gpar(col = NA),
  border = FALSE,
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 4,
  na_col = "#EFE8DF",
  width = unit(205, "mm"),
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_final <- ht_score + ht_prop + ht_gene


############################################################
## 16. Capture heatmap grob
##     Same bottom padding style as residual v12
############################################################
heatmap_grob <- grid.grabExpr(
  draw(
    ht_final,
    newpage = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    merge_legends = FALSE,
    padding = unit(c(3, 3, 28, 4), "mm")
  )
)


############################################################
## 17. Final draw function
##     Same title / heatmap / right legend position as residual v12
############################################################
draw_integrated_figure <- function(
    heatmap_grob,
    legend_pack,
    title_text
) {
  grid.newpage()
  
  lay <- grid.layout(
    nrow = 4,
    ncol = 5,
    heights = unit.c(
      unit(0.82, "in"),   ## title row
      unit(0.04, "in"),   ## tiny spacer
      unit(1.00, "null"), ## body row
      unit(0.34, "in")    ## bottom margin
    ),
    widths = unit.c(
      unit(0.16, "in"),   ## left outer margin
      unit(12.20, "in"),  ## heatmap panel
      unit(0.18, "in"),   ## gap
      unit(2.35, "in"),   ## legend panel
      unit(0.14, "in")    ## right margin
    )
  )
  
  pushViewport(viewport(layout = lay, name = "root"))
  
  ## Title
  pushViewport(viewport(
    layout.pos.row = 1,
    layout.pos.col = 1:5,
    clip = "off"
  ))
  
  grid.text(
    label = title_text,
    x = 0.5,
    y = 0.54,
    gp = gpar(
      fontsize = 19,
      fontface = "bold",
      family = global_font_family,
      col = col_text
    )
  )
  
  popViewport()
  
  ## Heatmap panel
  pushViewport(viewport(
    layout.pos.row = 3,
    layout.pos.col = 2,
    just = c("center", "top"),
    clip = "off"
  ))
  
  grid.draw(heatmap_grob)
  
  popViewport()
  
  ## Legend panel
  pushViewport(viewport(
    layout.pos.row = 3,
    layout.pos.col = 4,
    just = c("left", "top"),
    clip = "off"
  ))
  
  draw(
    legend_pack,
    x = unit(0, "npc"),
    y = unit(1, "npc"),
    just = c("left", "top")
  )
  
  popViewport()
  
  popViewport()
}


############################################################
## 18. Export main SCI figure
############################################################
out_pdf <- file.path(
  plots_dir,
  "Virus_tagged_host_factor_reversal_profiles_integrated_heatmap_SCI.pdf"
)

out_png <- file.path(
  plots_dir,
  "Virus_tagged_host_factor_reversal_profiles_integrated_heatmap_SCI.png"
)

out_tiff <- file.path(
  plots_dir,
  "Virus_tagged_host_factor_reversal_profiles_integrated_heatmap_SCI.tiff"
)

out_svg <- file.path(
  plots_dir,
  "Virus_tagged_host_factor_reversal_profiles_integrated_heatmap_SCI.svg"
)

## PDF
open_pdf_device(out_pdf, width = pdf_width, height = pdf_height)
draw_integrated_figure(
  heatmap_grob = heatmap_grob,
  legend_pack = legend_pack,
  title_text = main_title
)
dev.off()

## PNG
open_png_device(out_png, width = pdf_width, height = pdf_height, dpi = plot_dpi)
draw_integrated_figure(
  heatmap_grob = heatmap_grob,
  legend_pack = legend_pack,
  title_text = main_title
)
dev.off()

## TIFF
open_tiff_device(out_tiff, width = pdf_width, height = pdf_height, dpi = plot_dpi)
draw_integrated_figure(
  heatmap_grob = heatmap_grob,
  legend_pack = legend_pack,
  title_text = main_title
)
dev.off()

## SVG
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(
    file = out_svg,
    width = pdf_width,
    height = pdf_height,
    bg = "white"
  )
  draw_integrated_figure(
    heatmap_grob = heatmap_grob,
    legend_pack = legend_pack,
    title_text = main_title
  )
  dev.off()
}


############################################################
## 19. SCI-style metric summary plot
##     Backup / supplementary figure
############################################################
rank_df <- drug_gene_effect_summary %>%
  filter(
    !is.na(prop_reversed),
    !is.na(reversal_score),
    !is.na(heatmap_label)
  ) %>%
  arrange(desc(reversal_score)) %>%
  mutate(
    drug_name_plot = factor(
      heatmap_label,
      levels = rev(heatmap_label)
    )
  )

if (nrow(rank_df) > 0) {
  
  metric_df <- rank_df %>%
    select(
      drug_name_plot,
      heatmap_label,
      reversal_score,
      prop_reversed
    ) %>%
    tidyr::pivot_longer(
      cols = c(reversal_score, prop_reversed),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = case_when(
        metric == "reversal_score" ~ "Reversal score",
        metric == "prop_reversed" ~ "Reversal proportion",
        TRUE ~ metric
      ),
      metric = factor(
        metric,
        levels = c("Reversal score", "Reversal proportion")
      ),
      value_label = sprintf("%.3f", value)
    )
  
  x_max <- max(metric_df$value, na.rm = TRUE)
  
  if (!is.finite(x_max) || is.na(x_max) || x_max <= 0) {
    x_max <- 1
  } else {
    x_max <- x_max * 1.25
  }
  
  p_metric <- ggplot(
    metric_df,
    aes(
      x = value,
      y = drug_name_plot
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = value,
        y = drug_name_plot,
        yend = drug_name_plot
      ),
      linewidth = 0.65,
      color = "#BDBDBD"
    ) +
    geom_point(
      aes(color = metric),
      size = 4.0,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = value_label),
      hjust = -0.15,
      size = 3.6,
      family = global_font_family
    ) +
    facet_wrap(
      ~ metric,
      nrow = 1,
      scales = "fixed"
    ) +
    scale_color_manual(
      values = c(
        "Reversal score" = col_metric_score,
        "Reversal proportion" = col_metric_prop
      )
    ) +
    scale_x_continuous(
      limits = c(0, x_max),
      expand = expansion(mult = c(0, 0.10))
    ) +
    theme_pub_cmap(base_size = 13) +
    labs(
      title = "Drug reversal metrics for virus-tagged host factors",
      x = NULL,
      y = "Drug"
    ) +
    theme(
      strip.background = element_rect(
        fill = col_strip_fill,
        color = col_border,
        linewidth = 0.4
      ),
      strip.text = element_text(
        size = 12,
        face = "bold",
        color = col_text
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(
        size = 10.5,
        color = col_text
      ),
      axis.title.y = element_text(
        size = 12,
        face = "bold"
      ),
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )
    )
  
  metric_height <- max(4.8, nrow(rank_df) * 0.46 + 2.0)
  
  save_ggplot_sci(
    plot_obj = p_metric,
    filename_prefix = file.path(
      plots_dir,
      "Drug_reversal_metric_summary_virus_tagged_host_factors_SCI"
    ),
    width = 10.8,
    height = metric_height,
    dpi = plot_dpi
  )
}


############################################################
## 20. README
############################################################
readme_txt <- c(
  "Virus-tagged host-factor drug reversal analysis",
  "",
  "Working directory:",
  "/users/hzhang1/project/Cmap",
  "",
  "Input files:",
  paste0("1. ", gctx_file),
  paste0("2. ", geneinfo_file),
  paste0("3. ", drug_list_file),
  paste0("4. ", matched_file),
  paste0("5. ", up_file),
  paste0("6. ", down_file),
  "",
  "Output directory:",
  output_dir,
  "",
  "Biological convention:",
  "virus_up_genes.txt   -> HDFs",
  "virus_down_genes.txt -> HRFs",
  "",
  "Reversal value definition:",
  "HDFs: reversal_value = -zscore",
  "HRFs: reversal_value =  zscore",
  "reversal_value > 0 indicates drug-induced reversal",
  "reversal_value < 0 indicates same-direction perturbation",
  "",
  "Reversal strength:",
  "reversal_strength = max(reversal_value, 0)",
  "",
  "Reversal score:",
  "reversal_score = sum(reversal_strength) / n_total_genes_for_reversal",
  "",
  "Main SCI figure:",
  "plots/Virus_tagged_host_factor_reversal_profiles_integrated_heatmap_SCI.pdf/png/tiff/svg",
  "",
  "Main figure style:",
  "Same plotting style as residual_gene_reversal_integrated_heatmap_FINAL_SCI_v12.R",
  "",
  "Figure components:",
  "1. Left narrow strip: Reversal score",
  "2. Left narrow strip: Reversal proportion",
  "3. Main heatmap: gene-level reversal value clipped to [-1.5, 1.5]",
  "4. Top annotation: HDFs / HRFs",
  "5. Right legends: Host factor type and Reversal value",
  "",
  "Main output files:",
  "1. drug_gene_effect_summary.csv",
  "2. all_drug_gene_effect_long_table.csv",
  "3. matrices/heatmap_reversal_value_matrix_raw.csv",
  "4. matrices/heatmap_reversal_value_matrix_clipped.csv",
  "5. matrices/heatmap_gene_annotation.csv",
  "6. input_checks/",
  "7. per_signature_csv/",
  "",
  "Caption note:",
  "The two left-side columns show drug-level reversal score and reversal proportion,",
  "with darker colors indicating higher values."
)

writeLines(
  readme_txt,
  con = file.path(output_dir, "README_results.txt")
)


############################################################
## 21. Done
############################################################
end_time <- Sys.time()

run_time <- round(
  as.numeric(
    difftime(end_time, start_time, units = "mins")
  ),
  2
)

cat("\n========================================\n")
cat("病毒标签宿主因子药物逆转分析完成！\n")
cat("工作目录：", getwd(), "\n", sep = "")
cat("输入 matched 文件：", matched_file, "\n", sep = "")
cat("药物过滤文件：", drug_list_file, "\n", sep = "")
cat("输出目录：", output_dir, "\n", sep = "")
cat("图形目录：", plots_dir, "\n", sep = "")
cat("矩阵目录：", matrices_dir, "\n", sep = "")
cat("单 signature 结果目录：", per_sig_dir, "\n", sep = "")
cat("输入检查目录：", checks_dir, "\n", sep = "")
cat("主图 PDF：", out_pdf, "\n", sep = "")
cat("主图 PNG：", out_png, "\n", sep = "")
cat("主图 TIFF：", out_tiff, "\n", sep = "")
if (requireNamespace("svglite", quietly = TRUE)) {
  cat("主图 SVG：", out_svg, "\n", sep = "")
}
cat("总运行时间：", run_time, " 分钟\n", sep = "")
cat("========================================\n")




############################################################
## residual_gene_reversal_integrated_heatmap_FINAL_SCI_v12.R
##
## Purpose:
##   Final SCI-style integrated residual-gene reversal heatmap
##   with:
##     1. left narrow Reversal score strip
##     2. left narrow Reversal proportion strip
##     3. gene-level reversal-value heatmap
##     4. top gene-set annotation
##     5. fixed right legend panel
##
## 
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")
set.seed(123)

############################################################
## 1. Load packages
############################################################
required_pkgs <- c(
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ComplexHeatmap",
  "circlize",
  "grid"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "缺少 R 包：",
    paste(missing_pkgs, collapse = ", "),
    "\n请先在 r45 环境中安装这些包。"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

ht_opt$message <- FALSE

############################################################
## 2. Paths
##    路径保持不变
############################################################
setwd("/users/hzhang1/project/Cmap")

output_dir   <- "/users/hzhang1/project/Cmap/drug_reversal_residual_genes"
plots_dir    <- file.path(output_dir, "plots")
tables_dir   <- file.path(output_dir, "tables")
matrices_dir <- file.path(output_dir, "matrices")

summary_file <- file.path(output_dir, "drug_gene_effect_summary.csv")
long_file    <- file.path(output_dir, "all_drug_gene_effect_long_table.csv")

if (!file.exists(summary_file)) {
  stop("Missing input file: ", summary_file)
}
if (!file.exists(long_file)) {
  stop("Missing input file: ", long_file)
}

## 如果文件系统仍是只读，这里直接给出清晰报错
if (file.access(output_dir, 2) != 0) {
  stop(
    "输出目录当前不可写：\n",
    output_dir,
    "\n请先恢复该路径写权限后再运行。"
  )
}

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrices_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Global settings
############################################################
global_font_family <- "sans"
main_title <- "Residual-gene reversal profiles of candidate drugs"
heatmap_clip_value <- 1.5

pdf_width  <- 15.8
pdf_height <- 8.2
plot_dpi   <- 600

############################################################
## 4. Color system
############################################################
## Core palette
col_text <- "#222222"
col_bg0  <- "#F2ECE4"   # not pure white, softer center

## Gene-set annotation
col_residual_up   <- "#FF4D2D"   # brighter than before
col_residual_down <- "#2C7FB8"
col_both          <- "#7A7A7A"
col_unknown       <- "#BDBDBD"

## Heatmap main colors: wine red <-> soft neutral <-> navy blue
col_wine_dark  <- "#8E1B1B"
col_wine_mid   <- "#C13F2B"
col_wine_light <- "#F2B9A3"

col_navy_dark  <- "#083D77"
col_navy_mid   <- "#2C7FB8"
col_navy_light <- "#BFD7EA"

## Reversal score strip (blue scale)
score_col_fun <- circlize::colorRamp2(
  c(0.20, 0.35, 0.50, 0.65),
  c(
    "#A7C6DF",
    "#6FA6CC",
    "#2C7FB8",
    "#12508D"
  )
)

## Reversal proportion strip (red-orange scale)
prop_col_fun <- circlize::colorRamp2(
  c(0.40, 0.55, 0.70, 0.85),
  c(
    "#F4AE8C",
    "#EA7C58",
    "#D24A34",
    "#A91D1D"
  )
)

## Main gene-level reversal-value heatmap
## Use soft neutral instead of pure white at 0
reversal_col_fun <- circlize::colorRamp2(
  c(-1.5, -0.75, -0.25, 0, 0.25, 0.75, 1.5),
  c(
    col_navy_dark,
    col_navy_mid,
    col_navy_light,
    col_bg0,
    col_wine_light,
    col_wine_mid,
    col_wine_dark
  )
)

############################################################
## 5. Helper functions
############################################################
open_pdf_device <- function(file, width, height) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoPDF(
      file = file,
      width = width,
      height = height,
      family = global_font_family
    )
  } else if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = file,
      width = width,
      height = height,
      family = global_font_family
    )
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      family = global_font_family
    )
  }
}

open_tiff_device <- function(file, width, height, dpi = 600) {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoTIFF(
      filename = file,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw",
      bg = "white"
    )
  } else {
    grDevices::tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      type = if (capabilities("cairo")) "cairo" else "Xlib",
      bg = "white"
    )
  }
}

clip_matrix <- function(mat, clip_val = 1.5) {
  m <- mat
  m[m > clip_val] <- clip_val
  m[m < -clip_val] <- -clip_val
  m
}

cluster_cols_within_group <- function(mat, cols) {
  cols <- intersect(cols, colnames(mat))
  if (length(cols) <= 1) return(cols)
  
  submat <- mat[, cols, drop = FALSE]
  
  cor_mat <- suppressWarnings(
    stats::cor(
      submat,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  if (all(is.na(cor_mat))) return(cols)
  
  cor_mat[is.na(cor_mat)] <- 0
  d <- as.dist(1 - cor_mat)
  hc <- hclust(d, method = "average")
  cols[hc$order]
}

############################################################
## 6. Read input tables
############################################################
drug_gene_effect_summary <- readr::read_csv(
  summary_file,
  show_col_types = FALSE
)

all_gene_effect_df <- readr::read_csv(
  long_file,
  show_col_types = FALSE
)

required_summary_cols <- c(
  "drug_name",
  "sig_id",
  "reversal_score",
  "prop_reversed",
  "n_reversed",
  "n_total_genes_for_reversal"
)

required_long_cols <- c(
  "drug_name",
  "sig_id",
  "gene_symbol",
  "gene_set",
  "reversal_value",
  "reversed",
  "reversal_strength"
)

missing_summary_cols <- setdiff(required_summary_cols, colnames(drug_gene_effect_summary))
missing_long_cols    <- setdiff(required_long_cols, colnames(all_gene_effect_df))

if (length(missing_summary_cols) > 0) {
  stop("Summary file missing columns: ", paste(missing_summary_cols, collapse = ", "))
}
if (length(missing_long_cols) > 0) {
  stop("Long file missing columns: ", paste(missing_long_cols, collapse = ", "))
}

############################################################
## 7. Drug-level summary
############################################################
drug_summary <- drug_gene_effect_summary %>%
  mutate(
    reversal_score = as.numeric(reversal_score),
    prop_reversed  = as.numeric(prop_reversed),
    n_reversed = as.numeric(n_reversed),
    n_total_genes_for_reversal = as.numeric(n_total_genes_for_reversal)
  ) %>%
  group_by(drug_name) %>%
  summarise(
    reversal_score = mean(reversal_score, na.rm = TRUE),
    reversal_proportion = mean(prop_reversed, na.rm = TRUE),
    n_reversed = sum(n_reversed, na.rm = TRUE),
    n_total_genes_for_reversal = sum(n_total_genes_for_reversal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(reversal_score), desc(reversal_proportion))

drug_order <- drug_summary$drug_name

write.csv(
  drug_summary,
  file.path(tables_dir, "drug_level_reversal_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 8. Gene-level matrix
############################################################
all_gene_effect_df <- all_gene_effect_df %>%
  filter(drug_name %in% drug_order) %>%
  mutate(
    reversal_value    = as.numeric(reversal_value),
    reversal_strength = as.numeric(reversal_strength)
  )

gene_level_df <- all_gene_effect_df %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(reversal_value)
  ) %>%
  group_by(drug_name, gene_symbol) %>%
  summarise(
    reversal_value = mean(reversal_value, na.rm = TRUE),
    .groups = "drop"
  )

gene_level_df$reversal_value[is.nan(gene_level_df$reversal_value)] <- NA_real_

gene_level_mat <- gene_level_df %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = reversal_value
  ) %>%
  as.data.frame()

rownames(gene_level_mat) <- gene_level_mat$drug_name
gene_level_mat$drug_name <- NULL
gene_level_mat <- as.matrix(gene_level_mat)

gene_level_mat <- gene_level_mat[
  drug_order[drug_order %in% rownames(gene_level_mat)],
  ,
  drop = FALSE
]

if (nrow(gene_level_mat) == 0 || ncol(gene_level_mat) == 0) {
  stop("Gene-level reversal matrix is empty.")
}

############################################################
## 9. Gene-set annotation and column order
############################################################
gene_anno_df <- all_gene_effect_df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(
    raw_gene_set = case_when(
      any(gene_set == "up")   & any(gene_set == "down") ~ "both",
      any(gene_set == "up")   ~ "up",
      any(gene_set == "down") ~ "down",
      TRUE ~ "unknown"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    gene_set_display = case_when(
      raw_gene_set == "up" ~ "Residual up",
      raw_gene_set == "down" ~ "Residual down",
      raw_gene_set == "both" ~ "Both",
      TRUE ~ "Unknown"
    )
  ) %>%
  as.data.frame()

rownames(gene_anno_df) <- gene_anno_df$gene_symbol

annotation_col_df <- data.frame(
  gene_symbol = colnames(gene_level_mat),
  gene_set_display = gene_anno_df[colnames(gene_level_mat), "gene_set_display"],
  stringsAsFactors = FALSE
)

annotation_col_df$gene_set_display[is.na(annotation_col_df$gene_set_display)] <- "Unknown"
rownames(annotation_col_df) <- annotation_col_df$gene_symbol

up_cols      <- rownames(annotation_col_df)[annotation_col_df$gene_set_display == "Residual up"]
down_cols    <- rownames(annotation_col_df)[annotation_col_df$gene_set_display == "Residual down"]
both_cols    <- rownames(annotation_col_df)[annotation_col_df$gene_set_display == "Both"]
unknown_cols <- rownames(annotation_col_df)[annotation_col_df$gene_set_display == "Unknown"]

ordered_cols <- c(
  cluster_cols_within_group(gene_level_mat, up_cols),
  cluster_cols_within_group(gene_level_mat, down_cols),
  cluster_cols_within_group(gene_level_mat, both_cols),
  cluster_cols_within_group(gene_level_mat, unknown_cols)
)

ordered_cols <- intersect(ordered_cols, colnames(gene_level_mat))

gene_level_mat <- gene_level_mat[, ordered_cols, drop = FALSE]
annotation_col_df <- annotation_col_df[ordered_cols, , drop = FALSE]

gene_level_mat_plot <- clip_matrix(gene_level_mat, clip_val = heatmap_clip_value)

write.csv(
  gene_level_mat,
  file.path(matrices_dir, "integrated_gene_level_reversal_value_matrix_raw.csv"),
  quote = FALSE
)

write.csv(
  gene_level_mat_plot,
  file.path(matrices_dir, "integrated_gene_level_reversal_value_matrix_clipped.csv"),
  quote = FALSE
)

write.csv(
  annotation_col_df,
  file.path(matrices_dir, "integrated_gene_set_annotation.csv"),
  row.names = TRUE,
  quote = FALSE
)

############################################################
## 10. Narrow metric matrices
############################################################
score_mat <- matrix(drug_summary$reversal_score, ncol = 1)
rownames(score_mat) <- drug_summary$drug_name
colnames(score_mat) <- "Reversal score"

prop_mat <- matrix(drug_summary$reversal_proportion, ncol = 1)
rownames(prop_mat) <- drug_summary$drug_name
colnames(prop_mat) <- "Reversal proportion"

score_mat <- score_mat[rownames(gene_level_mat_plot), , drop = FALSE]
prop_mat  <- prop_mat[rownames(gene_level_mat_plot), , drop = FALSE]

############################################################
## 11. Top annotation
############################################################
annotation_col_df$gene_set_display <- factor(
  annotation_col_df$gene_set_display,
  levels = c("Residual up", "Residual down", "Both", "Unknown")
)

top_ha <- HeatmapAnnotation(
  `Gene set` = annotation_col_df$gene_set_display,
  col = list(
    `Gene set` = c(
      "Residual up"   = col_residual_up,
      "Residual down" = col_residual_down,
      "Both"          = col_both,
      "Unknown"       = col_unknown
    )
  ),
  show_annotation_name = FALSE,
  simple_anno_size = unit(5.0, "mm"),
  gp = gpar(col = NA),
  show_legend = FALSE
)

############################################################
## 12. Legends
############################################################
lgd_gene <- Legend(
  title = "Gene set",
  labels = c("Residual up genes", "Residual down genes"),
  legend_gp = gpar(fill = c(col_residual_up, col_residual_down), col = NA),
  title_gp = gpar(fontsize = 12, fontface = "bold", family = global_font_family),
  labels_gp = gpar(fontsize = 11, family = global_font_family),
  grid_width = unit(5.0, "mm"),
  grid_height = unit(5.0, "mm")
)

lgd_value <- Legend(
  title = "Reversal value",
  col_fun = reversal_col_fun,
  at = c(-1.5, 0, 1.5),
  labels = c("-1.5", "0", "1.5"),
  title_gp = gpar(fontsize = 12, fontface = "bold", family = global_font_family),
  labels_gp = gpar(fontsize = 11, family = global_font_family),
  legend_height = unit(45, "mm")
)

legend_pack <- packLegend(
  lgd_gene,
  lgd_value,
  direction = "vertical",
  gap = unit(10, "mm")
)

############################################################
## 13. Build heatmap object
############################################################
n_drug <- nrow(gene_level_mat_plot)

row_font_size <- ifelse(
  n_drug <= 8, 14,
  ifelse(n_drug <= 20, 10.5, 8.5)
)

body_height <- if (n_drug <= 8) {
  unit(86, "mm")
} else if (n_drug <= 20) {
  unit(max(95, n_drug * 5.6), "mm")
} else {
  unit(max(115, n_drug * 4.0), "mm")
}

## 左侧两个色柱稍微加宽
metric_strip_width <- unit(5.8, "mm")

ht_score <- Heatmap(
  score_mat,
  name = "Reversal score",
  col = score_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(
    fontsize = row_font_size,
    family = global_font_family
  ),
  row_names_max_width = unit(46, "mm"),
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 90,
  column_names_centered = TRUE,
  column_names_gp = gpar(
    fontsize = 10.5,
    fontface = "bold",
    family = global_font_family
  ),
  column_names_max_height = unit(40, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  border = FALSE,
  width = metric_strip_width,
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_prop <- Heatmap(
  prop_mat,
  name = "Reversal proportion",
  col = prop_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_rot = 90,
  column_names_centered = TRUE,
  column_names_gp = gpar(
    fontsize = 10.5,
    fontface = "bold",
    family = global_font_family
  ),
  column_names_max_height = unit(48, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.8),
  border = FALSE,
  width = metric_strip_width,
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_gene <- Heatmap(
  gene_level_mat_plot,
  name = "Reversal value",
  col = reversal_col_fun,
  top_annotation = top_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  rect_gp = gpar(col = NA),
  border = FALSE,
  use_raster = TRUE,
  raster_device = "png",
  raster_quality = 4,
  na_col = "#EFE8DF",
  width = unit(205, "mm"),
  height = body_height,
  show_heatmap_legend = FALSE
)

ht_final <- ht_score + ht_prop + ht_gene

############################################################
## 14. Capture heatmap grob
##     bottom padding increased to avoid clipping
############################################################
heatmap_grob <- grid.grabExpr(
  draw(
    ht_final,
    newpage = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    merge_legends = FALSE,
    padding = unit(c(3, 3, 28, 4), "mm")
  )
)

############################################################
## 15. Final draw function
############################################################
draw_integrated_figure <- function(
    heatmap_grob,
    legend_pack,
    title_text
) {
  grid.newpage()
  
  ## 整体布局：
  ## 标题更靠上，热图整体上移，右侧单独留足图例空间
  lay <- grid.layout(
    nrow = 4,
    ncol = 5,
    heights = unit.c(
      unit(0.82, "in"),  ## title row
      unit(0.04, "in"),  ## tiny spacer
      unit(1.00, "null"),## body row
      unit(0.34, "in")   ## bottom margin for vertical labels
    ),
    widths = unit.c(
      unit(0.16, "in"),  ## left outer margin
      unit(12.20, "in"), ## heatmap panel
      unit(0.18, "in"),  ## gap
      unit(2.35, "in"),  ## legend panel
      unit(0.14, "in")   ## right margin
    )
  )
  
  pushViewport(viewport(layout = lay, name = "root"))
  
  ## Title
  pushViewport(viewport(
    layout.pos.row = 1,
    layout.pos.col = 1:5,
    clip = "off"
  ))
  grid.text(
    label = title_text,
    x = 0.5,
    y = 0.54,
    gp = gpar(
      fontsize = 19,
      fontface = "bold",
      family = global_font_family,
      col = col_text
    )
  )
  popViewport()
  
  ## Heatmap panel
  pushViewport(viewport(
    layout.pos.row = 3,
    layout.pos.col = 2,
    just = c("center", "top"),
    clip = "off"
  ))
  grid.draw(heatmap_grob)
  popViewport()
  
  ## Legend panel
  pushViewport(viewport(
    layout.pos.row = 3,
    layout.pos.col = 4,
    just = c("left", "top"),
    clip = "off"
  ))
  draw(
    legend_pack,
    x = unit(0, "npc"),
    y = unit(1, "npc"),
    just = c("left", "top")
  )
  popViewport()
  
  popViewport()
}

############################################################
## 16. Export
############################################################
out_pdf <- file.path(
  plots_dir,
  "Residual_gene_reversal_profiles_of_candidate_drugs_integrated_heatmap_SCI.pdf"
)

out_tiff <- file.path(
  plots_dir,
  "Residual_gene_reversal_profiles_of_candidate_drugs_integrated_heatmap_SCI.tiff"
)

## PDF
open_pdf_device(out_pdf, width = pdf_width, height = pdf_height)
draw_integrated_figure(
  heatmap_grob = heatmap_grob,
  legend_pack = legend_pack,
  title_text = main_title
)
dev.off()

## TIFF
open_tiff_device(out_tiff, width = pdf_width, height = pdf_height, dpi = plot_dpi)
draw_integrated_figure(
  heatmap_grob = heatmap_grob,
  legend_pack = legend_pack,
  title_text = main_title
)
dev.off()

############################################################
## 17. Done
############################################################
cat("\n============================================================\n")
cat("Integrated residual-gene reversal heatmap finished.\n")
cat("PDF : ", out_pdf, "\n", sep = "")
cat("TIFF: ", out_tiff, "\n", sep = "")
cat("============================================================\n")
