# --------------------- 环境设置 ---------------------
options(warn = -1)  # 关闭警告（调试时建议开启警告）
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(scDblFinder)   # 双细胞检测
library(SingleR)       # 细胞类型注释
library(celldex)       # 细胞类型参考数据库
library(future)
library(future.apply)
library(AnnotationDbi)
library(RColorBrewer)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)          # 用于 KEGG 可视化
library(ggplot2)
library(ggdendro)
library(dplyr)
library(tidyr)
library(grid)
library(conflicted)
library(data.table) #用于在marker分析中得到真正想要的“长表”（每行一个基因 marker）
library(dorothea) # DoRothEA / VIPER
library(viper) # DoRothEA / VIPER

plan("multicore", workers = 4)  # 启用并行计算
options(future.globals.maxSize = 8000 * 1024^2) # 设置内存限制为8GB

# 设置 ggplot2 默认字体为 SimHei（确保系统已安装该字体）
theme_set(theme_gray(base_family = "SimHei"))
#
#
# # --------------------- 1. 读取多个组数据并构建 Seurat 对象 ---------------------
# cat("1. 开始【读取多个组别数据并构建Seurat对象】\n")
# group_paths <- list(control1 = "./data/control1/",
#                     s1_hIL = "./data/s1-hIL/",
#                     s40kD = "./data/s40kD/")
# seurat_list <- list()
# for (grp in names(group_paths)) {
#   path <- group_paths[[grp]]
#   cat("正在加载组别：", grp, "，数据路径：", path, "\n")
#
#   # 使用 Seurat 内置函数 ReadMtx 读取 matrix.mtx、features.tsv 和 barcodes.tsv
#   counts <- ReadMtx(
#     mtx = paste0(path, "matrix.mtx"),
#     features = paste0(path, "features.tsv"),
#     cells = paste0(path, "barcodes.tsv")
#   )
#   # 构建 Seurat 对象
#   seu_obj <- CreateSeuratObject(counts = counts, project = grp, min.cells = 3, min.features = 0)
#   seu_obj$group <- grp
#   seurat_list[[grp]] <- seu_obj
# }
# # 合并各组对象
# combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "CombinedData")
# cat("1. 完成【读取多个组别数据并构建Seurat对象】\n")
# saveRDS(combined, file = "./temp/combined_raw_data_1.rds")
# combined <- readRDS("./temp/combined_raw_data_1.rds")
#
# # --------------------- 2. 数据过滤 ---------------------
# cat("2. 开始【数据过滤】\n")
# # 基因过滤：保留表达大于0的细胞数大于3的基因
# counts_matrix <- LayerData(combined, assay = "RNA", layer = "counts")
# genes_cell_counts <- rowSums(counts_matrix > 0)
# genes_to_keep <- names(genes_cell_counts[genes_cell_counts > 3])
# cat("基因过滤前：", nrow(combined), "个基因\n")
# combined <- subset(combined, features = genes_to_keep)
# cat("基因过滤后：", nrow(combined), "个基因\n")
# # 计算线粒体基因比例
# combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
# # 血红蛋白基因
# combined[["percent.hb"]] <- PercentageFeatureSet(combined, pattern = "^Hba-|^Hbb-")
# # 核糖体基因
# combined[["percent.rb"]] <- PercentageFeatureSet(combined, pattern = "^Rps|^Rp")
# # # 风琴图，后续 QC 阈值提供参考
# # p1= VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),layer = "counts",pt.size = 0.1,raster = FALSE, ncol = 5)
# # p1_group= VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),layer = "counts",group.by = "group",pt.size = 0.1,raster = FALSE, ncol = 5)
# # p_qc <- p1 / p1_group
# # ggsave(filename = "temp/QC_violin_overall_and_group.pdf",plot = p_qc, width = 14, height = 10, dpi = 300)
# # 后续 QC 阈值提供参考，决定阈值的关键证据，nCount_RNA vs percent.mt（坏细胞通常：低 nCount + 高 mt）；nCount_RNA vs nFeature_RNA（doublet 通常：高 nCount + 高 nFeature）
# df_qc <- combined@meta.data
# p2 <- ggplot(df_qc, aes(x = nCount_RNA, y = percent.mt)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + labs(x = "nCount_RNA (UMI)", y = "percent.mt") + theme_bw()
# p3 <- ggplot(df_qc, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + labs(x = "nCount_RNA (UMI)", y = "nFeature_RNA (genes)") + theme_bw()
# p_qc2 <- p2 / p3
# ggsave(filename = "temp/QC_scatter_preQC.pdf",plot = p_qc2, width = 14, height = 10, dpi = 300)
# # ---- 2) Print pre-QC cell counts per group ----让你知道过滤前每组规模，便于判断过滤是否过度。
# cat("==== Pre-QC cell counts ====\n")
# print(table(combined$group))
# cat("Total cells:", ncol(combined), "\n\n")
# # 计算“percent.mt > 60% 的细胞占比”，快速感知你的数据里极端坏细胞多不多
# extreme_mt <- mean(combined$percent.mt > 60) * 100
# cat(sprintf("Diagnostic: percent.mt > 60%% cells = %.2f%%\n\n", extreme_mt))
# # ---- 3) Define medium thresholds ----设置“中等阈值”
# min_features <- 500
# max_features <- 6500
# min_counts   <- 1500   # 如果你发现过滤后细胞数太少，最常调的是 min_counts（1000→800）
# max_counts   <- 30000
# max_mt       <- 15     # 如果你发现过滤后细胞数太少，最常调的是 max_mt（30→35/40）。
# # ---- 4) HB/RB extreme cutoff (per-group 99th percentile) ----hb/rb 不要设死值，而是“每组只砍掉最极端的 1%”
# # This avoids one group dominating the global quantile.
# md <- combined@meta.data %>%
#   mutate(cell = rownames(.))
# hb_cut_by_group <- md %>%  #每组情况可能不同，用“每组各自 top1%”更公平，不会出现某组被误砍更多。
#   group_by(group) %>%
#   summarise(hb_q99 = quantile(percent.hb, probs = 0.99, na.rm = TRUE), rb_q99 = quantile(percent.rb, probs = 0.99, na.rm = TRUE), .groups = "drop")
# cat("==== Per-group 99% quantiles for HB/RB ====\n") #把上面每组阈值打印出来，让你知道 hb/rb “砍掉极端值”到底砍到哪面
# print(hb_cut_by_group)
# cat("\n")
# # Join thresholds back to metadata 这样每个细胞都知道“自己所在组的阈值是多少”，方便下一步筛选
# md2 <- md %>%
#   left_join(hb_cut_by_group, by = "group")
# cells_keep <- md2 %>%  #真正执行 QC 过滤：选出要保留的细胞3v
#   dplyr::filter(nFeature_RNA >= min_features, nFeature_RNA <= max_features, nCount_RNA >= min_counts, nCount_RNA   <= max_counts, percent.mt   <= max_mt, percent.hb   <= hb_q99, percent.rb   <= rb_q99) %>%
#   dplyr::pull(cell)  #把保留细胞的名字抽出来做成一个向量
# cat("Cells kept after QC:", length(cells_keep), " / ", ncol(combined), "\n\n") #看每组保留率是否合理、是否某组被砍得特别严重
# combined_qc <- subset(combined, cells = cells_keep)
# # ---- 5) Print post-QC cell counts per group ----
# cat("==== Post-QC cell counts ====\n")
# print(table(combined_qc$group))
# cat("Total cells:", ncol(combined_qc), "\n\n")
# # ---- 6) Post-QC plots (violin + scatter) ----
# # (A) Violin plots: set raster=FALSE to avoid rasterization-related internal errors
# p_vln_all <- VlnPlot(combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "percent.rb"), pt.size = 0.1, raster = FALSE, ncol = 5)
# p_vln_grp <- VlnPlot(combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "percent.rb"), group.by = "group", pt.size = 0.1, raster = FALSE, ncol = 5)
# p_vln <- p_vln_all / p_vln_grp
# ggsave(filename = "temp/QC_violin_postQC_overall_and_group.pdf", plot = p_vln, width = 14, height = 10)
# # (B) Scatter plots using ggplot facet (stable for large data)
# df_qc <- combined_qc@meta.data
# p_sc1 <- ggplot(df_qc, aes(x = nCount_RNA, y = percent.mt)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + scale_x_log10() + geom_hline(yintercept = max_mt, linetype = 2) + labs(x = "nCount_RNA (UMI, log10)", y = "percent.mt") + theme_bw()
# p_sc2 <- ggplot(df_qc, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + scale_x_log10() + geom_hline(yintercept = c(min_features, max_features), linetype = 2) + geom_vline(xintercept = c(min_counts, max_counts), linetype = 2) + labs(x = "nCount_RNA (UMI, log10)", y = "nFeature_RNA (genes)") + theme_bw()
# p_sc <- p_sc1/p_sc2
# ggsave("temp/QC_scatter_postQC.pdf", plot = p_sc, width = 12, height = 6)
# cat("Post-QC plots saved to:\n", file.path(getwd(), "temp/QC_violin_postQC_overall_and_group.pdf"), "\n", file.path(getwd(), "temp/QC_scatter_postQC_nCount_vs_percentMT.pdf"), "\n", file.path(getwd(), "temp/QC_scatter_postQC_nCount_vs_nFeature.pdf"), "\n")
# # 细胞过滤：保留 nFeature_RNA 在 [200,7500] 内、线粒体比例小于 25%
# # cat("细胞过滤前：", ncol(combined), "个细胞\n")
# # combined <- subset(combined, subset = nFeature_RNA >= 200 &
# #   nFeature_RNA <= 7500 &
# #   percent.mt < 25)
# # cat("细胞过滤后：", ncol(combined), "个细胞\n")
# cat("2. 完成【数据过滤】\n")
#
# ## --------------------- 2. 数据过滤 ---------------------
# #cat("2. 开始【数据过滤】\n")
# ## 基因过滤：保留表达大于0的细胞数大于3的基因
# #counts_matrix <- LayerData(combined, assay = "RNA", layer = "counts")
# #genes_cell_counts <- rowSums(counts_matrix > 0)
# #genes_to_keep <- names(genes_cell_counts[genes_cell_counts > 3])
# #cat("基因过滤前：", nrow(combined), "个基因\n")
# #combined <- subset(combined, features = genes_to_keep)
# #cat("基因过滤后：", nrow(combined), "个基因\n")
# ## 计算线粒体基因比例
# #combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
# ## 细胞过滤：保留 nFeature_RNA 在 [200,7500] 内、线粒体比例小于 25%
# #cat("细胞过滤前：", ncol(combined), "个细胞\n")
# #combined <- subset(combined, subset = nFeature_RNA >= 200 &
# #  nFeature_RNA <= 7500 &
# #  percent.mt < 25)
# #cat("细胞过滤后：", ncol(combined), "个细胞\n")
# #cat("2. 完成【数据过滤】\n")
#
#
# # --------------------- 3. 双细胞检测 ---------------------
# cat("3. 开始【双细胞检测及过滤】\n")
# # 预处理数据以便进行双细胞检测
# combined <- NormalizeData(combined)
# combined <- FindVariableFeatures(combined, nfeatures = 2000)
# # 仅缩放变量基因以节省内存
# combined <- ScaleData(combined, features = VariableFeatures(combined))
# combined <- RunPCA(combined, npcs = 30)
# # 使用 scDblFinder 进行双细胞检测
# counts <- LayerData(combined, assay = "RNA", layer = "counts")
# sce <- SingleCellExperiment(assays = list(counts = counts))
# sce <- scDblFinder(sce, dbr = 0.06)  # 设置预期双细胞比例为6%
# # 将双细胞信息添加回 Seurat 对象
# combined$scDblFinder.class <- sce$scDblFinder.class
# # 过滤掉被标记为双细胞的细胞
# combined <- subset(combined, subset = scDblFinder.class == "singlet")
# cat("3. 完成【双细胞检测及过滤】。剩余细胞数:", ncol(combined), "\n")
#
#
# # --------------------- 4. 数据预处理 (归一化) ---------------------
# cat("4. 开始【数据预处理】\n")
# # 标准归一化和对数转换
# combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 1e4)
# # 找出高变异基因
# combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
# # 优化内存使用的缩放数据
# # 1. 仅缩放变量基因而非所有基因
# # 2. 分批处理缩放数据
# var_genes <- VariableFeatures(combined)
# cat("正在缩放数据（仅变量基因），节省内存使用...\n")
# combined <- ScaleData(combined, features = var_genes, vars.to.regress = "percent.mt", clip.max = 10)
# cat("4. 完成【数据预处理】\n")
# saveRDS(combined, file = "./temp/combined_normalization_4.rds")
# combined <- readRDS("./temp/combined_normalization_4.rds")
#
#
# # --------------------- 5. 降维 & 聚类 ---------------------
# cat("5. 开始【降维与聚类分析】\n")
# # 使用PCA结果构建KNN图
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:15, k.param = 20)
# # Louvain社区聚类，resolution控制聚类粒度
# combined <- FindClusters(combined, resolution = 0.5)
# # 计算UMAP（非线性降维可视化）
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:15)
# # 绘制UMAP聚类结果
# p <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
#   ggtitle("UMAP 聚类结果")
# ggsave("./temp/UMAP_louvain.pdf", plot = p, width = 10, height = 8)
# cat("5. 完成【降维与聚类分析】\n")
# saveRDS(combined, file = "./temp/combined_umap_5.rds")
# combined <- readRDS("./temp/combined_umap_5.rds")
#
#
# # --------------------- 6. Marker 基因分析（差异表达） ---------------------
# cat("6. 开始【Marker 基因分析】\n")
# # 确保使用的是正确的标识符
# Idents(combined) <- combined$seurat_clusters
# # 使用 FindMarkers 逐个群体查找
# all_markers <- list()
# clusters <- levels(Idents(combined))
# combined <- JoinLayers(object = combined, assay = "RNA", layers = c("counts", "data"))
# for (cluster in clusters) {
#   markers <- FindMarkers(combined, ident.1 = cluster, ident.2 = NULL, only.pos = TRUE,
#                          min.pct = 0.25, logfc.threshold = 0.25)
#   if (nrow(markers) > 0) {
#     markers$cluster <- cluster
#     all_markers[[cluster]] <- markers
#   }
# }
# # 合并所有结果
# if (length(all_markers) > 0) {
#   all_markers_df <- do.call(rbind, all_markers)
# }
# # 转换为data.table以备后续使用
# markers_dt <- if (length(all_markers) > 0) data.table::rbindlist(all_markers, idcol = "cluster") else data.table::data.table()
# cat("6. 完成【Marker 基因分析】\n")
# saveRDS(all_markers, file = "./temp/marker_analysis_6.rds")
# all_markers <- readRDS("./temp/marker_analysis_6.rds")
# # 保存 combined（含当前 Idents、layers 等）
# saveRDS(combined, file = "./temp/combined_after_marker_6.rds")
# cat("combined 已保存到：", normalizePath("./temp/combined_after_marker_6.rds"), "\n")

combined <- readRDS("./temp/combined_after_marker_6.rds")
all_markers <- readRDS("./temp/marker_analysis_6.rds")

# --------------------- 7.1 提取 T cells 并重聚类 --------------------
# --------------------- 7.1 提取 T cells 并重聚类 ---------------------
if (file.exists("./temp/t_cells_reclustered_7.1.rds")) {
  t_cells <- readRDS("./temp/t_cells_reclustered_7.1.rds")
  message("Loaded ./temp/t_cells_reclustered_7.1.rds")
} else {
  t_cells <- subset(combined, subset = Trac > 0 | Cd3e > 0 | Cd3d > 0)
  t_cells$group <- factor(t_cells$group)

  DefaultAssay(t_cells) <- "RNA"
  t_cells <- NormalizeData(t_cells)
  t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
  t_cells <- ScaleData(t_cells, features = VariableFeatures(t_cells), vars.to.regress = "percent.mt")
  t_cells <- RunPCA(t_cells, npcs = 30)
  t_cells <- FindNeighbors(t_cells, reduction = "pca", dims = 1:20, k.param = 20)
  t_cells <- FindClusters(t_cells, resolution = 0.6)
  t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:20)

  saveRDS(t_cells, file = "./temp/t_cells_reclustered_7.1.rds")
  message("Saved ./temp/t_cells_reclustered_7.1.rds")
}

# =========================================================
# 8. CD8 Tex analysis - full coherent workflow
# =========================================================
dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
dir.create("./figs", showWarnings = FALSE, recursive = TRUE)
# =========================================================
# DotPlot 全局颜色覆盖（只影响基因表达量 DotPlot）
# =========================================================

# 定义你的渐变色（低 → 中 → 高）
DOT_GRADIENT_COLS <- c("#377EB8", "#A6A6A6", "#F46D09")

# 不覆盖 Seurat 原函数，单独定义一个带固定配色的 DotPlot 包装器
dotplot3c <- function(...) {
  p <- Seurat::DotPlot(...)
  p + ggplot2::scale_colour_gradientn(colors = DOT_GRADIENT_COLS)
}

# 统一 group 命名，避免 s1-hIL / s1.hIL / s1_hIL 混乱
standardize_group_labels <- function(x) {
  x <- as.character(x)
  x <- gsub("-", "_", x)
  x <- gsub("\\.", "_", x)
  x <- trimws(x)

  x[x %in% c("control", "ctrl", "control_1", "control1")] <- "control1"
  x[x %in% c("s1_hIL", "s1_hil", "s1_hIL_1")] <- "s1_hIL"
  x[x %in% c("s40kD", "s40kd", "s40_kD", "s40_kd")] <- "s40kD"

  factor(x, levels = c("control1", "s1_hIL", "s40kD"))
}

# 强制检查 group，并补一个 sample_id 接口；如果没有真实 sample_id，只给一个临时占位
ensure_group_and_sample <- function(obj, obj_name = deparse(substitute(obj))) {
  if (!("group" %in% colnames(obj@meta.data))) {
    stop(obj_name, "@meta.data 缺少 group 列")
  }

  obj$group <- standardize_group_labels(obj$group)

  if (!("sample_id" %in% colnames(obj@meta.data))) {
    if ("orig.ident" %in% colnames(obj@meta.data)) {
      obj$sample_id <- as.character(obj$orig.ident)
    } else {
      obj$sample_id <- paste0(as.character(obj$group), "_sample1")
      message("[WARN] ", obj_name, " 没有 sample_id/orig.ident；已临时创建 sample_id。")
      message("[WARN] 这说明当前脚本还不满足严格 pseudobulk DE 的发表级条件。")
    }
  }

  obj$sample_id <- as.character(obj$sample_id)
  return(obj)
}

# 检查映射表是否完整，防止 composition/function 映射悄悄出错
assert_mapping_complete <- function(values, mapper, mapper_name = "mapper") {
  miss_in_mapper <- base::setdiff(unique(as.character(values)), names(mapper))
  if (length(miss_in_mapper) > 0) {
    stop(sprintf("%s 缺少以下键：%s", mapper_name, paste(miss_in_mapper, collapse = ", ")))
  }
}
# -------------------------
# 8.1 工具函数
# -------------------------
intersect2 <- function(a, b) base::intersect(a, b)

pick_gene <- function(obj, candidates) {
  hit <- candidates[candidates %in% rownames(obj)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

present_only <- function(obj, genes) {
  unique(genes[genes %in% rownames(obj)])
}

score_by_genesum <- function(obj, genes) {
  genes <- present_only(obj, genes)
  if (length(genes) == 0) {
    out <- rep(0, ncol(obj))
    names(out) <- colnames(obj)
    return(out)
  }
  df <- FetchData(obj, vars = genes)
  s <- rowSums(df, na.rm = TRUE)
  names(s) <- rownames(df)
  return(s)
}

bootstrap_prop_subsample <- function(group_cells_list, cd8_cell_set, frac = 0.8, B = 500, seed = 123) {
  set.seed(seed)
  groups <- names(group_cells_list)
  n_min <- min(sapply(group_cells_list, length))
  sample_size <- max(10, floor(n_min * frac))

  res <- lapply(groups, function(g) {
    cells <- group_cells_list[[g]]
    props <- replicate(B, {
      samp <- sample(cells, size = sample_size, replace = FALSE)
      mean(samp %in% cd8_cell_set)
    })
    data.frame(group = g, prop = props, sample_size = sample_size)
  })
  do.call(rbind, res)
}

recluster_and_export <- function(
  obj,
  prefix,
  group_var = "group",
  dims_use = 1:20,
  resolution_use = 0.6,
  nfeatures_use = 2000
) {
  obj <- ensure_group_and_sample(obj, prefix)
  DefaultAssay(obj) <- "RNA"

  obj <- FindVariableFeatures(obj, nfeatures = nfeatures_use)
  saveRDS(obj, paste0("./temp/", prefix, "_after_HVG.rds"))

  obj <- ScaleData(obj, verbose = FALSE)
  saveRDS(obj, paste0("./temp/", prefix, "_after_ScaleData.rds"))

  obj <- RunPCA(obj, npcs = max(dims_use), verbose = FALSE)
  saveRDS(obj, paste0("./temp/", prefix, "_after_PCA.rds"))

  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  saveRDS(obj, paste0("./temp/", prefix, "_after_Neighbors.rds"))

  obj <- FindClusters(obj, resolution = resolution_use, verbose = FALSE)
  saveRDS(obj, paste0("./temp/", prefix, "_after_Clusters.rds"))

  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  saveRDS(obj, paste0("./temp/", prefix, ".rds"))

  Idents(obj) <- "seurat_clusters"

  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  write.csv(markers, paste0("./temp/", prefix, "_markers.csv"), row.names = FALSE)

  avg_expr <- AverageExpression(
    obj,
    group.by = "seurat_clusters",
    assays = "RNA",
    slot = "data"
  )$RNA
  write.csv(avg_expr, paste0("./temp/", prefix, "_avgexpr.csv"))

  comp <- as.data.frame(table(
    group = obj[[group_var, drop = TRUE]],
    cluster = obj$seurat_clusters
  )) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = Freq / sum(Freq)) %>%
    dplyr::ungroup()
  write.csv(comp, paste0("./temp/", prefix, "_composition.csv"), row.names = FALSE)

  return(obj)
}
# -------------------------
# 8.1b QC / identity helper functions
# -------------------------
get_doublet_score_col <- function(obj) {
  cand <- c("scDblFinder.score", "doublet_score", "DoubletScore", "doubletScore")
  hit <- cand[cand %in% colnames(obj@meta.data)]
  if (length(hit) > 0) return(hit[1])
  return(NA_character_)
}

safe_meta_feature_panel <- function(obj, include_doublet = TRUE) {
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  dbl_col <- get_doublet_score_col(obj)
  if (include_doublet && !is.na(dbl_col)) feats <- c(feats, dbl_col)
  feats[feats %in% colnames(obj@meta.data)]
}

make_qc_violin <- function(obj, group_by = NULL, title = NULL, pt_size = 0.05) {
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  dbl_col <- get_doublet_score_col(obj)
  if (!is.na(dbl_col)) feats <- c(feats, dbl_col)
  feats <- feats[feats %in% colnames(obj@meta.data)]
  if (length(feats) == 0) return(NULL)

  md <- obj@meta.data
  md$cell <- rownames(md)

  if (is.null(group_by) || !(group_by %in% colnames(md))) {
    md$plot_group <- "All cells"
    group_col <- "plot_group"
  } else {
    md$plot_group <- as.character(md[[group_by]])
    group_col <- "plot_group"
  }

  keep_cols <- unique(c(group_col, feats))
  df <- md[, keep_cols, drop = FALSE]

  df_long <- tidyr::pivot_longer(
    data = df,
    cols = dplyr::all_of(feats),
    names_to = "feature",
    values_to = "value"
  )
  df_long <- df_long[!is.na(df_long$value), , drop = FALSE]

  p <- ggplot(df_long, aes(x = .data[[group_col]], y = value)) +
    geom_violin(fill = "grey85", color = "black", scale = "width", trim = TRUE) +
    geom_jitter(width = 0.15, size = pt_size, alpha = 0.2) +
    facet_wrap(~ feature, scales = "free_y", ncol = min(4, length(feats))) +
    theme_bw() +
    labs(x = "", y = "Value")

  if (length(unique(df_long[[group_col]])) > 1) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  if (!is.null(title)) p <- p + ggtitle(title)
  return(p)
}

make_qc_scatter_pair <- function(obj, group_by = "group", title_prefix = "QC") {
  md <- obj@meta.data
  md$cell <- rownames(md)

  p1 <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point(size = 0.25, alpha = 0.5) +
    facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
    scale_x_log10() +
    labs(x = "nCount_RNA (log10)", y = "nFeature_RNA", title = paste0(title_prefix, ": nCount vs nFeature")) +
    theme_bw()

  p2 <- ggplot(md, aes(x = nCount_RNA, y = percent.mt)) +
    geom_point(size = 0.25, alpha = 0.5) +
    facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
    scale_x_log10() +
    labs(x = "nCount_RNA (log10)", y = "percent.mt", title = paste0(title_prefix, ": nCount vs percent.mt")) +
    theme_bw()

  dbl_col <- get_doublet_score_col(obj)
  if (!is.na(dbl_col)) {
    p3 <- ggplot(md, aes(x = nCount_RNA, y = .data[[dbl_col]])) +
      geom_point(size = 0.25, alpha = 0.5) +
      facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
      scale_x_log10() +
      labs(x = "nCount_RNA (log10)", y = dbl_col, title = paste0(title_prefix, ": nCount vs doublet score")) +
      theme_bw()
    return(list(p1 = p1, p2 = p2, p3 = p3))
  } else {
    return(list(p1 = p1, p2 = p2))
  }
}

save_qc_bundle <- function(obj, prefix, group_by = "group") {
  dir.create("./figs", showWarnings = FALSE, recursive = TRUE)
  dir.create("./temp", showWarnings = FALSE, recursive = TRUE)

  md_export_cols <- unique(c(
    "group", "sample_id", "orig.ident",
    "nFeature_RNA", "nCount_RNA", "percent.mt",
    "scDblFinder.class", get_doublet_score_col(obj)
  ))
  md_export_cols <- md_export_cols[md_export_cols %in% colnames(obj@meta.data)]
  md_export <- obj@meta.data[, md_export_cols, drop = FALSE]
  write.csv(md_export, paste0("./temp/", prefix, "_metadata_qc_table.csv"))

  p_vln_all <- tryCatch(
    make_qc_violin(obj, group_by = NULL, title = paste0(prefix, ": overall QC")),
    error = function(e) NULL
  )
  if (!is.null(p_vln_all)) {
    ggsave(paste0("./figs/", prefix, "_QC_violin_overall.pdf"), p_vln_all, width = 12, height = 5)
  }

  p_vln_group <- tryCatch(
    make_qc_violin(obj, group_by = group_by, title = paste0(prefix, ": QC by ", group_by)),
    error = function(e) NULL
  )
  if (!is.null(p_vln_group)) {
    ggsave(paste0("./figs/", prefix, "_QC_violin_by_", group_by, ".pdf"), p_vln_group, width = 12, height = 5)
  }

  md <- obj@meta.data
  if (all(c("nCount_RNA", "nFeature_RNA", group_by) %in% colnames(md))) {
    p1 <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(size = 0.25, alpha = 0.5) +
      facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
      scale_x_log10() +
      labs(x = "nCount_RNA (log10)", y = "nFeature_RNA", title = paste0(prefix, ": nCount vs nFeature")) +
      theme_bw()
    ggsave(paste0("./figs/", prefix, "_p1.pdf"), p1, width = 10, height = 6)
  }

  if (all(c("nCount_RNA", "percent.mt", group_by) %in% colnames(md))) {
    p2 <- ggplot(md, aes(x = nCount_RNA, y = percent.mt)) +
      geom_point(size = 0.25, alpha = 0.5) +
      facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
      scale_x_log10() +
      labs(x = "nCount_RNA (log10)", y = "percent.mt", title = paste0(prefix, ": nCount vs percent.mt")) +
      theme_bw()
    ggsave(paste0("./figs/", prefix, "_p2.pdf"), p2, width = 10, height = 6)
  }

  dbl_col <- get_doublet_score_col(obj)
  if (!is.na(dbl_col) && all(c("nCount_RNA", group_by, dbl_col) %in% colnames(md))) {
    p3 <- ggplot(md, aes(x = nCount_RNA, y = .data[[dbl_col]])) +
      geom_point(size = 0.25, alpha = 0.5) +
      facet_wrap(stats::as.formula(paste("~", group_by)), ncol = 3, scales = "free") +
      scale_x_log10() +
      labs(x = "nCount_RNA (log10)", y = dbl_col, title = paste0(prefix, ": nCount vs doublet score")) +
      theme_bw()
    ggsave(paste0("./figs/", prefix, "_p3.pdf"), p3, width = 10, height = 6)
  }

  if ("scDblFinder.class" %in% colnames(md) && group_by %in% colnames(md)) {
    df_cls <- as.data.frame(table(group = md[[group_by]], dbl_class = md$scDblFinder.class))
    p_cls <- ggplot(df_cls, aes(x = group, y = Freq, fill = dbl_class)) +
      geom_col(position = "fill") +
      theme_bw() +
      labs(title = paste0(prefix, ": scDblFinder.class proportion"), x = group_by, y = "Fraction", fill = "Doublet class")
    ggsave(paste0("./figs/", prefix, "_doublet_class_fraction.pdf"), p_cls, width = 8, height = 5)
  }

  invisible(TRUE)
}

plot_pre_post_qc_compare <- function(obj_pre, obj_post, prefix = "CD8_pre_post_QC") {
  dir.create("./figs", showWarnings = FALSE, recursive = TRUE)
  dir.create("./temp", showWarnings = FALSE, recursive = TRUE)

  md_pre <- obj_pre@meta.data
  md_post <- obj_post@meta.data

  md_pre$qc_stage <- "pre_filter"
  md_post$qc_stage <- "post_filter"

  common_cols <- base::intersect(colnames(md_pre), colnames(md_post))

  feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  dbl_pre <- get_doublet_score_col(obj_pre)
  dbl_post <- get_doublet_score_col(obj_post)
  if (!is.na(dbl_pre) && !is.na(dbl_post) && dbl_pre == dbl_post) feats <- c(feats, dbl_pre)
  feats <- feats[feats %in% common_cols]

  keep_cols <- base::unique(c("group", "qc_stage", "scDblFinder.class", feats))
  keep_cols <- keep_cols[keep_cols %in% common_cols]

  md_pre2 <- md_pre[, keep_cols, drop = FALSE]
  md_post2 <- md_post[, keep_cols, drop = FALSE]
  md_all <- rbind(md_pre2, md_post2)

  write.csv(md_all, paste0("./temp/", prefix, "_metadata_combined.csv"), row.names = TRUE)

  if (length(feats) > 0) {
    df_long_stage <- tidyr::pivot_longer(
      data = md_all,
      cols = dplyr::all_of(feats),
      names_to = "feature",
      values_to = "value"
    )
    df_long_stage <- df_long_stage[!is.na(df_long_stage$value), , drop = FALSE]

    p_stage <- ggplot(df_long_stage, aes(x = qc_stage, y = value)) +
      geom_violin(fill = "grey85", color = "black", scale = "width", trim = TRUE) +
      geom_jitter(width = 0.15, size = 0.03, alpha = 0.2) +
      facet_wrap(~ feature, scales = "free_y", ncol = min(4, length(feats))) +
      theme_bw() +
      labs(title = prefix, x = "", y = "Value")
    ggsave(paste0("./figs/", prefix, "_violin_by_stage.pdf"), p_stage, width = 12, height = 5)
  }

  if (length(feats) > 0 && "group" %in% colnames(md_all)) {
    df_long_group <- tidyr::pivot_longer(
      data = md_all,
      cols = dplyr::all_of(feats),
      names_to = "feature",
      values_to = "value"
    )
    df_long_group <- df_long_group[!is.na(df_long_group$value), , drop = FALSE]

    p_stage_group <- ggplot(df_long_group, aes(x = qc_stage, y = value, fill = qc_stage)) +
      geom_violin(scale = "width", trim = TRUE, color = "black") +
      facet_grid(feature ~ group, scales = "free_y") +
      theme_bw() +
      labs(title = paste0(prefix, ": grouped QC stage comparison"), x = "", y = "Value") +
      theme(legend.position = "none")
    ggsave(paste0("./figs/", prefix, "_violin_by_stage_and_group.pdf"), p_stage_group, width = 12, height = 8)
  }

  if (all(c("nCount_RNA", "nFeature_RNA", "group", "qc_stage") %in% colnames(md_all))) {
    p1 <- ggplot(md_all, aes(x = nCount_RNA, y = nFeature_RNA, color = qc_stage)) +
      geom_point(size = 0.2, alpha = 0.5) +
      facet_wrap(~ group, ncol = 3, scales = "free") +
      scale_x_log10() +
      theme_bw() +
      labs(title = paste0(prefix, ": nCount vs nFeature"), x = "nCount_RNA (log10)", y = "nFeature_RNA")
    ggsave(paste0("./figs/", prefix, "_scatter_nCount_nFeature.pdf"), p1, width = 10, height = 6)
  }

  if (all(c("nCount_RNA", "percent.mt", "group", "qc_stage") %in% colnames(md_all))) {
    p2 <- ggplot(md_all, aes(x = nCount_RNA, y = percent.mt, color = qc_stage)) +
      geom_point(size = 0.2, alpha = 0.5) +
      facet_wrap(~ group, ncol = 3, scales = "free") +
      scale_x_log10() +
      theme_bw() +
      labs(title = paste0(prefix, ": nCount vs percent.mt"), x = "nCount_RNA (log10)", y = "percent.mt")
    ggsave(paste0("./figs/", prefix, "_scatter_nCount_percentMT.pdf"), p2, width = 10, height = 6)
  }

  dbl_col <- if (!is.na(dbl_pre) && !is.na(dbl_post) && dbl_pre == dbl_post) dbl_pre else NA_character_
  if (!is.na(dbl_col) && all(c("nCount_RNA", "group", "qc_stage", dbl_col) %in% colnames(md_all))) {
    p3 <- ggplot(md_all, aes(x = nCount_RNA, y = .data[[dbl_col]], color = qc_stage)) +
      geom_point(size = 0.2, alpha = 0.5) +
      facet_wrap(~ group, ncol = 3, scales = "free") +
      scale_x_log10() +
      theme_bw() +
      labs(title = paste0(prefix, ": nCount vs doublet score"), x = "nCount_RNA (log10)", y = dbl_col)
    ggsave(paste0("./figs/", prefix, "_scatter_nCount_doubletScore.pdf"), p3, width = 10, height = 6)
  }

  if ("scDblFinder.class" %in% colnames(md_all) && all(c("group", "qc_stage") %in% colnames(md_all))) {
    df_cls <- as.data.frame(table(group = md_all$group, qc_stage = md_all$qc_stage, dbl_class = md_all$scDblFinder.class))
    p_cls <- ggplot(df_cls, aes(x = qc_stage, y = Freq, fill = dbl_class)) +
      geom_col(position = "fill") +
      facet_wrap(~ group, ncol = 3) +
      theme_bw() +
      labs(title = paste0(prefix, ": scDblFinder.class fraction"), x = "", y = "Fraction", fill = "Doublet class")
    ggsave(paste0("./figs/", prefix, "_doublet_class_fraction.pdf"), p_cls, width = 10, height = 6)
  }

  invisible(TRUE)
}

# -------------------------
# 8.2 基础检查
# -------------------------
combined <- ensure_group_and_sample(combined, "combined")
t_cells  <- ensure_group_and_sample(t_cells,  "t_cells")

DefaultAssay(combined) <- "RNA"
DefaultAssay(t_cells) <- "RNA"

# -------------------------
# 8.3 构建 immune_cells
# -------------------------
g_cd45 <- if ("Ptprc" %in% rownames(combined)) "Ptprc" else if ("PTPRC" %in% rownames(combined)) "PTPRC" else NA

if (!is.na(g_cd45)) {
  cd45_df <- FetchData(combined, vars = g_cd45)
  keep_cells_immune <- rownames(cd45_df)[cd45_df[[g_cd45]] > 0]
  immune_cells <- subset(combined, cells = keep_cells_immune)
} else {
  stop("combined 中找不到 Ptprc/PTPRC，无法自动构建 immune_cells")
}

immune_cells$group <- factor(immune_cells$group)
saveRDS(immune_cells, "./temp/immune_cells.rds")

# -------------------------
# 8.3b QC overview for immune_cells
# -------------------------
immune_cells <- ensure_group_and_sample(immune_cells, "immune_cells")
save_qc_bundle(immune_cells, prefix = "Immune_cells_QC", group_by = "group")

# -------------------------
# 8.4 提取并过滤 CD8
# -------------------------
g_trac <- pick_gene(t_cells, c("TRAC", "Trac"))
g_cd3d <- pick_gene(t_cells, c("CD3D", "Cd3d"))
g_cd3e <- pick_gene(t_cells, c("CD3E", "Cd3e"))
g_cd8a <- pick_gene(t_cells, c("CD8A", "Cd8a"))
g_cd8b <- pick_gene(t_cells, c("CD8B", "Cd8b", "CD8B1", "Cd8b1"))
g_cd4  <- pick_gene(t_cells, c("CD4", "Cd4"))

panel_myeloid <- present_only(t_cells, c("Lyz2","Lyz","Itgam","Itgax","Tyrobp","Lst1","S100a8","S100a9"))
panel_nk <- present_only(t_cells, c("Ncr1","Klrd1","Klrk1","Klra1","Klra8","Trdc","Trgc1","Trgc2","Fcgr3","Tyrobp"))
panel_treg <- present_only(t_cells, c("Foxp3","Il2ra","Ikzf2","Lrrc32","Tnfrsf4","Tnfrsf18"))
panel_tumor <- present_only(t_cells, c("Mlana","Slc45a2","Pmel","Tyr","Dct","Mc1r"))

vars_gate <- unique(na.omit(c(g_trac, g_cd3d, g_cd3e, g_cd8a, g_cd8b, g_cd4)))
df_gate <- FetchData(t_cells, vars = vars_gate)

t_score <- rep(0, nrow(df_gate))
if (!is.na(g_trac)) t_score <- t_score + df_gate[[g_trac]]
if (!is.na(g_cd3d)) t_score <- t_score + df_gate[[g_cd3d]]
if (!is.na(g_cd3e)) t_score <- t_score + df_gate[[g_cd3e]]

cd8_score <- rep(0, nrow(df_gate))
if (!is.na(g_cd8a)) cd8_score <- cd8_score + df_gate[[g_cd8a]]
if (!is.na(g_cd8b)) cd8_score <- cd8_score + df_gate[[g_cd8b]]

cd4_expr <- if (!is.na(g_cd4)) df_gate[[g_cd4]] else rep(0, nrow(df_gate))
use_cd4_exclude <- !is.na(g_cd4) && (sum(cd4_expr > 0) > 0)

is_cd8_raw <- (t_score > 0) & (cd8_score > 0) & (if (use_cd4_exclude) cd4_expr == 0 else TRUE)

cd8_cells_raw <- subset(t_cells, cells = rownames(df_gate)[is_cd8_raw])
cd8_cells_raw$group <- droplevels(cd8_cells_raw$group)
saveRDS(cd8_cells_raw, "./temp/cd8_cells_gated_raw.rds")

cd8_cells_raw$score_myeloid <- score_by_genesum(cd8_cells_raw, panel_myeloid)
cd8_cells_raw$score_nk      <- score_by_genesum(cd8_cells_raw, panel_nk)
cd8_cells_raw$score_treg    <- score_by_genesum(cd8_cells_raw, panel_treg)
cd8_cells_raw$score_tumor   <- score_by_genesum(cd8_cells_raw, panel_tumor)

has_ncount <- "nCount_RNA" %in% colnames(cd8_cells_raw@meta.data)
has_nfeat  <- "nFeature_RNA" %in% colnames(cd8_cells_raw@meta.data)

thr_my <- quantile(cd8_cells_raw$score_myeloid, probs = 0.98, na.rm = TRUE)
thr_nk <- quantile(cd8_cells_raw$score_nk, probs = 0.99, na.rm = TRUE)
thr_tr <- quantile(cd8_cells_raw$score_treg, probs = 0.995, na.rm = TRUE)
thr_tu <- quantile(cd8_cells_raw$score_tumor, probs = 0.995, na.rm = TRUE)

keep_flag <- rep(TRUE, ncol(cd8_cells_raw))
names(keep_flag) <- colnames(cd8_cells_raw)

if (sd(cd8_cells_raw$score_myeloid, na.rm = TRUE) > 0) keep_flag <- keep_flag & (cd8_cells_raw$score_myeloid <= thr_my)
if (sd(cd8_cells_raw$score_nk, na.rm = TRUE) > 0)      keep_flag <- keep_flag & (cd8_cells_raw$score_nk <= thr_nk)
if (sd(cd8_cells_raw$score_treg, na.rm = TRUE) > 0)    keep_flag <- keep_flag & (cd8_cells_raw$score_treg <= thr_tr)
if (sd(cd8_cells_raw$score_tumor, na.rm = TRUE) > 0)   keep_flag <- keep_flag & (cd8_cells_raw$score_tumor <= thr_tu)

if (has_ncount) {
  thr_nc <- quantile(cd8_cells_raw$nCount_RNA, probs = 0.99, na.rm = TRUE)
  keep_flag <- keep_flag & (cd8_cells_raw$nCount_RNA <= thr_nc)
}
if (has_nfeat) {
  thr_nf <- quantile(cd8_cells_raw$nFeature_RNA, probs = 0.99, na.rm = TRUE)
  keep_flag <- keep_flag & (cd8_cells_raw$nFeature_RNA <= thr_nf)
}

keep_cells <- names(keep_flag)[keep_flag]
keep_cells <- intersect2(keep_cells, colnames(cd8_cells_raw))

filter_stat <- data.frame(
  stage = c("raw_gate", "after_filter"),
  n_cells = c(ncol(cd8_cells_raw), length(keep_cells))
)
write.csv(filter_stat, "./temp/CD8_filtering_cell_numbers.csv", row.names = FALSE)

by_group_raw <- as.data.frame(table(group = cd8_cells_raw$group))
colnames(by_group_raw) <- c("group", "n_raw")

cd8_cells_raw@graphs <- list()
cd8_cells_raw@neighbors <- list()

cd8_cells <- subset(cd8_cells_raw, cells = keep_cells)
cd8_cells$group <- droplevels(cd8_cells$group)

by_group_flt <- as.data.frame(table(group = cd8_cells$group))
colnames(by_group_flt) <- c("group", "n_filtered")

by_group_filter <- dplyr::left_join(by_group_raw, by_group_flt, by = "group") %>%
  dplyr::mutate(
    n_filtered = ifelse(is.na(n_filtered), 0, n_filtered),
    removed = n_raw - n_filtered,
    removed_rate = removed / n_raw
  )
write.csv(by_group_filter, "./temp/CD8_filtering_by_group.csv", row.names = FALSE)

saveRDS(cd8_cells, "./temp/cd8_cells_gated_filtered.rds")

# -------------------------
# 8.4b QC overview for CD8 raw / filtered objects
# -------------------------
cd8_cells_raw <- ensure_group_and_sample(cd8_cells_raw, "cd8_cells_raw")
cd8_cells     <- ensure_group_and_sample(cd8_cells, "cd8_cells")

# 导出 raw 和 filtered 的 QC 图
save_qc_bundle(cd8_cells_raw, prefix = "CD8_raw_gated", group_by = "group")
save_qc_bundle(cd8_cells,     prefix = "CD8_filtered_final", group_by = "group")

# 导出 pre/post 对照图
plot_pre_post_qc_compare(
  obj_pre = cd8_cells_raw,
  obj_post = cd8_cells,
  prefix = "CD8_pre_post_filter_QC"
)

# 额外输出过滤前后数量表（更适合论文方法补充）
qc_count_table <- data.frame(
  object = c("cd8_cells_raw", "cd8_cells_filtered"),
  n_cells = c(ncol(cd8_cells_raw), ncol(cd8_cells)),
  n_genes = c(nrow(cd8_cells_raw), nrow(cd8_cells))
)
write.csv(qc_count_table, "./temp/CD8_QC_object_size_summary.csv", row.names = FALSE)

# 每组过滤前后保留率表
group_qc_before <- as.data.frame(table(group = cd8_cells_raw$group))
colnames(group_qc_before) <- c("group", "n_pre")

group_qc_after <- as.data.frame(table(group = cd8_cells$group))
colnames(group_qc_after) <- c("group", "n_post")

group_qc_keep <- dplyr::left_join(group_qc_before, group_qc_after, by = "group") %>%
  dplyr::mutate(
    n_post = ifelse(is.na(n_post), 0, n_post),
    keep_rate = n_post / n_pre,
    remove_rate = 1 - keep_rate
  )
write.csv(group_qc_keep, "./temp/CD8_QC_group_keep_rate.csv", row.names = FALSE)

# -------------------------
# 8.5 数量与比例分析
# -------------------------
cd8_ids <- colnames(cd8_cells)

t_by_group <- as.data.frame(table(group = t_cells$group))
colnames(t_by_group) <- c("group", "t_n")

cd8_by_group <- as.data.frame(table(group = cd8_cells$group))
colnames(cd8_by_group) <- c("group", "cd8_n")

qty_t <- dplyr::left_join(cd8_by_group, t_by_group, by = "group") %>%
  dplyr::mutate(cd8_prop_in_T = cd8_n / t_n)

imm_by_group <- as.data.frame(table(group = immune_cells$group))
colnames(imm_by_group) <- c("group", "imm_n")

imm_cells_by_group <- split(colnames(immune_cells), immune_cells$group)

cd8_in_imm_n <- sapply(names(imm_cells_by_group), function(g) {
  length(intersect2(imm_cells_by_group[[g]], cd8_ids))
})

qty_imm <- data.frame(
  group = names(cd8_in_imm_n),
  cd8_in_imm_n = as.integer(cd8_in_imm_n)
) %>%
  dplyr::left_join(imm_by_group, by = "group") %>%
  dplyr::mutate(cd8_prop_in_immune = cd8_in_imm_n / imm_n)

qty_all <- dplyr::full_join(qty_t, qty_imm, by = "group")
write.csv(qty_all, "./temp/CD8_quantity_and_props_point_estimate.csv", row.names = FALSE)

boot_imm <- bootstrap_prop_subsample(
  group_cells_list = imm_cells_by_group,
  cd8_cell_set = cd8_ids,
  frac = 0.8,
  B = 500,
  seed = 123
)
write.csv(boot_imm, "./temp/CD8_prop_in_immune_bootstrap_raw.csv", row.names = FALSE)

boot_sum <- boot_imm %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean_prop = mean(prop),
    sd_prop = sd(prop),
    sample_size = unique(sample_size)[1],
    n_draw = dplyr::n(),
    .groups = "drop"
  )
write.csv(boot_sum, "./temp/CD8_prop_in_immune_bootstrap_summary.csv", row.names = FALSE)

# 8.5b CD8/T bootstrap
t_cells_by_group <- split(colnames(t_cells), t_cells$group)

boot_t <- bootstrap_prop_subsample(
  group_cells_list = t_cells_by_group,
  cd8_cell_set = cd8_ids,
  frac = 0.8,
  B = 500,
  seed = 123
)
write.csv(boot_t, "./temp/CD8_prop_in_T_bootstrap_raw.csv", row.names = FALSE)

boot_t_sum <- boot_t %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    mean_prop = mean(prop),
    sd_prop = sd(prop),
    sample_size = unique(sample_size)[1],
    n_draw = dplyr::n(),
    .groups = "drop"
  )
write.csv(boot_t_sum, "./temp/CD8_prop_in_T_bootstrap_summary.csv", row.names = FALSE)

# -------------------------
# 8.6 first-pass clustering
# -------------------------
cd8_sub_firstpass <- recluster_and_export(
  obj = cd8_cells,
  prefix = "cd8_sub_firstpass",
  group_var = "group",
  dims_use = 1:20,
  resolution_use = 0.6,
  nfeatures_use = 2000
)

# -------------------------
# 8.7 删除 first-pass confirmed contaminants: 4 / 10 / 11
# -------------------------
clusters_exclude_firstpass <- c("4","10","11")
cells_keep_main <- colnames(cd8_sub_firstpass)[!(as.character(cd8_sub_firstpass$seurat_clusters) %in% clusters_exclude_firstpass)]

cd8_sub_firstpass@graphs <- list()
cd8_sub_firstpass@neighbors <- list()

cd8_tex_clean <- subset(cd8_sub_firstpass, cells = cells_keep_main)
cd8_tex_clean$group <- droplevels(cd8_tex_clean$group)

cd8_tex_clean <- recluster_and_export(
  obj = cd8_tex_clean,
  prefix = "cd8_tex_clean_main_revised",
  group_var = "group",
  dims_use = 1:20,
  resolution_use = 0.6,
  nfeatures_use = 2000
)

# -------------------------
# 8.8 删除 revised main 中的 cluster 8（Tumor/Melanocyte-like contamination）
# -------------------------
clusters_exclude_final <- c("8")
cells_keep_final <- colnames(cd8_tex_clean)[!(as.character(cd8_tex_clean$seurat_clusters) %in% clusters_exclude_final)]

cd8_tex_clean@graphs <- list()
cd8_tex_clean@neighbors <- list()

cd8_tex_final <- subset(cd8_tex_clean, cells = cells_keep_final)
cd8_tex_final$group <- droplevels(cd8_tex_final$group)

cd8_tex_final <- recluster_and_export(
  obj = cd8_tex_final,
  prefix = "cd8_tex_final_clean",
  group_var = "group",
  dims_use = 1:20,
  resolution_use = 0.6,
  nfeatures_use = 2000
)

# -------------------------
# 8.9 final clean 注释
# -------------------------
cluster_map_final_clean <- data.frame(
  cluster = c("0","1","2","3","4","5","6","7","8"),
  final_label = c(
    "CD8 activated / memory-like",
    "CD8 Pdcd1hi chronic-stimulated / Tex-like",
    "CD8 cycling-G2/M",
    "CD8 cycling-S phase",
    "CD8 naive/TCM-like",
    "CD8 cytotoxic effector-1",
    "CD8 atypical / stress-like",
    "CD8 activated effector / Xcl1+",
    "CD8 intermediate / transitional"
  ),
  short_label = c(
    "Act/Memory",
    "PD1hi Tex-like",
    "Cycling-G2M",
    "Cycling-S",
    "Naive/TCM",
    "CTL Eff-1",
    "Atypical",
    "Xcl1+ Eff",
    "Transitional"
  ),
  stringsAsFactors = FALSE
)
write.csv(cluster_map_final_clean, "./temp/CD8_cluster_final_annotation_table_final_clean.csv", row.names = FALSE)

export_final_cd8_bundle <- function(obj, cluster_map, prefix = "CD8_final_clean") {
  obj <- ensure_group_and_sample(obj, "cd8_final_obj")
  DefaultAssay(obj) <- "RNA"

  cluster_map$cluster <- as.character(cluster_map$cluster)
  cluster_ids <- as.character(obj$seurat_clusters)

  map_label <- setNames(cluster_map$final_label, cluster_map$cluster)
  map_short <- setNames(cluster_map$short_label, cluster_map$cluster)

  assert_mapping_complete(cluster_ids, map_label, "cluster_map_final_clean")
  assert_mapping_complete(cluster_ids, map_short, "cluster_map_final_clean")

  obj$cluster_label <- base::unname(map_label[cluster_ids])
  obj$cluster_short <- factor(base::unname(map_short[cluster_ids]), levels = cluster_map$short_label)

  saveRDS(obj, "./temp/cd8_final_obj.rds")
  saveRDS(obj, "./temp/cd8_tex_final_clean_annotated.rds")

  comp_labeled <- as.data.frame(table(group = obj$group, cluster = obj$cluster_short)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = Freq / sum(Freq)) %>%
    dplyr::ungroup()
  write.csv(comp_labeled, paste0("./temp/", prefix, "_composition_labeled.csv"), row.names = FALSE)

  cluster_to_function <- c(
    "Act/Memory" = "Activated/Memory-like",
    "PD1hi Tex-like" = "Tex-like",
    "Cycling-G2M" = "Cycling",
    "Cycling-S" = "Cycling",
    "Naive/TCM" = "Naive/TCM",
    "CTL Eff-1" = "Effector",
    "Atypical" = "Other/Transitional",
    "Xcl1+ Eff" = "Effector",
    "Transitional" = "Other/Transitional"
  )
  assert_mapping_complete(as.character(comp_labeled$cluster), cluster_to_function, "cluster_to_function")

  comp_function <- comp_labeled %>%
    dplyr::mutate(function_class = base::unname(cluster_to_function[as.character(cluster)])) %>%
    dplyr::group_by(group, function_class) %>%
    dplyr::summarise(Freq = sum(Freq), .groups = "drop") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(prop = Freq / sum(Freq)) %>%
    dplyr::ungroup()
  write.csv(comp_function, paste0("./temp/", prefix, "_composition_function_merged.csv"), row.names = FALSE)

  Idents(obj) <- "cluster_short"
  p_umap_labeled <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "cluster_short",
    label = TRUE,
    repel = TRUE
  ) + ggtitle("CD8 final clean object (annotated)") + theme_bw()
  ggsave(paste0("./figs/", prefix, "_UMAP_annotated.pdf"), p_umap_labeled, width = 8, height = 6)

  p_comp_labeled <- ggplot(comp_labeled, aes(x = group, y = prop, fill = cluster)) +
    geom_col() + theme_bw() +
    labs(title = "CD8 final clean composition", x = "Group", y = "Proportion", fill = "Cluster")
  ggsave(paste0("./figs/", prefix, "_composition_labeled.pdf"), p_comp_labeled, width = 8, height = 5)

  p_comp_function <- ggplot(comp_function, aes(x = group, y = prop, fill = function_class)) +
    geom_col() + theme_bw() +
    labs(title = "CD8 final clean composition (merged functional classes)", x = "Group", y = "Proportion", fill = "Functional class")
  ggsave(paste0("./figs/", prefix, "_composition_function_merged.pdf"), p_comp_function, width = 8, height = 5)

  marker_panel_final <- present_only(obj, c(
    "Pdcd1","Tox","Nr4a2","Nr4a3","Lag3","Havcr2","Tigit","Ctla4",
    "Entpd1","Layn","Cd101","Cxcl13",
    "Tcf7","Slamf6","Il7r","Lef1","Ccr7","Sell",
    "Cd69","Cxcr6","Itgb7","S1pr1",
    "Gzma","Gzmb","Prf1","Nkg7","Xcl1","Ifng","Ccl5","Gzmk",
    "Trgc2"
  ))

  avg_final_clean <- AverageExpression(obj, group.by = "seurat_clusters", assays = "RNA", slot = "data")$RNA
  avg_final_clean_sub <- avg_final_clean[rownames(avg_final_clean) %in% marker_panel_final, , drop = FALSE]
  write.csv(avg_final_clean_sub, paste0("./temp/", prefix, "_avgexpr_marker_panel.csv"))

  p_dot_final <- dotplot3c(obj, features = marker_panel_final, group.by = "seurat_clusters") +
    RotatedAxis() + theme_bw() +
    labs(title = "CD8 final clean object: marker evidence")
  ggsave(paste0("./figs/", prefix, "_marker_dotplot.pdf"), p_dot_final, width = 16, height = 6)

  # -------------------------
  # 8.11b cluster identity validation panel
  # -------------------------
  identity_panel <- present_only(obj, c(
    "Cd3d","Cd3e","Trac","Cd247","Lck",
    "Cd8a","Cd8b1",
    "Tcf7","Lef1","Il7r","Ccr7","Sell",
    "Trdc","Trgc1","Trgc2",
    "Ncr1","Klrd1","Klrk1","Tyrobp","Fcgr3",
    "Lyz2","Itgam","Lst1","S100a8","S100a9",
    "Clec9a","Flt3","Xcr1","Itgax",
    "Foxp3","Il2ra","Ikzf2",
    "Pmel","Mlana","Tyr","Dct","Slc45a2"
  ))

  avg_identity <- AverageExpression(obj, group.by = "cluster_short", assays = "RNA", slot = "data")$RNA
  avg_identity_sub <- avg_identity[rownames(avg_identity) %in% identity_panel, , drop = FALSE]
  write.csv(avg_identity_sub, "./temp/CD8_identity_validation_avgexpr.csv")

  p_identity_dot <- dotplot3c(obj, features = identity_panel, group.by = "cluster_short") +
    RotatedAxis() + theme_bw() +
    labs(title = "CD8 final object: cluster identity validation panel",
         x = "Cluster", y = "Identity / contamination markers")
  ggsave("./figs/CD8_identity_validation_dotplot.pdf", p_identity_dot, width = 18, height = 7)

  identity_panel_core <- present_only(obj, c(
    "Cd3d","Cd3e","Trac","Cd8a","Cd8b1",
    "Trdc","Trgc1","Trgc2",
    "Ncr1","Klrd1","Tyrobp","Lyz2"
  ))

  p_identity_dot_core <- dotplot3c(obj, features = identity_panel_core, group.by = "cluster_short") +
    RotatedAxis() + theme_bw() +
    labs(title = "CD8 final object: core identity validation panel",
         x = "Cluster", y = "Core lineage markers")
  ggsave("./figs/CD8_identity_validation_core_dotplot.pdf", p_identity_dot_core, width = 14, height = 6)

  identity_group_panel <- present_only(obj, c(
    "Cd3d","Cd3e","Trac","Cd8a","Cd8b1",
    "Trdc","Trgc2",
    "Ncr1","Klrd1","Tyrobp",
    "Lyz2","Clec9a",
    "Pmel","Mlana"
  ))

  p_identity_group <- dotplot3c(obj, features = identity_group_panel, group.by = "group") +
    RotatedAxis() + theme_bw() +
    labs(title = "Group-wise lineage contamination check",
         x = "Group", y = "Lineage / contamination markers")
  ggsave("./figs/CD8_identity_groupwise_dotplot.pdf", p_identity_group, width = 14, height = 6)

  score_signature <- function(obj2, genes, score_name) {
    genes_use <- base::intersect(genes, rownames(obj2))
    if (length(genes_use) == 0) {
      obj2[[score_name]] <- 0
      return(obj2)
    }
    obj2 <- AddModuleScore(obj2, features = list(genes_use), name = score_name)
    return(obj2)
  }

  obj <- score_signature(obj, c("Cd3d","Cd3e","Trac","Lck"), "PanTScore")
  obj <- score_signature(obj, c("Cd8a","Cd8b1"), "CD8Score")
  obj <- score_signature(obj, c("Trdc","Trgc1","Trgc2"), "GDTScore")
  obj <- score_signature(obj, c("Ncr1","Klrd1","Klrk1","Tyrobp"), "NKScore")
  obj <- score_signature(obj, c("Lyz2","Itgam","Lst1","S100a8","S100a9"), "MyeloidScore")
  obj <- score_signature(obj, c("Clec9a","Flt3","Xcr1","Itgax"), "DCScore")
  obj <- score_signature(obj, c("Foxp3","Il2ra","Ikzf2"), "TregScore")
  obj <- score_signature(obj, c("Pmel","Mlana","Tyr","Dct","Slc45a2"), "TumorScore")

  contam_cols <- c(
    "PanTScore1","CD8Score1","GDTScore1","NKScore1",
    "MyeloidScore1","DCScore1","TregScore1","TumorScore1"
  )
  contam_cols <- contam_cols[contam_cols %in% colnames(obj@meta.data)]

  cluster_contam_summary <- obj@meta.data %>%
    dplyr::group_by(cluster_short) %>%
    dplyr::summarise(
      dplyr::across(all_of(contam_cols), ~ mean(.x, na.rm = TRUE)),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
  write.csv(cluster_contam_summary, "./temp/CD8_cluster_identity_score_summary.csv", row.names = FALSE)
  saveRDS(obj, "./temp/cd8_final_obj_identity_scored.rds")

  return(obj)
}

cd8_final_obj <- export_final_cd8_bundle(
  obj = cd8_tex_final,
  cluster_map = cluster_map_final_clean,
  prefix = "CD8_final_clean"
)

# 为了兼容后面旧代码，这里把 cd8_tex_final 直接替换成唯一的 final object
cd8_tex_final <- cd8_final_obj

# =========================================================
# 8.11f DoRothEA TF activity analysis
# 放置位置：
# 在 cd8_final_obj <- export_final_cd8_bundle(...) 之后
# 在 cd8_tex_final <- cd8_final_obj 之后
# 在 8.12 cluster1 之前
# =========================================================

# ---------- 0) 基础检查 ----------
if (!requireNamespace("dorothea", quietly = TRUE)) {
  stop("缺少 dorothea 包。请先安装：BiocManager::install('dorothea')")
}
if (!requireNamespace("viper", quietly = TRUE)) {
  stop("缺少 viper 包。请先安装：BiocManager::install('viper')")
}

# 使用 mouse regulons；先用高置信度 A/B/C，避免太噪
data(dorothea_mm, package = "dorothea")

regulon_mm_abc <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

write.csv(regulon_mm_abc, "./temp/dorothea_mm_ABC_regulons.csv", row.names = FALSE)

# ---------- 1) 在 final CD8 object 上运行 DoRothEA ----------
# run_viper 官方支持 Seurat 对象输入；会在对象中添加 dorothea assay
cd8_final_obj <- dorothea::run_viper(
  input = cd8_final_obj,
  regulons = regulon_mm_abc,
  options = list(
    method = "scale",
    minsize = 4,
    eset.filter = FALSE,
    cores = 1,
    verbose = FALSE
  ),
  assay_key = "RNA"
)

saveRDS(cd8_final_obj, "./temp/cd8_final_obj_dorothea.rds")

# ---------- 2) 提取 TF activity matrix ----------
if (!("dorothea" %in% names(cd8_final_obj@assays))) {
  stop("DoRothEA 运行后没有生成 dorothea assay，请检查 run_viper 是否成功。")
}

tf_mat <- GetAssayData(cd8_final_obj, assay = "dorothea", slot = "data")
write.csv(as.matrix(tf_mat), "./temp/CD8_final_dorothea_activity_matrix.csv")

# ---------- 3) 导出：按 cluster_short 的平均 TF 活性 ----------
avg_tf_cluster <- AverageExpression(
  cd8_final_obj,
  assays = "dorothea",
  group.by = "cluster_short",
  slot = "data"
)$dorothea
write.csv(avg_tf_cluster, "./temp/CD8_dorothea_avgTF_by_cluster.csv")

# 选最值得看的 TF
tf_focus <- c(
  "Tcf7","Lef1","Bach2","Stat5a","Stat5b",
  "Tox","Nr4a1","Nr4a2","Nr4a3",
  "Tbx21","Eomes","Jun","Fos","Hif1a"
)
tf_focus <- tf_focus[tf_focus %in% rownames(tf_mat)]

if (length(tf_focus) > 0) {
  avg_tf_cluster_sub <- avg_tf_cluster[rownames(avg_tf_cluster) %in% tf_focus, , drop = FALSE]
  write.csv(avg_tf_cluster_sub, "./temp/CD8_dorothea_avgTF_by_cluster_focus.csv")

  p_tf_cluster <- dotplot3c(
    cd8_final_obj,
    assay = "dorothea",
    features = tf_focus,
    group.by = "cluster_short"
  ) +
    RotatedAxis() +
    theme_bw() +
    labs(
      title = "DoRothEA TF activity by cluster",
      x = "Cluster",
      y = "TF activity"
    )

  ggsave("./figs/CD8_dorothea_TF_dotplot_by_cluster.pdf", p_tf_cluster, width = 14, height = 6)
}

# ---------- 4) 导出：按 group + cluster_short 的 TF activity 长表 ----------
tf_meta <- t(as.matrix(tf_mat))
tf_meta <- as.data.frame(tf_meta)
tf_meta$cell <- rownames(tf_meta)

meta_use <- cd8_final_obj@meta.data[, c("group", "cluster_short"), drop = FALSE]
meta_use$cell <- rownames(meta_use)

tf_long <- dplyr::left_join(meta_use, tf_meta, by = "cell") %>%
  tidyr::pivot_longer(
    cols = -c(cell, group, cluster_short),
    names_to = "TF",
    values_to = "activity"
  )

write.csv(tf_long, "./temp/CD8_dorothea_TF_activity_long.csv", row.names = FALSE)

# ---------- 5) cluster0 和 cluster1 专门导出 ----------
for (cl in c("0", "1")) {
  obj_sub <- subset(cd8_final_obj, subset = seurat_clusters == cl)
  prefix_sub <- ifelse(cl == "0", "cluster0", "cluster1")

  tf_sub <- GetAssayData(obj_sub, assay = "dorothea", slot = "data")
  tf_sub_df <- t(as.matrix(tf_sub))
  tf_sub_df <- as.data.frame(tf_sub_df)
  tf_sub_df$cell <- rownames(tf_sub_df)

  meta_sub <- obj_sub@meta.data[, c("group"), drop = FALSE]
  meta_sub$cell <- rownames(meta_sub)

  tf_sub_long <- dplyr::left_join(meta_sub, tf_sub_df, by = "cell") %>%
    tidyr::pivot_longer(
      cols = -c(cell, group),
      names_to = "TF",
      values_to = "activity"
    )

  write.csv(tf_sub_long, paste0("./temp/", prefix_sub, "_dorothea_TF_activity_long.csv"), row.names = FALSE)

  tf_focus_sub <- tf_focus[tf_focus %in% colnames(tf_sub_df)]
  if (length(tf_focus_sub) > 0) {
    tf_sub_long_focus <- tf_sub_long %>% dplyr::filter(TF %in% tf_focus_sub)

    p_tf_sub <- ggplot(tf_sub_long_focus, aes(x = group, y = activity)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      facet_wrap(~ TF, scales = "free_y", ncol = 3) +
      theme_bw() +
      labs(
        title = paste0(prefix_sub, ": DoRothEA TF activity by group"),
        x = "Group",
        y = "TF activity"
      )

    ggsave(paste0("./figs/", prefix_sub, "_dorothea_TF_activity_by_group.pdf"), p_tf_sub, width = 12, height = 8)
  }
}

# ---------- 6) cluster1 subcluster / pseudotime 预留接口 ----------
# 注意：这一步不强行运行，只是预留。等 cluster1_sub 生成后可以追加：
# - 对 cluster1_sub 重新 run_viper，或直接 subset final object 的 dorothea assay
# - 将 TF activity 与 pseudotime 关联分析



# -------------------------
# 8.11e pathway-like program scores for pooled-sample analysis
# -------------------------
score_programs_for_object <- function(obj, program_list, prefix = "Prog_") {
  for (nm in names(program_list)) {
    genes_use <- base::intersect(program_list[[nm]], rownames(obj))
    if (length(genes_use) == 0) {
      obj[[paste0(prefix, nm)]] <- 0
    } else {
      obj <- AddModuleScore(obj, features = list(genes_use), name = paste0(prefix, nm))
    }
  }
  obj
}

tcell_programs <- list(
  IL2_STAT5 = c("Il2ra","Il2rb","Stat5a","Stat5b","Bcl2","Cish","Pim1","Socs2"),
  TCF7_LEF1 = c("Tcf7","Lef1","Il7r","Ccr7","Sell","Bach2","S1pr1"),
  TOX_NR4A = c("Tox","Nr4a1","Nr4a2","Nr4a3","Pdcd1","Ctla4","Havcr2","Tigit"),
  CYTOTOXICITY = c("Prf1","Gzmb","Gzma","Ifng","Nkg7","Xcl1","Ccl5"),
  HYPOXIA_GLYCOLYSIS = c("Hif1a","Slc2a1","Ldha","Pgk1","Tpi1","Aldoa","Gapdh","Eno1")
)

cd8_final_obj <- score_programs_for_object(cd8_final_obj, tcell_programs, prefix = "Prog_")
saveRDS(cd8_final_obj, "./temp/cd8_final_obj_program_scored.rds")
cd8_tex_final <- cd8_final_obj

prog_cols <- c("Prog_IL2_STAT51","Prog_TCF7_LEF11","Prog_TOX_NR4A1","Prog_CYTOTOXICITY1","Prog_HYPOXIA_GLYCOLYSIS1")
prog_cols <- prog_cols[prog_cols %in% colnames(cd8_final_obj@meta.data)]

if (length(prog_cols) > 0) {
  df_prog <- cd8_final_obj@meta.data %>%
    dplyr::select(group, cluster_short, dplyr::all_of(prog_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(prog_cols), names_to = "program", values_to = "score")
  write.csv(df_prog, "./temp/CD8_final_program_scores_long.csv", row.names = FALSE)

  p_prog <- ggplot(df_prog, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_grid(program ~ cluster_short, scales = "free_y") +
    theme_bw() +
    labs(title = "CD8 final object: pathway-like program scores by group and cluster",
         x = "Group", y = "Module score")
  ggsave("./figs/CD8_final_program_scores_by_group_cluster.pdf", p_prog, width = 18, height = 10)
}

cluster0_obj <- subset(cd8_final_obj, subset = seurat_clusters == "0")
cluster0_obj <- score_programs_for_object(cluster0_obj, tcell_programs, prefix = "Prog_")
saveRDS(cluster0_obj, "./temp/cluster0_obj_program_scored.rds")

cluster1_obj_for_prog <- subset(cd8_final_obj, subset = seurat_clusters == "1")
cluster1_obj_for_prog <- score_programs_for_object(cluster1_obj_for_prog, tcell_programs, prefix = "Prog_")
saveRDS(cluster1_obj_for_prog, "./temp/cluster1_obj_program_scored.rds")

for (nm in c("cluster0_obj", "cluster1_obj_for_prog")) {
  objx <- get(nm)
  prog_cols2 <- c("Prog_IL2_STAT51","Prog_TCF7_LEF11","Prog_TOX_NR4A1","Prog_CYTOTOXICITY1","Prog_HYPOXIA_GLYCOLYSIS1")
  prog_cols2 <- prog_cols2[prog_cols2 %in% colnames(objx@meta.data)]

  if (length(prog_cols2) > 0) {
    dfx <- objx@meta.data %>%
      dplyr::select(group, dplyr::all_of(prog_cols2)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(prog_cols2), names_to = "program", values_to = "score")
    write.csv(dfx, paste0("./temp/", nm, "_program_scores_long.csv"), row.names = FALSE)

    px <- ggplot(dfx, aes(x = group, y = score)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      facet_wrap(~ program, scales = "free_y", ncol = 3) +
      theme_bw() +
      labs(title = paste0(nm, ": pathway-like program scores by group"),
           x = "Group", y = "Module score")
    ggsave(paste0("./figs/", nm, "_program_scores_by_group.pdf"), px, width = 12, height = 7)
  }
}
# =========================================================
# 8.11e pathway-like / regulon-like program scoring for cluster0 & cluster1
# 放置位置：
# 在 cd8_final_obj <- export_final_cd8_bundle(...) 之后
# 在 cd8_tex_final <- cd8_final_obj 之后
# 在 8.12 cluster1 之前
# =========================================================

score_programs_for_object <- function(obj, program_list, prefix = "Prog_") {
  for (nm in names(program_list)) {
    genes_use <- base::intersect(program_list[[nm]], rownames(obj))
    if (length(genes_use) == 0) {
      obj[[paste0(prefix, nm)]] <- 0
    } else {
      obj <- AddModuleScore(obj, features = list(genes_use), name = paste0(prefix, nm))
    }
  }
  return(obj)
}

# -------- 1) 五条主线程序 --------
tcell_programs <- list(
  IL2_STAT5 = c("Il2ra","Il2rb","Stat5a","Stat5b","Bcl2","Cish","Pim1","Socs2"),
  TCF7_LEF1 = c("Tcf7","Lef1","Il7r","Ccr7","Sell","Bach2","S1pr1"),
  TOX_NR4A = c("Tox","Nr4a1","Nr4a2","Nr4a3","Pdcd1","Ctla4","Havcr2","Tigit"),
  CYTOTOXICITY = c("Prf1","Gzmb","Gzma","Ifng","Nkg7","Xcl1","Ccl5"),
  HYPOXIA_GLYCOLYSIS = c("Hif1a","Slc2a1","Ldha","Pgk1","Tpi1","Aldoa","Gapdh","Eno1")
)

# -------- 2) 在 final CD8 object 上打分 --------
cd8_final_obj <- score_programs_for_object(cd8_final_obj, tcell_programs, prefix = "Prog_")
saveRDS(cd8_final_obj, "./temp/cd8_final_obj_program_scored.rds")

prog_cols <- c(
  "Prog_IL2_STAT51",
  "Prog_TCF7_LEF11",
  "Prog_TOX_NR4A1",
  "Prog_CYTOTOXICITY1",
  "Prog_HYPOXIA_GLYCOLYSIS1"
)
prog_cols <- prog_cols[prog_cols %in% colnames(cd8_final_obj@meta.data)]

# -------- 3) 全局导出：按 group + cluster_short 比较 --------
if (length(prog_cols) > 0) {
  df_prog <- cd8_final_obj@meta.data %>%
    dplyr::select(group, cluster_short, dplyr::all_of(prog_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(prog_cols),
      names_to = "program",
      values_to = "score"
    )

  write.csv(df_prog, "./temp/CD8_final_program_scores_long.csv", row.names = FALSE)

  p_prog <- ggplot(df_prog, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_grid(program ~ cluster_short, scales = "free_y") +
    theme_bw() +
    labs(
      title = "CD8 final object: pathway-like program scores by group and cluster",
      x = "Group",
      y = "Module score"
    )

  ggsave("./figs/CD8_final_program_scores_by_group_cluster.pdf", p_prog, width = 18, height = 10)
}

# -------- 4) cluster0 专门分析 --------
cluster0_obj_prog <- subset(cd8_final_obj, subset = seurat_clusters == "0")
cluster0_obj_prog$group <- droplevels(cluster0_obj_prog$group)
cluster0_obj_prog <- score_programs_for_object(cluster0_obj_prog, tcell_programs, prefix = "Prog_")
saveRDS(cluster0_obj_prog, "./temp/cluster0_obj_program_scored.rds")

prog_cols0 <- c(
  "Prog_IL2_STAT51",
  "Prog_TCF7_LEF11",
  "Prog_TOX_NR4A1",
  "Prog_CYTOTOXICITY1",
  "Prog_HYPOXIA_GLYCOLYSIS1"
)
prog_cols0 <- prog_cols0[prog_cols0 %in% colnames(cluster0_obj_prog@meta.data)]

if (length(prog_cols0) > 0) {
  df0 <- cluster0_obj_prog@meta.data %>%
    dplyr::select(group, dplyr::all_of(prog_cols0)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(prog_cols0),
      names_to = "program",
      values_to = "score"
    )

  write.csv(df0, "./temp/cluster0_program_scores_long.csv", row.names = FALSE)

  p0 <- ggplot(df0, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ program, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster 0 (Act/Memory): pathway-like program scores by group",
      x = "Group",
      y = "Module score"
    )

  ggsave("./figs/cluster0_program_scores_by_group.pdf", p0, width = 12, height = 7)
}

# -------- 5) cluster1 专门分析 --------
cluster1_obj_prog <- subset(cd8_final_obj, subset = seurat_clusters == "1")
cluster1_obj_prog$group <- droplevels(cluster1_obj_prog$group)
cluster1_obj_prog <- score_programs_for_object(cluster1_obj_prog, tcell_programs, prefix = "Prog_")
saveRDS(cluster1_obj_prog, "./temp/cluster1_obj_program_scored.rds")

prog_cols1 <- c(
  "Prog_IL2_STAT51",
  "Prog_TCF7_LEF11",
  "Prog_TOX_NR4A1",
  "Prog_CYTOTOXICITY1",
  "Prog_HYPOXIA_GLYCOLYSIS1"
)
prog_cols1 <- prog_cols1[prog_cols1 %in% colnames(cluster1_obj_prog@meta.data)]

if (length(prog_cols1) > 0) {
  df1 <- cluster1_obj_prog@meta.data %>%
    dplyr::select(group, dplyr::all_of(prog_cols1)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(prog_cols1),
      names_to = "program",
      values_to = "score"
    )

  write.csv(df1, "./temp/cluster1_program_scores_long.csv", row.names = FALSE)

  p1 <- ggplot(df1, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ program, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster 1 (PD1hi Tex-like): pathway-like program scores by group",
      x = "Group",
      y = "Module score"
    )

  ggsave("./figs/cluster1_program_scores_by_group.pdf", p1, width = 12, height = 7)
}

# =========================================================
# 8.11f simplified TF-axis scoring (regulon-like, lightweight version)
# 放在 pathway-like module score 后面
# =========================================================

tf_axis_list <- list(
  TCF7_LEF1_TF = c("Tcf7","Lef1","Bach2","Klf2","S1pr1","Il7r"),
  TOX_NR4A_TF = c("Tox","Nr4a1","Nr4a2","Nr4a3","Pdcd1","Ctla4","Tigit"),
  AP1_ACTIVATION = c("Jun","Fos","Fosb","Jund","Junb"),
  TBX21_EOMES = c("Tbx21","Eomes","Ifng","Cxcr3","Prf1"),
  HIF1A_GLYCOLYSIS = c("Hif1a","Slc2a1","Ldha","Pgk1","Aldoa")
)

score_programs_for_object <- function(obj, program_list, prefix = "TF_") {
  for (nm in names(program_list)) {
    genes_use <- base::intersect(program_list[[nm]], rownames(obj))
    if (length(genes_use) == 0) {
      obj[[paste0(prefix, nm)]] <- 0
    } else {
      obj <- AddModuleScore(obj, features = list(genes_use), name = paste0(prefix, nm))
    }
  }
  return(obj)
}

cluster0_obj_prog <- score_programs_for_object(cluster0_obj_prog, tf_axis_list, prefix = "TF_")
cluster1_obj_prog <- score_programs_for_object(cluster1_obj_prog, tf_axis_list, prefix = "TF_")

tf_cols <- c(
  "TF_TCF7_LEF1_TF1",
  "TF_TOX_NR4A_TF1",
  "TF_AP1_ACTIVATION1",
  "TF_TBX21_EOMES1",
  "TF_HIF1A_GLYCOLYSIS1"
)

tf_cols0 <- tf_cols[tf_cols %in% colnames(cluster0_obj_prog@meta.data)]
if (length(tf_cols0) > 0) {
  df_tf0 <- cluster0_obj_prog@meta.data %>%
    dplyr::select(group, dplyr::all_of(tf_cols0)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(tf_cols0),
      names_to = "tf_axis",
      values_to = "score"
    )

  write.csv(df_tf0, "./temp/cluster0_tf_axis_scores_long.csv", row.names = FALSE)

  p_tf0 <- ggplot(df_tf0, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ tf_axis, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster 0: simplified TF-axis scores by group",
      x = "Group",
      y = "Module score"
    )

  ggsave("./figs/cluster0_tf_axis_scores_by_group.pdf", p_tf0, width = 12, height = 7)
}

tf_cols1 <- tf_cols[tf_cols %in% colnames(cluster1_obj_prog@meta.data)]
if (length(tf_cols1) > 0) {
  df_tf1 <- cluster1_obj_prog@meta.data %>%
    dplyr::select(group, dplyr::all_of(tf_cols1)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(tf_cols1),
      names_to = "tf_axis",
      values_to = "score"
    )

  write.csv(df_tf1, "./temp/cluster1_tf_axis_scores_long.csv", row.names = FALSE)

  p_tf1 <- ggplot(df_tf1, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ tf_axis, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster 1: simplified TF-axis scores by group",
      x = "Group",
      y = "Module score"
    )

  ggsave("./figs/cluster1_tf_axis_scores_by_group.pdf", p_tf1, width = 12, height = 7)
}
# -------------------------
# 8.12 聚焦 cluster 1：Tex/Tpex/terminal Tex 打磨
# -------------------------
genes_tex_core <- c("Pdcd1","Tox","Nr4a2","Nr4a3","Lag3","Havcr2","Tigit","Ctla4")
genes_terminal <- c("Entpd1","Layn","Cd101","Cxcl13","Tnfrsf9")
genes_tpex <- c("Tcf7","Slamf6","Cxcr5","Il7r","Lef1","Bach2")
genes_memory_act <- c("Cd69","Cxcr6","S1pr1","Itgb7","Gzma")
genes_effector <- c("Gzmb","Prf1","Nkg7","Ifng","Xcl1")
genes_misc <- c("Trgc2")

marker_panel_tex_refine <- unique(c(
  genes_tex_core, genes_terminal, genes_tpex, genes_memory_act, genes_effector, genes_misc
))
marker_panel_tex_refine <- marker_panel_tex_refine[marker_panel_tex_refine %in% rownames(cd8_tex_final)]

avg_tex_refine <- AverageExpression(
  cd8_tex_final,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"
)$RNA

avg_tex_refine_sub <- avg_tex_refine[
  rownames(avg_tex_refine) %in% marker_panel_tex_refine,
  ,
  drop = FALSE
]
write.csv(avg_tex_refine_sub, "./temp/CD8_cluster1_tex_refine_avgexpr.csv")

p_dot_tex_refine <- dotplot3c(
  cd8_tex_final,
  features = marker_panel_tex_refine,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme_bw() +
  labs(title = "CD8 final clean: Tex/Tpex refinement markers")

ggsave("./figs/CD8_final_clean_Tex_Tpex_dotplot.pdf", p_dot_tex_refine, width = 16, height = 6)

cluster1_cells <- colnames(cd8_tex_final)[as.character(cd8_tex_final$seurat_clusters) == "1"]

if (length(cluster1_cells) > 0) {
  cd8_tex_final$cluster1_flag <- ifelse(colnames(cd8_tex_final) %in% cluster1_cells, "cluster1", "others")

  avg_cluster1_vs_others <- AverageExpression(
    cd8_tex_final,
    group.by = "cluster1_flag",
    assays = "RNA",
    slot = "data"
  )$RNA

  avg_cluster1_vs_others_sub <- avg_cluster1_vs_others[
    rownames(avg_cluster1_vs_others) %in% marker_panel_tex_refine,
    ,
    drop = FALSE
  ]
  write.csv(avg_cluster1_vs_others_sub, "./temp/CD8_cluster1_vs_others_tex_refine_avgexpr.csv")

  Idents(cd8_tex_final) <- "cluster1_flag"
  markers_cluster1_vs_others <- FindMarkers(
    cd8_tex_final,
    ident.1 = "cluster1",
    ident.2 = "others",
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.10,
    logfc.threshold = 0.10
  )
  write.csv(markers_cluster1_vs_others, "./temp/CD8_cluster1_vs_others_markers.csv")
}

cluster1_obj <- subset(cd8_tex_final, subset = seurat_clusters == "1")

if (ncol(cluster1_obj) > 0) {
  cluster1_obj$group <- factor(cluster1_obj$group)

  avg_cluster1_by_group <- AverageExpression(
    cluster1_obj,
    group.by = "group",
    assays = "RNA",
    slot = "data"
  )$RNA

  avg_cluster1_by_group_sub <- avg_cluster1_by_group[
    rownames(avg_cluster1_by_group) %in% marker_panel_tex_refine,
    ,
    drop = FALSE
  ]
  write.csv(avg_cluster1_by_group_sub, "./temp/CD8_cluster1_by_group_tex_refine_avgexpr.csv")

  p_dot_cluster1_group <- dotplot3c(
    cluster1_obj,
    features = marker_panel_tex_refine,
    group.by = "group"
  ) +
    RotatedAxis() +
    theme_bw() +
    labs(title = "Cluster 1 only: Tex/Tpex markers by group")

  ggsave("./figs/CD8_cluster1_by_group_Tex_Tpex_dotplot.pdf", p_dot_cluster1_group, width = 12, height = 6)

  tex_axis <- base::intersect(c("Pdcd1","Tox","Nr4a2","Nr4a3","Ctla4","Havcr2","Tigit"), rownames(cluster1_obj))
  tpex_axis <- base::intersect(c("Tcf7","Slamf6","Il7r","Lef1","Bach2"), rownames(cluster1_obj))
  terminal_axis <- base::intersect(c("Entpd1","Layn","Cd101","Cxcl13","Tnfrsf9"), rownames(cluster1_obj))
  actmem_axis <- base::intersect(c("Gzma","S1pr1","Il7r","Cd69","Cxcr6","Itgb7"), rownames(cluster1_obj))
  effector_axis <- base::intersect(c("Prf1","Gzmb","Nkg7","Ifng","Xcl1","Ccl5"), rownames(cluster1_obj))

  if (length(tex_axis) > 0) cluster1_obj <- AddModuleScore(cluster1_obj, features = list(tex_axis), name = "TexAxis")
  if (length(tpex_axis) > 0) cluster1_obj <- AddModuleScore(cluster1_obj, features = list(tpex_axis), name = "TpexAxis")
  if (length(terminal_axis) > 0) cluster1_obj <- AddModuleScore(cluster1_obj, features = list(terminal_axis), name = "TerminalAxis")
  if (length(actmem_axis) > 0) cluster1_obj <- AddModuleScore(cluster1_obj, features = list(actmem_axis), name = "ActMemAxis")
  if (length(effector_axis) > 0) cluster1_obj <- AddModuleScore(cluster1_obj, features = list(effector_axis), name = "EffectorAxis")

  saveRDS(cluster1_obj, "./temp/cluster1_obj_group_scored.rds")

  score_cols <- base::intersect(c("TexAxis1","TpexAxis1","TerminalAxis1","ActMemAxis1","EffectorAxis1"), colnames(cluster1_obj@meta.data))

  if (length(score_cols) > 0) {
    cluster1_score_long <- cluster1_obj@meta.data %>%
      dplyr::select(group, all_of(score_cols)) %>%
      tidyr::pivot_longer(cols = all_of(score_cols), names_to = "axis", values_to = "score")

    write.csv(cluster1_score_long, "./temp/CD8_cluster1_function_axes_scores_long.csv", row.names = FALSE)

    p_cluster1_scores <- ggplot(cluster1_score_long, aes(x = group, y = score)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA) +
      facet_wrap(~ axis, scales = "free_y", ncol = 3) +
      theme_bw() +
      labs(title = "Cluster 1 functional axes by group", x = "Group", y = "Module score")

      ggsave("./figs/CD8_cluster1_function_axes_by_group.pdf", p_cluster1_scores, width = 12, height = 7)
  }
}

# -------------------------
# 8.13 聚焦 cycling：组间比较
# -------------------------
if (file.exists("./temp/cycling_obj_final_clean_scored.rds")) {
  cycling_obj <- readRDS("./temp/cycling_obj_final_clean_scored.rds")
  message("Loaded ./temp/cycling_obj_final_clean_scored.rds")
} else if (file.exists("./temp/cycling_obj_final_clean.rds")) {
  cycling_obj <- readRDS("./temp/cycling_obj_final_clean.rds")
  message("Loaded ./temp/cycling_obj_final_clean.rds")
} else {
  cycling_obj <- subset(cd8_tex_final, subset = seurat_clusters %in% c("2", "3"))
  cycling_obj$group <- droplevels(cycling_obj$group)
  saveRDS(cycling_obj, "./temp/cycling_obj_final_clean.rds")
}

cycling_comp <- as.data.frame(table(group = cycling_obj$group, cluster = cycling_obj$seurat_clusters)) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(prop = Freq / sum(Freq)) %>%
  dplyr::ungroup()

write.csv(cycling_comp, "./temp/CD8_cycling_composition_by_group.csv", row.names = FALSE)

p_cycling_comp <- ggplot(cycling_comp, aes(x = group, y = prop, fill = cluster)) +
  geom_col() +
  theme_bw() +
  labs(title = "Cycling pool composition by group", x = "Group", y = "Proportion", fill = "Cycling cluster")
ggsave("./figs/CD8_cycling_composition_by_group.pdf", p_cycling_comp, width = 7, height = 5)

cycling_marker_panel <- present_only(cycling_obj, c(
  "Mki67","Top2a","Ccna2","Cdk1","Pbk","Mcm2","Mcm5","Dhfr","Ung",
  "Pdcd1","Tox","Nr4a2","Nr4a3","Ctla4","Havcr2",
  "Tcf7","Slamf6","Il7r","Lef1","Ccr7","Sell",
  "Gzma","S1pr1","Cxcr6",
  "Prf1","Gzmb","Nkg7","Ifng","Xcl1"
))

avg_cycling_by_group <- AverageExpression(
  cycling_obj,
  group.by = "group",
  assays = "RNA",
  slot = "data"
)$RNA

avg_cycling_by_group_sub <- avg_cycling_by_group[
  rownames(avg_cycling_by_group) %in% cycling_marker_panel,
  ,
  drop = FALSE
]
write.csv(avg_cycling_by_group_sub, "./temp/CD8_cycling_by_group_avgexpr.csv")

p_dot_cycling_group <- dotplot3c(
  cycling_obj,
  features = cycling_marker_panel,
  group.by = "group"
) +
  RotatedAxis() +
  theme_bw() +
  labs(title = "Cycling pool: marker comparison by group")

ggsave("./figs/CD8_cycling_by_group_dotplot.pdf", p_dot_cycling_group, width = 13, height = 6)

need_cycling_scores <- !all(c("CyclingTexAxis1","CyclingTpexAxis1","CyclingActMemAxis1","CyclingEffectorAxis1","CyclingCellCycleAxis1") %in%
                             colnames(cycling_obj@meta.data))

if (need_cycling_scores) {
  cycling_tex_axis <- base::intersect(c("Pdcd1","Tox","Nr4a2","Nr4a3","Ctla4","Havcr2"), rownames(cycling_obj))
  cycling_tpex_axis <- base::intersect(c("Tcf7","Slamf6","Il7r","Lef1","Bach2"), rownames(cycling_obj))
  cycling_actmem_axis <- base::intersect(c("Gzma","S1pr1","Cxcr6","Il7r"), rownames(cycling_obj))
  cycling_effector_axis <- base::intersect(c("Prf1","Gzmb","Nkg7","Ifng","Xcl1"), rownames(cycling_obj))
  cycling_cycle_axis <- base::intersect(c("Mki67","Top2a","Ccna2","Cdk1","Pbk","Mcm2","Mcm5","Dhfr"), rownames(cycling_obj))

  if (length(cycling_tex_axis) > 0) cycling_obj <- AddModuleScore(cycling_obj, features = list(cycling_tex_axis), name = "CyclingTexAxis")
  if (length(cycling_tpex_axis) > 0) cycling_obj <- AddModuleScore(cycling_obj, features = list(cycling_tpex_axis), name = "CyclingTpexAxis")
  if (length(cycling_actmem_axis) > 0) cycling_obj <- AddModuleScore(cycling_obj, features = list(cycling_actmem_axis), name = "CyclingActMemAxis")
  if (length(cycling_effector_axis) > 0) cycling_obj <- AddModuleScore(cycling_obj, features = list(cycling_effector_axis), name = "CyclingEffectorAxis")
  if (length(cycling_cycle_axis) > 0) cycling_obj <- AddModuleScore(cycling_obj, features = list(cycling_cycle_axis), name = "CyclingCellCycleAxis")

  saveRDS(cycling_obj, "./temp/cycling_obj_final_clean_scored.rds")
}

cycling_score_cols <- base::intersect(
  c("CyclingTexAxis1","CyclingTpexAxis1","CyclingActMemAxis1","CyclingEffectorAxis1","CyclingCellCycleAxis1"),
  colnames(cycling_obj@meta.data)
)

if (length(cycling_score_cols) > 0) {
  cycling_score_long <- cycling_obj@meta.data %>%
    dplyr::select(group, all_of(cycling_score_cols)) %>%
    tidyr::pivot_longer(cols = all_of(cycling_score_cols), names_to = "axis", values_to = "score")

  write.csv(cycling_score_long, "./temp/CD8_cycling_function_axes_scores_long.csv", row.names = FALSE)

  p_cycling_scores <- ggplot(cycling_score_long, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ axis, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(title = "Cycling pool functional axes by group", x = "Group", y = "Module score")

  ggsave("./figs/CD8_cycling_function_axes_by_group.pdf", p_cycling_scores, width = 12, height = 7)
}

# =========================================================
# 9. cluster 1 伪时间分析（完整版，保留原逻辑）
# =========================================================
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("缺少 SingleCellExperiment 包。请先安装：BiocManager::install('SingleCellExperiment')")
}
if (!requireNamespace("slingshot", quietly = TRUE)) {
  stop("缺少 slingshot 包。请先安装：BiocManager::install('slingshot')")
}

if (file.exists("./temp/cluster1_obj_with_pseudotime.rds")) {
  cluster1_sub <- readRDS("./temp/cluster1_obj_with_pseudotime.rds")
  message("Loaded ./temp/cluster1_obj_with_pseudotime.rds")
} else if (file.exists("./temp/cluster1_obj_subclustered_scored.rds")) {
  cluster1_sub <- readRDS("./temp/cluster1_obj_subclustered_scored.rds")
  message("Loaded ./temp/cluster1_obj_subclustered_scored.rds")
} else if (file.exists("./temp/cluster1_obj_subclustered.rds")) {
  cluster1_sub <- readRDS("./temp/cluster1_obj_subclustered.rds")
  message("Loaded ./temp/cluster1_obj_subclustered.rds")
} else {
  cluster1_obj <- subset(cd8_tex_final, subset = seurat_clusters == "1")
  cluster1_obj$group <- droplevels(cluster1_obj$group)
  saveRDS(cluster1_obj, "./temp/cluster1_obj_final_clean.rds")

  cluster1_sub <- cluster1_obj
  DefaultAssay(cluster1_sub) <- "RNA"

  cluster1_sub <- FindVariableFeatures(cluster1_sub, nfeatures = 1500)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_after_HVG.rds")

  cluster1_sub <- ScaleData(cluster1_sub, verbose = FALSE)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_after_ScaleData.rds")

  cluster1_sub <- RunPCA(cluster1_sub, npcs = 20, verbose = FALSE)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_after_PCA.rds")

  cluster1_sub <- FindNeighbors(cluster1_sub, dims = 1:10, verbose = FALSE)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_after_Neighbors.rds")

  cluster1_sub <- FindClusters(cluster1_sub, resolution = 0.3, verbose = FALSE)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_after_Clusters.rds")

  cluster1_sub <- RunUMAP(cluster1_sub, dims = 1:10, verbose = FALSE)
  saveRDS(cluster1_sub, "./temp/cluster1_obj_subclustered.rds")
}

DefaultAssay(cluster1_sub) <- "RNA"
cluster1_sub$group <- factor(cluster1_sub$group)

Idents(cluster1_sub) <- "seurat_clusters"

if (!file.exists("./temp/cluster1_subcluster_markers.csv")) {
  cluster1_sub_markers <- FindAllMarkers(
    cluster1_sub,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.10,
    logfc.threshold = 0.10
  )
  write.csv(cluster1_sub_markers, "./temp/cluster1_subcluster_markers.csv", row.names = FALSE)
}

cluster1_sub_comp <- as.data.frame(table(group = cluster1_sub$group, subcluster = cluster1_sub$seurat_clusters)) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(prop = Freq / sum(Freq)) %>%
  dplyr::ungroup()
write.csv(cluster1_sub_comp, "./temp/cluster1_subcluster_composition_by_group.csv", row.names = FALSE)

p_cluster1_sub_umap <- DimPlot(
  cluster1_sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + theme_bw() + labs(title = "Cluster 1 internal subclusters")

ggsave("./figs/cluster1_subclusters_UMAP.pdf", p_cluster1_sub_umap, width = 7, height = 5)

need_scores <- !all(c("Cluster1TpexAxis1","Cluster1TexAxis1","Cluster1TerminalAxis1","Cluster1ActMemAxis1") %in%
                      colnames(cluster1_sub@meta.data))

if (need_scores) {
  cluster1_tpex_axis <- base::intersect(c("Tcf7","Slamf6","Il7r","Lef1","Bach2","Cxcr5"), rownames(cluster1_sub))
  cluster1_tex_axis <- base::intersect(c("Pdcd1","Tox","Nr4a2","Nr4a3","Ctla4","Havcr2","Tigit"), rownames(cluster1_sub))
  cluster1_terminal_axis <- base::intersect(c("Entpd1","Layn","Cd101","Cxcl13","Tnfrsf9"), rownames(cluster1_sub))
  cluster1_actmem_axis <- base::intersect(c("Gzma","S1pr1","Il7r","Cd69","Cxcr6","Itgb7","Lef1"), rownames(cluster1_sub))

  if (length(cluster1_tpex_axis) > 0) cluster1_sub <- AddModuleScore(cluster1_sub, features = list(cluster1_tpex_axis), name = "Cluster1TpexAxis")
  if (length(cluster1_tex_axis) > 0) cluster1_sub <- AddModuleScore(cluster1_sub, features = list(cluster1_tex_axis), name = "Cluster1TexAxis")
  if (length(cluster1_terminal_axis) > 0) cluster1_sub <- AddModuleScore(cluster1_sub, features = list(cluster1_terminal_axis), name = "Cluster1TerminalAxis")
  if (length(cluster1_actmem_axis) > 0) cluster1_sub <- AddModuleScore(cluster1_sub, features = list(cluster1_actmem_axis), name = "Cluster1ActMemAxis")

  saveRDS(cluster1_sub, "./temp/cluster1_obj_subclustered_scored.rds")
}

cluster1_sub_scores <- cluster1_sub@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    mean_tpex = if ("Cluster1TpexAxis1" %in% colnames(cluster1_sub@meta.data)) mean(Cluster1TpexAxis1, na.rm = TRUE) else NA_real_,
    mean_tex = if ("Cluster1TexAxis1" %in% colnames(cluster1_sub@meta.data)) mean(Cluster1TexAxis1, na.rm = TRUE) else NA_real_,
    mean_terminal = if ("Cluster1TerminalAxis1" %in% colnames(cluster1_sub@meta.data)) mean(Cluster1TerminalAxis1, na.rm = TRUE) else NA_real_,
    mean_actmem = if ("Cluster1ActMemAxis1" %in% colnames(cluster1_sub@meta.data)) mean(Cluster1ActMemAxis1, na.rm = TRUE) else NA_real_,
    n_cells = dplyr::n(),
    .groups = "drop"
  )

write.csv(cluster1_sub_scores, "./temp/cluster1_subcluster_axis_summary.csv", row.names = FALSE)

if (all(is.na(cluster1_sub_scores$mean_tpex))) {
  start_cluster <- as.character(cluster1_sub_scores$seurat_clusters[1])
} else {
  tmp <- cluster1_sub_scores %>%
    dplyr::filter(mean_tpex == max(mean_tpex, na.rm = TRUE)) %>%
    dplyr::arrange(mean_tex)
  start_cluster <- as.character(tmp$seurat_clusters[1])
}
write.csv(data.frame(start_cluster = start_cluster), "./temp/cluster1_pseudotime_start_cluster.csv", row.names = FALSE)

if (!("pseudotime" %in% colnames(cluster1_sub@meta.data))) {
  if (file.exists("./temp/cluster1_slingshot_sce.rds")) {
    sce <- readRDS("./temp/cluster1_slingshot_sce.rds")
    message("Loaded ./temp/cluster1_slingshot_sce.rds")
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(logcounts = as.matrix(GetAssayData(cluster1_sub, slot = "data")))
    )

    SingleCellExperiment::reducedDims(sce)$PCA <- Embeddings(cluster1_sub, "pca")
    SingleCellExperiment::reducedDims(sce)$UMAP <- Embeddings(cluster1_sub, "umap")
    SummarizedExperiment::colData(sce)$cluster <- as.character(cluster1_sub$seurat_clusters)
    SummarizedExperiment::colData(sce)$group <- as.character(cluster1_sub$group)


        sce <- slingshot::slingshot(
      sce,
      clusterLabels = "cluster",
      reducedDim = "PCA",
      start.clus = start_cluster
    )

    saveRDS(sce, "./temp/cluster1_slingshot_sce.rds")
  }

  pst <- slingshot::slingPseudotime(sce)
  pseudotime_vec <- if (is.matrix(pst)) pst[,1] else as.numeric(pst)

  cluster1_sub$pseudotime <- pseudotime_vec
  saveRDS(cluster1_sub, "./temp/cluster1_obj_with_pseudotime.rds")
}

pseudotime_df <- data.frame(
  cell = colnames(cluster1_sub),
  pseudotime = cluster1_sub$pseudotime,
  group = cluster1_sub$group,
  subcluster = cluster1_sub$seurat_clusters
)
write.csv(pseudotime_df, "./temp/cluster1_pseudotime_values.csv", row.names = FALSE)

umap_df <- as.data.frame(Embeddings(cluster1_sub, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$group <- cluster1_sub$group
umap_df$subcluster <- cluster1_sub$seurat_clusters
umap_df$pseudotime <- cluster1_sub$pseudotime
write.csv(umap_df, "./temp/cluster1_pseudotime_umap_df.csv", row.names = FALSE)

p_pseudotime_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.7) +
  theme_bw() +
  labs(title = "Cluster 1 pseudotime (UMAP)", color = "Pseudotime")
ggsave("./figs/cluster1_pseudotime_umap.pdf", p_pseudotime_umap, width = 7, height = 5)

p_pseudotime_group <- ggplot(pseudotime_df, aes(x = group, y = pseudotime)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  theme_bw() +
  labs(title = "Cluster 1 pseudotime by group", x = "Group", y = "Pseudotime")
ggsave("./figs/cluster1_pseudotime_by_group.pdf", p_pseudotime_group, width = 7, height = 5)

axis_time_df <- data.frame(
  pseudotime = cluster1_sub$pseudotime,
  group = cluster1_sub$group,
  subcluster = cluster1_sub$seurat_clusters
)
if ("Cluster1TexAxis1" %in% colnames(cluster1_sub@meta.data)) axis_time_df$TexAxis <- cluster1_sub$Cluster1TexAxis1
if ("Cluster1TpexAxis1" %in% colnames(cluster1_sub@meta.data)) axis_time_df$TpexAxis <- cluster1_sub$Cluster1TpexAxis1
if ("Cluster1TerminalAxis1" %in% colnames(cluster1_sub@meta.data)) axis_time_df$TerminalAxis <- cluster1_sub$Cluster1TerminalAxis1
if ("Cluster1ActMemAxis1" %in% colnames(cluster1_sub@meta.data)) axis_time_df$ActMemAxis <- cluster1_sub$Cluster1ActMemAxis1

write.csv(axis_time_df, "./temp/cluster1_pseudotime_axes_df.csv", row.names = FALSE)

axis_cols <- base::intersect(c("TexAxis","TpexAxis","TerminalAxis","ActMemAxis"), colnames(axis_time_df))

if (length(axis_cols) > 0) {
  axis_time_long <- axis_time_df %>%
    tidyr::pivot_longer(cols = all_of(axis_cols), names_to = "axis", values_to = "score")

  p_axis_time <- ggplot(axis_time_long, aes(x = pseudotime, y = score)) +
    geom_point(size = 0.4, alpha = 0.5) +
    geom_smooth(se = FALSE) +
    facet_wrap(~ axis, scales = "free_y", ncol = 2) +
    theme_bw() +
    labs(title = "Cluster 1 functional axes along pseudotime", x = "Pseudotime", y = "Score")
  ggsave("./figs/cluster1_pseudotime_axes.pdf", p_axis_time, width = 10, height = 6)

  p_axis_time_group <- ggplot(axis_time_long, aes(x = pseudotime, y = score, color = group)) +
    geom_point(size = 0.35, alpha = 0.5) +
    geom_smooth(se = FALSE) +
    facet_wrap(~ axis, scales = "free_y", ncol = 2) +
    theme_bw() +
    labs(title = "Cluster 1 functional axes along pseudotime by group", x = "Pseudotime", y = "Score")
  ggsave("./figs/cluster1_pseudotime_axes_by_group.pdf", p_axis_time_group, width = 10, height = 6)
}
# =========================================================
# 9.3b cluster1_sub DoRothEA activity linked to pseudotime
# 放置位置：
# 在 ggsave("./figs/cluster1_pseudotime_axes_by_group.pdf", ...) 之后
# 在 9.4 cluster 0 功能轴比较之前
# 即：第 2130 行后面插入
# =========================================================

# ---------- 0) 基础检查 ----------
if (!("dorothea" %in% names(cd8_final_obj@assays))) {
  stop("cd8_final_obj 中没有 dorothea assay。请先确认前面的 8.11f DoRothEA 已成功运行。")
}

# ---------- 1) 从带 dorothea assay 的 final object 中取出 cluster1 ----------
cluster1_doro <- subset(cd8_final_obj, subset = seurat_clusters == "1")

# 只保留和 cluster1_sub 共有的细胞，避免对象来源不一致
common_cells_c1 <- base::intersect(colnames(cluster1_doro), colnames(cluster1_sub))
if (length(common_cells_c1) == 0) {
  stop("cluster1_doro 和 cluster1_sub 没有共同细胞，无法关联 DoRothEA 与 pseudotime。")
}

cluster1_doro <- subset(cluster1_doro, cells = common_cells_c1)
cluster1_sub2 <- subset(cluster1_sub, cells = common_cells_c1)

# ---------- 2) 提取 dorothea TF activity ----------
tf_mat_c1 <- GetAssayData(cluster1_doro, assay = "dorothea", slot = "data")
tf_mat_c1 <- as.matrix(tf_mat_c1)

# 你最值得关注的 TF
tf_focus_c1 <- c(
  "Tcf7","Lef1","Bach2","Stat5a","Stat5b",
  "Tox","Nr4a1","Nr4a2","Nr4a3",
  "Tbx21","Eomes","Jun","Fos","Hif1a"
)
tf_focus_c1 <- tf_focus_c1[tf_focus_c1 %in% rownames(tf_mat_c1)]

if (length(tf_focus_c1) > 0) {
  tf_focus_df <- as.data.frame(t(tf_mat_c1[tf_focus_c1, , drop = FALSE]))
  tf_focus_df$cell <- rownames(tf_focus_df)

  meta_c1 <- data.frame(
    cell = colnames(cluster1_sub2),
    group = cluster1_sub2$group,
    subcluster = cluster1_sub2$seurat_clusters,
    pseudotime = cluster1_sub2$pseudotime
  )

  tf_focus_long <- dplyr::left_join(meta_c1, tf_focus_df, by = "cell") %>%
    tidyr::pivot_longer(
      cols = -c(cell, group, subcluster, pseudotime),
      names_to = "TF",
      values_to = "activity"
    )

  write.csv(tf_focus_long, "./temp/cluster1_sub_dorothea_TF_activity_with_pseudotime.csv", row.names = FALSE)

  # ---------- 3) 按 subcluster 比较 TF activity ----------
  p_tf_subcluster <- ggplot(tf_focus_long, aes(x = subcluster, y = activity, fill = subcluster)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ TF, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster1 subclusters: DoRothEA TF activity",
      x = "cluster1 subcluster",
      y = "TF activity"
    ) +
    theme(legend.position = "none")

  ggsave("./figs/cluster1_subcluster_dorothea_TF_activity.pdf", p_tf_subcluster, width = 12, height = 8)

  # ---------- 4) 按 group 比较 TF activity ----------
  p_tf_group <- ggplot(tf_focus_long, aes(x = group, y = activity)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ TF, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster1: DoRothEA TF activity by group",
      x = "Group",
      y = "TF activity"
    )

  ggsave("./figs/cluster1_group_dorothea_TF_activity.pdf", p_tf_group, width = 12, height = 8)

  # ---------- 5) TF activity 随 pseudotime 变化 ----------
  p_tf_pseudotime <- ggplot(tf_focus_long, aes(x = pseudotime, y = activity, color = group)) +
    geom_point(size = 0.35, alpha = 0.5) +
    geom_smooth(se = FALSE) +
    facet_wrap(~ TF, scales = "free_y", ncol = 3) +
    theme_bw() +
    labs(
      title = "Cluster1: DoRothEA TF activity along pseudotime by group",
      x = "Pseudotime",
      y = "TF activity"
    )

  ggsave("./figs/cluster1_dorothea_TF_activity_along_pseudotime_by_group.pdf", p_tf_pseudotime, width = 12, height = 8)

  # ---------- 6) 导出每个 subcluster 的平均 TF 活性 ----------
  avg_tf_subcluster <- tf_focus_long %>%
    dplyr::group_by(subcluster, TF) %>%
    dplyr::summarise(mean_activity = mean(activity, na.rm = TRUE), .groups = "drop")

  write.csv(avg_tf_subcluster, "./temp/cluster1_subcluster_dorothea_avgTF.csv", row.names = FALSE)
}
# -------------------------
# 9.4 cluster 0 功能轴比较
# -------------------------
if (file.exists("./temp/cluster0_obj_final_clean_scored.rds")) {
  cluster0_obj <- readRDS("./temp/cluster0_obj_final_clean_scored.rds")
  message("Loaded ./temp/cluster0_obj_final_clean_scored.rds")
} else if (file.exists("./temp/cluster0_obj_final_clean.rds")) {
  cluster0_obj <- readRDS("./temp/cluster0_obj_final_clean.rds")
  message("Loaded ./temp/cluster0_obj_final_clean.rds")
} else {
  cluster0_obj <- subset(cd8_tex_final, subset = seurat_clusters == "0")
  cluster0_obj$group <- droplevels(cluster0_obj$group)
  saveRDS(cluster0_obj, "./temp/cluster0_obj_final_clean.rds")
}

cluster0_obj$group <- factor(cluster0_obj$group)
DefaultAssay(cluster0_obj) <- "RNA"

genes_cluster0_refine <- present_only(cluster0_obj, c(
  "Gzma","Gzmk","Ccl5","Il7r","Lef1","S1pr1","Itgb7",
  "Cd69","Cxcr6",
  "Pdcd1","Tox","Nr4a3","Ctla4","Havcr2",
  "Prf1","Gzmb","Ifng","Xcl1"
))

avg_cluster0_by_group <- AverageExpression(
  cluster0_obj,
  group.by = "group",
  assays = "RNA",
  slot = "data"
)$RNA

avg_cluster0_by_group_sub <- avg_cluster0_by_group[
  rownames(avg_cluster0_by_group) %in% genes_cluster0_refine,
  ,
  drop = FALSE
]
write.csv(avg_cluster0_by_group_sub, "./temp/CD8_cluster0_by_group_actmem_avgexpr.csv")

p_dot_cluster0_group <- DotPlot(
  cluster0_obj,
  features = genes_cluster0_refine,
  group.by = "group"
) +
  RotatedAxis() +
  theme_bw() +
  labs(title = "Cluster 0 only: activated/memory markers by group")
ggsave("./figs/CD8_cluster0_by_group_ActMem_dotplot.pdf", p_dot_cluster0_group, width = 11, height = 5.5)

need_cluster0_scores <- !all(c("Cluster0ActMemAxis1","Cluster0ResidentAxis1","Cluster0TexAxis1","Cluster0EffectorAxis1") %in%
                              colnames(cluster0_obj@meta.data))

if (need_cluster0_scores) {
  cluster0_actmem_axis <- base::intersect(c("Gzma","Gzmk","Ccl5","Il7r","Lef1","S1pr1","Itgb7"), rownames(cluster0_obj))
  cluster0_resident_axis <- base::intersect(c("Cd69","Cxcr6","Itgb7"), rownames(cluster0_obj))
  cluster0_tex_axis <- base::intersect(c("Pdcd1","Tox","Nr4a3","Ctla4","Havcr2"), rownames(cluster0_obj))
  cluster0_effector_axis <- base::intersect(c("Prf1","Gzmb","Ifng","Xcl1"), rownames(cluster0_obj))

  if (length(cluster0_actmem_axis) > 0) cluster0_obj <- AddModuleScore(cluster0_obj, features = list(cluster0_actmem_axis), name = "Cluster0ActMemAxis")
  if (length(cluster0_resident_axis) > 0) cluster0_obj <- AddModuleScore(cluster0_obj, features = list(cluster0_resident_axis), name = "Cluster0ResidentAxis")
  if (length(cluster0_tex_axis) > 0) cluster0_obj <- AddModuleScore(cluster0_obj, features = list(cluster0_tex_axis), name = "Cluster0TexAxis")
  if (length(cluster0_effector_axis) > 0) cluster0_obj <- AddModuleScore(cluster0_obj, features = list(cluster0_effector_axis), name = "Cluster0EffectorAxis")

  saveRDS(cluster0_obj, "./temp/cluster0_obj_final_clean_scored.rds")
}

cluster0_score_cols <- base::intersect(
  c("Cluster0ActMemAxis1","Cluster0ResidentAxis1","Cluster0TexAxis1","Cluster0EffectorAxis1"),
  colnames(cluster0_obj@meta.data)
)

if (length(cluster0_score_cols) > 0) {
  cluster0_score_long <- cluster0_obj@meta.data %>%
    dplyr::select(group, all_of(cluster0_score_cols)) %>%
    tidyr::pivot_longer(cols = all_of(cluster0_score_cols), names_to = "axis", values_to = "score")

  write.csv(cluster0_score_long, "./temp/CD8_cluster0_function_axes_scores_long.csv", row.names = FALSE)

  p_cluster0_scores <- ggplot(cluster0_score_long, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ axis, scales = "free_y", ncol = 2) +
    theme_bw() +
    labs(title = "Cluster 0 functional axes by group", x = "Group", y = "Module score")
  ggsave("./figs/CD8_cluster0_function_axes_by_group.pdf", p_cluster0_scores, width = 10, height = 6)
}

cat("\n==== CD8 full workflow finished ====\n")
cat("核心表格：\n")
cat("./temp/CD8_quantity_and_props_point_estimate.csv\n")
cat("./temp/CD8_prop_in_immune_bootstrap_summary.csv\n")
cat("./temp/CD8_prop_in_T_bootstrap_summary.csv\n")
cat("./temp/CD8_final_clean_composition_labeled.csv\n")
cat("./temp/CD8_final_clean_composition_function_merged.csv\n")
cat("./temp/CD8_cluster1_tex_refine_avgexpr.csv\n")
cat("./temp/CD8_cluster1_vs_others_tex_refine_avgexpr.csv\n")
cat("./temp/CD8_cluster1_by_group_tex_refine_avgexpr.csv\n")
cat("./temp/CD8_cluster1_function_axes_scores_long.csv\n")
cat("./temp/CD8_cycling_by_group_avgexpr.csv\n")
cat("./temp/CD8_cycling_function_axes_scores_long.csv\n")
cat("./temp/cluster1_subcluster_axis_summary.csv\n")
cat("./temp/cluster1_pseudotime_values.csv\n")
cat("./temp/cluster1_pseudotime_axes_df.csv\n")
cat("./temp/CD8_cluster0_by_group_actmem_avgexpr.csv\n")
cat("./temp/CD8_cluster0_function_axes_scores_long.csv\n")

cat("\n核心图片：\n")
cat("./figs/CD8_final_clean_UMAP_annotated.pdf\n")
cat("./figs/CD8_final_clean_composition_labeled.pdf\n")
cat("./figs/CD8_final_clean_composition_function_merged.pdf\n")
cat("./figs/CD8_final_clean_marker_dotplot.pdf\n")
cat("./figs/CD8_final_clean_Tex_Tpex_dotplot.pdf\n")
cat("./figs/CD8_cluster1_by_group_Tex_Tpex_dotplot.pdf\n")
cat("./figs/CD8_cluster1_function_axes_by_group.pdf\n")
cat("./figs/CD8_cycling_by_group_dotplot.pdf\n")
cat("./figs/CD8_cycling_function_axes_by_group.pdf\n")
cat("./figs/cluster1_subclusters_UMAP.pdf\n")
cat("./figs/cluster1_pseudotime_umap.pdf\n")
cat("./figs/cluster1_pseudotime_by_group.pdf\n")
cat("./figs/cluster1_pseudotime_axes.pdf\n")
cat("./figs/cluster1_pseudotime_axes_by_group.pdf\n")
cat("./figs/CD8_cluster0_by_group_ActMem_dotplot.pdf\n")
cat("./figs/CD8_cluster0_function_axes_by_group.pdf\n")

# # --------------------- 7 关键质控：确认是否存在T细胞信号（Tex分析前必做） ---------------------
# combined <- readRDS("./temp/combined_after_marker_6.rds")
# all_markers <- readRDS("./temp/marker_analysis_6.rds")
# cat("7 开始【检查T细胞标志基因】\n")
# DefaultAssay(combined) <- "RNA"
#
# p_tcheck <- FeaturePlot(
#   combined,
#   features = c("Ptprc","Cd3d","Cd3e","Trac","Cd8a","Cd4","Nkg7","Ms4a1","Lyz2"),
#   reduction = "umap"
# )
#
# ggsave("./temp/QC_Tcell_markers_on_umap.pdf", plot = p_tcheck, width = 14, height = 8)
# cat("7 完成【检查T细胞标志基因】，已保存: ./temp/QC_Tcell_markers_on_umap.pdf\n")
#
# #--------------------- 7.1 Tex分析准备：提取T细胞并子集重聚类 ---------------------
# cat("7.1 开始【提取T细胞子集并重聚类】\n")
# DefaultAssay(combined) <- "RNA"
#
# # 1) 提取T细胞（用Trac/Cd3e/Cd3d门控；先宽松一点，后面再剔除NK等）
# t_cells <- subset(combined, subset = Trac > 0 | Cd3e > 0 | Cd3d > 0)
# cat("T cell subset cells:", ncol(t_cells), "\n")
#
# # 2) 在T细胞子集中重新走一遍标准流程（关键：不要沿用全免疫的图结构）
# t_cells <- NormalizeData(t_cells)
# t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2000)
# t_cells <- ScaleData(t_cells, features = VariableFeatures(t_cells), vars.to.regress = "percent.mt")
# t_cells <- RunPCA(t_cells, npcs = 30)
#
# # 建图/聚类/UMAP：Tex研究通常从 dims=1:20, k=20, res=0.6 起步
# t_cells <- FindNeighbors(t_cells, reduction = "pca", dims = 1:20, k.param = 20)
# t_cells <- FindClusters(t_cells, resolution = 0.6)
# t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:20)
#
# cat("7.1完成【T细胞子集重聚类】\n")
# # 把 t_cells 序列化保存到本地（RDS 格式最合适），后续分析直接 readRDS() 读回来即可，这样就不需要再次运行 6.6（除非你更改了门控/参数想重新生成）
# dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# saveRDS(t_cells, file = "./temp/t_cells_reclustered_7.1.rds")
# #读取
# t_cells <- readRDS("./temp/t_cells_reclustered_7.1.rds")
# cat("已读取 t_cells：", ncol(t_cells), "cells\n")

## # --------------------- 7.2. 分析Treg细胞 ---------------------
# cat ("7.2. 开始【分析Treg细胞】\n")
# # #============（7.2）分析耗竭细胞，同时也需要分析Treg细胞的=================
# # #确保对象与 assay
# # DefaultAssay(t_cells) <- "RNA"
# # Idents(t_cells) <- "seurat_clusters"
# # dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# # #Treg 分析（marker门控）
# # cat("（7.2） 开始【Treg分析】\n")
# # treg <- subset(
# #   t_cells,
# #   subset = (Foxp3 > 0 | Il2ra > 0) & (Trac > 0 | Cd3e > 0 | Cd3d > 0))
# # cat("Treg cells:", ncol(treg), "\n")
# # saveRDS(treg, "./temp/treg_subset_raw.rds")
#
# # treg <- readRDS("./temp/treg_subset_raw.rds")
# # cat("已读取 treg：", ncol(treg), "cells\n")
#
# # #Treg 子集重聚类（避免被 CD8/Tex 结构干扰）
# # cat("（7.2） 开始【Treg 子集重聚类】\n")
# # treg <- NormalizeData(treg)
# # treg <- FindVariableFeatures(treg, nfeatures = 2000)
# # treg <- ScaleData(treg, features = VariableFeatures(treg), vars.to.regress = "percent.mt")
# # treg <- RunPCA(treg, npcs = 30)
# # treg <- FindNeighbors(treg, dims = 1:20, k.param = 20)
# # treg <- FindClusters(treg, resolution = 0.4)
# # treg <- RunUMAP(treg, dims = 1:20)
# # p_treg_umap <- DimPlot(treg, reduction = "umap", label = TRUE) + ggtitle("Treg subset UMAP (reclustered)")
# # ggsave("./temp/Treg_subset_UMAP.pdf", p_treg_umap, width = 10, height = 8)
# # saveRDS(treg, "./temp/treg_subset_reclustered.rds")
#
#
# # treg <- readRDS("./temp/treg_subset_reclustered.rds")
# # # cat("已读取 treg：", ncol(treg), "cells\n")
#
# # # 先做“真 Treg 纯度表”（一眼定位混入的 cluster）
# # DefaultAssay(treg) <- "RNA"
# # Idents(treg) <- "seurat_clusters"
# #
# # # 关键：用 Foxp3 或 Ikzf2 定义“真 Treg”（比单用 Il2ra 稳）
# # treg$flag_trueTreg <- (FetchData(treg, "Foxp3")[,1] > 0) | (FetchData(treg, "Ikzf2")[,1] > 0)
# # treg$flag_cd25hi   <- (FetchData(treg, "Il2ra")[,1] > 0)
# #
# # # 每个 subcluster 的真Treg比例/Il2ra比例
# # purity_tbl <- data.frame(cluster = levels(Idents(treg)))
# # purity_tbl$n <- as.integer(table(Idents(treg))[purity_tbl$cluster])
# #
# # purity_tbl$frac_trueTreg <- tapply(treg$flag_trueTreg, Idents(treg), mean)[purity_tbl$cluster]
# # purity_tbl$frac_cd25hi   <- tapply(treg$flag_cd25hi,   Idents(treg), mean)[purity_tbl$cluster]
# #
# # purity_tbl <- purity_tbl[order(purity_tbl$frac_trueTreg), ]
# # print(purity_tbl)
# # write.csv(purity_tbl, "./temp/Treg_cluster_purity_table.csv", row.names = FALSE)
# #
# # # # Treg 子集 UMAP 上标记 Foxp3/Il2ra（可选）
# # # cat("（7.2） 开始【Treg 子集 UMAP 上标记 Foxp3/Il2ra】\n")
# # # for (g in intersect(c("Foxp3","Il2ra"), rownames(treg))) {
# # #   p <- FeaturePlot(treg, features = g) + ggtitle(paste("Treg subset:", g))
# # #   ggsave(paste0("./temp/Treg_", g, "_FeaturePlot.pdf"), p, width = 8, height = 6)
# # # }
# # #
# # # #  Treg 关键 marker 图（审稿人爱看的图）
# # # cat("（7.2） 开始【Treg 关键 marker 图】\n")
# # # treg_markers <- intersect(
# # #   c("Foxp3","Il2ra","Ikzf2","Ctla4","Tnfrsf18","Tnfrsf4","Tnfrsf9","Lag3","Pdcd1"),
# # #   rownames(treg)
# # # )
# # #
# # # p_treg_dot <- DotPlot(treg, features = treg_markers) + RotatedAxis() + ggtitle("Treg markers (DotPlot)")
# # # ggsave("./temp/Treg_markers_DotPlot.pdf", p_treg_dot, width = 12, height = 6)
# # #
# # # p_treg_fp <- FeaturePlot(treg, features = treg_markers, ncol = 5) + ggtitle("Treg markers (FeaturePlot)")
# # # ggsave("./temp/Treg_markers_FeaturePlot.pdf", p_treg_fp, width = 16, height = 9)
# # #
# # # #Treg 功能分数（抑制/激活）
# # # cat("（7.2） 开始【Treg 功能分数（抑制/激活）】\n")
# # # treg_sets <- list(
# # #   TregCore = intersect(c("Foxp3","Il2ra","Ikzf2"), rownames(treg)),
# # #   Suppressive = intersect(c("Ctla4","Il10","Lgals1","Tgfb1"), rownames(treg)),
# # #   Activation = intersect(c("Tnfrsf4","Tnfrsf9","Icos"), rownames(treg))
# # # )
# # # treg_sets <- treg_sets[sapply(treg_sets, length) >= 1]
# # #
# # # for (nm in names(treg_sets)) {
# # #   treg <- AddModuleScore(treg, features = list(treg_sets[[nm]]), name = paste0(nm, "_Score"))
# # # }
# # #
# # # treg_score_cols <- paste0(names(treg_sets), "_Score1")
# # # p_treg_scores <- VlnPlot(treg, features = treg_score_cols, pt.size = 0, ncol = 3) +
# # #   ggtitle("Treg module scores by Treg subclusters")
# # # ggsave("./temp/Treg_module_scores.pdf", p_treg_scores, width = 14, height = 7)
# # # cat("（7.2） 完成【Treg初步分析】\n")
# #
# # #按你这张 purity 表“选高纯度 cluster”（更贴合你当前结果）
# # purity <- read.csv("./temp/Treg_cluster_purity_table.csv")
# #
# # keep_clust_strict <- purity$cluster[purity$frac_trueTreg >= 0.80]   # 推荐
# # keep_clust_loose  <- purity$cluster[purity$frac_trueTreg >= 0.60]   # 可选探索
# #
# # Idents(treg) <- "seurat_clusters"
# #
# # treg_strict <- subset(treg, idents = as.character(keep_clust_strict))
# # cat("Strict cluster-based Treg cells:", ncol(treg_strict), "\n")
# # saveRDS(treg_strict, "./temp/treg_strict_cluster_based.rds")
#
# # treg_strict <- readRDS("./temp/treg_strict_cluster_based.rds")
# # # cat("已读取 treg_strict：", ncol(treg_strict), "cells\n")
# # # 拆成两个对象：真 Treg vs CD25hi 非Treg（激活 Teff）
# # #真 Treg（推荐用于论文的 Treg 主分析）
# # treg_strict <- subset(treg, subset = flag_trueTreg == TRUE)
# # cat("Strict Treg cells:", ncol(treg_strict), "\n")
# # saveRDS(treg_strict, "./temp/treg_strict.rds")
# # #CD25hi 但非 Treg（可作为“激活 CD4/Teff”单独分析）
# # cd25hi_nonTreg <- subset(treg, subset = flag_cd25hi == TRUE & flag_trueTreg == FALSE)
# # cat("CD25hi non-Treg cells:", ncol(cd25hi_nonTreg), "\n")
# # saveRDS(cd25hi_nonTreg, "./temp/cd25hi_nonTreg.rds")
# #
# # # 对“真 Treg”做发表级亚群注释（resting / activated / checkpoint-high / cycling）
# # DefaultAssay(treg_strict) <- "RNA"
# # Idents(treg_strict) <- "seurat_clusters"
# #
# # treg_strict <- NormalizeData(treg_strict)
# # treg_strict <- FindVariableFeatures(treg_strict, nfeatures = 2000)
# # treg_strict <- ScaleData(treg_strict, features = VariableFeatures(treg_strict), vars.to.regress = "percent.mt")
# # treg_strict <- RunPCA(treg_strict, npcs = 30)
# # treg_strict <- FindNeighbors(treg_strict, dims = 1:20, k.param = 20)
# # treg_strict <- FindClusters(treg_strict, resolution = 0.4)
# # treg_strict <- RunUMAP(treg_strict, dims = 1:20)
# #
# # # Cluster 3（孤岛）：当前面板不足以解释其独立性 → “Uncertain/Artifact?”，必须补 QC 与污染排查
# # treg_strict <- readRDS("./temp/treg_strict.rds")  # 或你保存的重聚类版本rds
# # Idents(treg_strict) <- "seurat_clusters"
# #
# # # 立刻做 QC（按 cluster）排查：是否低质量/高线粒体/低基因数
# # Idents(treg_strict) <- "seurat_clusters"
# # p_qc <- VlnPlot(
# #   treg_strict,
# #   features = c("nCount_RNA","nFeature_RNA","percent.mt"),
# #   pt.size = 0
# # ) + ggtitle("Strict Treg QC by cluster")
# # ggsave("./temp/Treg_strict_QC_by_cluster.pdf", p_qc, width = 12, height = 6)
# #
# # # 立刻做污染检查（DotPlot）：是否混入 NK / myeloid / doublet
# # genes_check <- base::intersect(
# #   c("Trac","Cd3e","Cd3d",          # T
# #     "Ncr1","Klrb1c","Tyrobp",      # NK-like / innate
# #     "LST1","S100a8","S100a9",      # myeloid
# #     "Col1a1","Pdgfra",             # fibro
# #     "Mki67","Top2a"),              # cycling
# #   rownames(treg_strict)
# # )
# #
# # p_dot_check <- DotPlot(treg_strict, features = genes_check) + RotatedAxis() +
# #   ggtitle("Strict Treg contamination/QC markers")
# # ggsave("./temp/Treg_strict_contamination_check_DotPlot.pdf", p_dot_check, width = 16, height = 6)
# # # Treg混有其他细胞验证结束：如果 cluster 3 显示 高 percent.mt / 低 nFeature 或出现明显 非T marker，建议你在论文主分析里将其标为 low-quality / contaminant 并剔除；如果没有这些问题，再考虑它是否为真实生物学状态
#
# # # 目的确定cluster 3 的 Tyrobp 确实代表髓系污染
# # # 首先固定 cluster 身份并指定要查的 cluster（这里是 “3”）
# # DefaultAssay(treg_strict) <- "RNA"
# # Idents(treg_strict) <- "seurat_clusters"   # 如果你锁定的是 Treg_res0.4，就改成 Idents(treg_strict) <- "Treg_res0.4"
# #
# # cl_target <- "3"
# #
# # # 然后，群体层证据：cluster 3 的 “T vs Myeloid” 面板（DotPlot + VlnPlot）
# # # DotPlot：一眼看 cluster 3 是否同时带 T 与髓系轴
# # genes_T <- c("Trac","Cd3e","Cd3d","Lck")
# # genes_Treg <- c("Foxp3","Ikzf2","Il2ra","Ctla4","Icos","Tnfrsf18")
# # genes_myeloid <- c("Tyrobp","Lyz2","Csf1r","LST1","Itgam","Fcgr3","S100a8","S100a9")
# #
# # panel <- base::intersect(c(genes_T, genes_Treg, genes_myeloid), rownames(treg_strict))
# #
# # p_dot <- DotPlot(treg_strict, features = panel, group.by = "seurat_clusters") +
# #   RotatedAxis() + ggtitle(paste0("Check cluster ", cl_target, ": T/Treg vs Myeloid markers"))
# # ggsave("./temp/Check_cluster3_Treg_vs_Myeloid_DotPlot.pdf", p_dot, width = 18, height = 7)
# #
# # # VlnPlot：看 cluster 3 内部髓系基因分布是否整体上移
# # vln_genes <- base::intersect(c("Trac","Foxp3","Il2ra","Ctla4","Tyrobp","Lyz2","Csf1r","Itgam"), rownames(treg_strict))
# #
# # p_vln <- VlnPlot(treg_strict, features = vln_genes, group.by = "seurat_clusters", pt.size = 0, ncol = 4) +
# #   ggtitle(paste0("Cluster ", cl_target, " distribution check (VlnPlot)"))
# # ggsave("./temp/Check_cluster3_Treg_vs_Myeloid_VlnPlot.pdf", p_vln, width = 18, height = 10)
# #
# # # 空间层证据：FeaturePlot 叠加（看是否“同一区域同时亮”）
# # fp_genes <- base::intersect(c("Trac","Foxp3","Il2ra","Ctla4","Tyrobp","Lyz2","Csf1r"), rownames(treg_strict))
# #
# # p_fp <- FeaturePlot(treg_strict, features = fp_genes, ncol = 4) +
# #   ggtitle(paste0("FeaturePlot: co-localization check for cluster ", cl_target))
# # ggsave("./temp/Check_cluster3_FeaturePlot_Treg_Myeloid.pdf", p_fp, width = 16, height = 10)
# #
# # # 最后，细胞层“硬证据”：cluster 3 内共表达比例（最关键）
# # # 取出 cluster 3 的细胞并计算共表达比例
# # cells3 <- WhichCells(treg_strict, idents = cl_target)
# #
# # df3 <- FetchData(
# #   treg_strict,
# #   vars = base::intersect(c("Trac","Foxp3","Il2ra","Ctla4","Tyrobp","Lyz2","Csf1r","LST1","Itgam"), rownames(treg_strict)),
# #   cells = cells3
# # )
# #
# # # --- 安全工具：如果基因列不存在，返回全FALSE；存在则返回 df[[gene]] > cutoff 的逻辑向量 ---
# # has_expr <- function(df, gene, cutoff = 0) {
# #   if (gene %in% colnames(df)) {
# #     return(df[[gene]] > cutoff)
# #   } else {
# #     return(rep(FALSE, nrow(df)))
# #   }
# # }
# #
# # # 1) T 细胞标记
# # flag_T <- has_expr(df3, "Trac", 0)
# #
# # # 2) Treg 标记（向量级，用 | 而不是 ||；用 has_expr 避免缺列）
# # flag_Treg <- has_expr(df3, "Foxp3", 0) | has_expr(df3, "Ikzf2", 0) | has_expr(df3, "Ctla4", 0)
# #
# # # 3) 髓系标记（同理）
# # flag_myeloid <- has_expr(df3, "Tyrobp", 0) |
# #   has_expr(df3, "Lyz2", 0) |
# #   has_expr(df3, "Csf1r", 0) |
# #   has_expr(df3, "LST1", 0) |
# #   has_expr(df3, "Itgam", 0)
# #
# # cat("Cluster 3 cells:", nrow(df3), "\n")
# # cat("T (Trac>0) fraction:", mean(flag_T), "\n")
# # cat("Treg (Foxp3/Ikzf2/Ctla4>0) fraction:", mean(flag_Treg), "\n")
# # cat("Myeloid(any myeloid marker >0) fraction:", mean(flag_myeloid), "\n")
# # cat("T & Myeloid coexpression fraction:", mean(flag_T & flag_myeloid), "\n")
# # cat("Treg & Myeloid coexpression fraction:", mean(flag_Treg & flag_myeloid), "\n")
# #
# # # “铁证如山”：加一条髓系强证据（Lyz2/Csf1r/Lst1/Itgam 多个同时阳性）
# # # 至少2个髓系基因同时 >0 的细胞，视为“强髓系”
# # my_genes <- c("Lyz2","Csf1r","LST1","Itgam","S100a8","S100a9")
# # present <- my_genes[my_genes %in% colnames(df3)]
# # my_mat <- as.matrix(df3[, present, drop = FALSE] > 0)
# # flag_myeloid_strong <- rowSums(my_mat) >= 2
# #
# # cat("Strong myeloid fraction (>=2 myeloid genes >0):", mean(flag_myeloid_strong), "\n")
# # cat("T & strong-myeloid coexpression fraction:", mean(flag_T & flag_myeloid_strong), "\n")
# # # 目的确定cluster 3 的 Tyrobp 确实代表髓系污染结束
#
# # ============================================================
# # Strict Treg 清理（针对 cluster 3 的 myeloid/doublet）→ 重聚类 → 输出发表用图
# # 你只需要改 3 个地方：
# #   (1) 读入对象路径（treg_strict_rds）
# #   (2) 你要用的聚类列名（cluster_col，默认 seurat_clusters）
# #   (3) 目标可疑簇编号（target_cluster，默认 "3"）
# # ============================================================
# #
# # # ---------- 0. 读入 / 基本设置 ----------
# # treg_strict_rds <- "./temp/treg_strict.rds"   # <- 按你的实际路径改
# # treg_strict <- readRDS(treg_strict_rds)
# #
# # DefaultAssay(treg_strict) <- "RNA"
# #
# # cluster_col <- "seurat_clusters"   # 如果你锁定过例如 "Treg_res0.4"，就改成那个
# # target_cluster <- "3"
# #
# # # 强制把 Identity 设为你要评估的聚类列
# # Idents(treg_strict) <- cluster_col
# #
# # dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# #
# # # # ====== 全对象细胞级剔除：T & Myeloid 共表达（推荐） ======
# # # DefaultAssay(treg_strict) <- "RNA"
# # #
# # # genes_need <- c(
# # #   # T
# # #   "Trac","Cd3e","Cd3d",
# # #   # myeloid / neutrophil / macrophage
# # #   "Tyrobp","Lyz2","Csf1r","LST1","Itgam","Fcgr3","S100a8","S100a9","C1qa","H2-Ab1"
# # # )
# # # genes_need <- base::intersect(genes_need, rownames(treg_strict))
# # #
# # # df_all <- FetchData(treg_strict, vars = genes_need)
# # #
# # # has_expr <- function(df, gene, cutoff = 0) {
# # #   if (gene %in% colnames(df)) df[[gene]] > cutoff else rep(FALSE, nrow(df))
# # # }
# # #
# # # flag_T <- has_expr(df_all, "Trac", 0)
# # #
# # # flag_myeloid <- has_expr(df_all, "Tyrobp", 0) |
# # #   has_expr(df_all, "Lyz2", 0) |
# # #   has_expr(df_all, "Csf1r", 0) |
# # #   has_expr(df_all, "LST1", 0) |
# # #   has_expr(df_all, "Itgam", 0) |
# # #   has_expr(df_all, "Fcgr3", 0) |
# # #   has_expr(df_all, "S100a8", 0) |
# # #   has_expr(df_all, "S100a9", 0) |
# # #   has_expr(df_all, "C1qa", 0) |
# # #   has_expr(df_all, "H2-Ab1", 0)
# # #
# # # bad_cells <- rownames(df_all)[flag_T & flag_myeloid]
# # #
# # # cat("All cells:", nrow(df_all), "\n")
# # # cat("Bad (T & Myeloid) cells:", length(bad_cells), "\n")
# # # cat("Bad fraction:", length(bad_cells) / nrow(df_all), "\n")
# # #
# # # treg_clean <- subset(treg_strict, cells = base::setdiff(colnames(treg_strict), bad_cells))
# # # cat("Before clean:", ncol(treg_strict), "After clean:", ncol(treg_clean), "\n")
# # #
# # # saveRDS(treg_clean, "./temp/treg_strict_clean_removed_T_myeloid_doublets_ALL.rds")
# #
# # # 结果：Trac > 0 且 任意 1 个髓系/APC marker > 0 → 就删（因此你删了 38%） 因此需要重新调整清洗代码
# #
# # # ====== 新版：方案A（T & >=2 个髓系/APC基因阳性才删）======
# # DefaultAssay(treg_strict) <- "RNA"
# #
# # genes_my <- c("Tyrobp","Lyz2","Csf1r","LST1","Itgam","Fcgr3","S100a8","S100a9","C1qa","H2-Ab1")
# # genes_my <- base::intersect(genes_my, rownames(treg_strict))
# #
# # if (length(genes_my) < 2) {
# #   stop("髓系基因在对象中可用数量不足（<2），无法按方案A判定。")
# # }
# #
# # df_all <- FetchData(treg_strict, vars = base::intersect(c("Trac", genes_my), rownames(treg_strict)))
# #
# # flag_T <- df_all$Trac > 0
# #
# # my_mat <- as.matrix(df_all[, genes_my, drop = FALSE] > 0)
# # flag_my_strong <- rowSums(my_mat) >= 2
# #
# # bad_cells <- rownames(df_all)[flag_T & flag_my_strong]
# #
# # cat("All cells:", nrow(df_all), "\n")
# # cat("Bad (T & >=2 myeloid/APC genes) cells:", length(bad_cells), "\n")
# # cat("Bad fraction:", length(bad_cells) / nrow(df_all), "\n")
# #
# # treg_clean <- subset(treg_strict, cells = base::setdiff(colnames(treg_strict), bad_cells))
# # cat("Before clean:", ncol(treg_strict), "After clean:", ncol(treg_clean), "\n")
# # saveRDS(treg_clean, "./temp/treg_strict_clean_removed_T_myeloid_doublets_ALL.rds")
# #
# # #
# # # # 发现cluster6有髓系细胞污染：先量化 cluster 6 内到底有多少是 “T & Myeloid 共表达”
# # # DefaultAssay(treg_clean) <- "RNA"
# # # Idents(treg_clean) <- "seurat_clusters"   # 或你 clean 后锁定的列名
# # #
# # # cl <- "6"
# # # cells6 <- WhichCells(treg_clean, idents = cl)
# # #
# # # genes_check <- base::intersect(
# # #   c("Trac","Foxp3","Il2ra","Ctla4","Tyrobp","Lyz2","Csf1r","Itgam","H2-Ab1"),
# # #   rownames(treg_clean)
# # # )
# # #
# # # df6 <- FetchData(treg_clean, vars = genes_check, cells = cells6)
# # #
# # # has_expr <- function(df, gene, cutoff = 0) {
# # #   if (gene %in% colnames(df)) df[[gene]] > cutoff else rep(FALSE, nrow(df))
# # # }
# # #
# # # flag_T <- has_expr(df6, "Trac", 0)
# # # flag_my <- has_expr(df6, "Tyrobp", 0) | has_expr(df6, "Lyz2", 0) | has_expr(df6, "Csf1r", 0) |
# # #            has_expr(df6, "Itgam", 0)  | has_expr(df6, "H2-Ab1", 0)
# # #
# # # cat("Cluster 6 cells:", nrow(df6), "\n")
# # # cat("T fraction:", mean(flag_T), "\n")
# # # cat("Myeloid/APC fraction:", mean(flag_my), "\n")
# # # cat("T & Myeloid/APC coexpression fraction:", mean(flag_T & flag_my), "\n")
# # #
# # # # 发现cluster6有髓系细胞污染：量化 cluster 6 内到底有多少是 “T & Myeloid 共表达”结束
# # # # 判定标准：T & Myeloid/APC coexpression fraction 如果 ≥ 0.2：非常像 doublet/污染，建议剔除该类细胞或整簇剔除
# # #
# # # # 只删 cluster 6 内 “Trac+ 且 myeloid/APC+” 的细胞，再重聚类
# # # bad6 <- rownames(df6)[flag_T & flag_my]
# # #
# # # treg_clean2 <- subset(treg_clean, cells = base::setdiff(colnames(treg_clean), bad6))
# # # saveRDS(treg_clean2, "./temp/treg_clean2_removed_cluster6_T_myeloid_cells.rds")
# # # 只删 cluster 6 内 “Trac+ 且 myeloid/APC+” 的细胞结束
# #
# #
# #
# # # ---------- 2. 清理后：重新跑 PCA/聚类/UMAP（重要：图结构会变化） ----------
# # DefaultAssay(treg_clean) <- "RNA"
# #
# # # 若你之前已经 Normalize/Scale 过，也可以直接复用；但为了“最稳可复现”，这里重新跑一遍
# # treg_clean <- NormalizeData(treg_clean)
# # treg_clean <- FindVariableFeatures(treg_clean, selection.method = "vst", nfeatures = 2000)
# # treg_clean <- ScaleData(treg_clean, features = VariableFeatures(treg_clean), vars.to.regress = "percent.mt")
# # treg_clean <- RunPCA(treg_clean, npcs = 30)
# #
# # # 发表常用起步：dims 1:20, k=20, res=0.4（你可按需要改）
# # dims_use <- 1:20
# # k_use <- 20
# # res_use <- 0.4
# #
# # treg_clean <- FindNeighbors(treg_clean, reduction = "pca", dims = dims_use, k.param = k_use)
# # treg_clean <- FindClusters(treg_clean, resolution = res_use)
# # treg_clean <- RunUMAP(treg_clean, reduction = "pca", dims = dims_use)
# #
# # # 锁定聚类列（防止后面被覆盖）
# # clean_cluster_col <- paste0("Treg_clean_res", res_use)
# # treg_clean[[clean_cluster_col]] <- Idents(treg_clean)
# # Idents(treg_clean) <- clean_cluster_col
# #
# # saveRDS(treg_clean, "./temp/treg_strict_clean_reclustered.rds")
# #
# # # ---------- 3. 输出“审稿人友好”的 QC + 污染检查 ----------
# # # 3.1 QC violin
# # p_qc <- VlnPlot(
# #   treg_clean,
# #   features = c("nCount_RNA","nFeature_RNA","percent.mt"),
# #   group.by = clean_cluster_col,
# #   pt.size = 0,
# #   ncol = 3
# # ) + ggtitle("Strict Treg (clean) QC by cluster")
# # ggsave("./temp/Treg_clean_QC_by_cluster.pdf", p_qc, width = 16, height = 5)
# #
# # # 3.2 污染检查 DotPlot
# # markers_contam <- c("Trac","Cd3e","Cd3d","Ncr1","Klrb1c","Tyrobp","Lyz2","Csf1r","LST1","Itgam","S100a8","S100a9","C1qa","H2-Ab1","Mki67","Top2a")
# # markers_contam <- base::intersect(markers_contam, rownames(treg_clean))
# #
# # p_contam <- DotPlot(treg_clean, features = markers_contam, group.by = clean_cluster_col) +
# #   RotatedAxis() +
# #   ggtitle("Strict Treg (clean) contamination/QC markers")
# # ggsave("./temp/Treg_clean_contamination_check_DotPlot.pdf", p_contam, width = 18, height = 7)
# # # QC + 污染检查 结束
# #
# # # # 生成“每个 cluster 的 Treg core 纯度”和“髓系污染率”表
# # # DefaultAssay(treg_clean) <- "RNA"
# # # Idents(treg_clean) <- "seurat_clusters"  # 或你锁定的 clean_cluster_col
# # #
# # # # 需要的基因（存在则用，不存在自动忽略）
# # # genes <- c("Trac","Foxp3","Ikzf2","Il2ra","Ctla4",
# # #            "Tyrobp","Lyz2","Csf1r","Itgam","S100a8","S100a9","C1qa","H2-Ab1")
# # # genes <- base::intersect(genes, rownames(treg_clean))
# # #
# # # df <- FetchData(treg_clean, vars = genes)
# # # df$cluster <- as.character(Idents(treg_clean))
# # #
# # # has <- function(g) if (g %in% colnames(df)) df[[g]] > 0 else rep(FALSE, nrow(df))
# # #
# # # flag_T      <- has("Trac")
# # # flag_Treg   <- has("Foxp3") | has("Ikzf2") | has("Il2ra") | has("Ctla4")
# # #
# # # # “强髓系”：≥2 个髓系/APC 基因阳性（与你方案A一致）
# # # my_genes <- base::intersect(c("Tyrobp","Lyz2","Csf1r","Itgam","C1qa","H2-Ab1","S100a8","S100a9"), colnames(df))
# # # my_mat <- as.matrix(df[, my_genes, drop = FALSE] > 0)
# # # flag_my_strong <- rowSums(my_mat) >= 2
# # #
# # # tab <- aggregate(
# # #   cbind(T_frac = flag_T, Treg_frac = flag_Treg, MyStrong_frac = flag_my_strong) ~ cluster,
# # #   data = df,
# # #   FUN = mean
# # # )
# # # tab$n <- as.integer(table(df$cluster)[tab$cluster])
# # # tab <- tab[order(-tab$Treg_frac), ]
# # # print(tab)
# # #
# # # write.csv(tab, "./temp/Treg_clean_cluster_purity_table.csv", row.names = FALSE)
# # # # 生成“每个 cluster 的 Treg core 纯度”和“髓系污染率”表结束
# # # 结果这个Table中Treg的比例都是1，说明前面Treg的条件太宽松。不利于后面Treg的分类。
# #
# # #“真正发表级”的纯度表代码（Foxp3/Il2ra 双阳 + MyStrong）
# # DefaultAssay(treg_clean) <- "RNA"
# # Idents(treg_clean) <- "seurat_clusters"   # 或你的最终聚类列名
# #
# # genes <- c("Trac","Foxp3","Il2ra","Ikzf2","Ctla4",
# #            "Tyrobp","Lyz2","Csf1r","LST1","Itgam","Fcgr3","S100a8","S100a9","C1qa","H2-Ab1")
# # genes <- base::intersect(genes, rownames(treg_clean))
# #
# # df <- FetchData(treg_clean, vars = genes)
# # df$cluster <- as.character(Idents(treg_clean))
# #
# # has <- function(g) if (g %in% colnames(df)) df[[g]] > 0 else rep(FALSE, nrow(df))
# #
# # flag_T <- has("Trac")
# # flag_Foxp3 <- has("Foxp3")
# # flag_Il2ra <- has("Il2ra")
# # flag_Foxp3Il2ra <- flag_Foxp3 & flag_Il2ra
# #
# # my_genes <- base::intersect(c("Tyrobp","Lyz2","Csf1r","LST1","Itgam","Fcgr3","S100a8","S100a9","C1qa","H2-Ab1"), colnames(df))
# # my_mat <- as.matrix(df[, my_genes, drop = FALSE] > 0)
# # flag_my_strong <- rowSums(my_mat) >= 2
# #
# # tab <- aggregate(
# #   cbind(T_frac = flag_T,
# #         Foxp3_frac = flag_Foxp3,
# #         Il2ra_frac = flag_Il2ra,
# #         Foxp3_Il2ra_frac = flag_Foxp3Il2ra,
# #         MyStrong_frac = flag_my_strong) ~ cluster,
# #   data = df, FUN = mean
# # )
# # tab$n <- as.integer(table(df$cluster)[tab$cluster])
# # tab <- tab[order(-tab$Foxp3_Il2ra_frac), ]
# # print(tab)
# #
# # write.csv(tab, "./temp/Treg_clean_cluster_purity_table_publishable.csv", row.names = FALSE)
# # #“真正发表级”的纯度表代码（Foxp3/Il2ra 双阳 + MyStrong）结束。
# #
# # # 给每个 cluster 打标签（HighConf / Review / NonTreg / Contaminated）
# # Idents(treg_clean) <- "seurat_clusters"
# #
# # treg_clean$Treg_QC_flag <- "Review"
# #
# # # 高置信Treg：0,3（你表里Foxp3&Il2ra最高且MyStrong≈0）
# # treg_clean$Treg_QC_flag[Idents(treg_clean) %in% c("0","3")] <- "HighConf_Treg"
# #
# # # 可能Treg（需要Ikzf2/Ctla4等辅证）
# # treg_clean$Treg_QC_flag[Idents(treg_clean) %in% c("4")] <- "Probable_Treg"
# #
# # # Il2ra高但Foxp3极低：更像activated conventional T
# # treg_clean$Treg_QC_flag[Idents(treg_clean) %in% c("1","2")] <- "NonTreg_Activated_T"
# #
# # # 混杂/污染
# # treg_clean$Treg_QC_flag[Idents(treg_clean) %in% c("5")] <- "Mixed_Review"
# # treg_clean$Treg_QC_flag[Idents(treg_clean) %in% c("6")] <- "Contaminated"
# #
# # p <- DimPlot(treg_clean, group.by = "Treg_QC_flag", label = TRUE) +
# #   ggtitle("Treg QC flags (publishable purity table)")
# # ggsave("./temp/Treg_QC_flag_UMAP.pdf", p, width = 11, height = 8)
# #
# # # 给每个 cluster 打标签（HighConf / Review / NonTreg / Contaminated）完成
# #
# # # ====== 生成 Treg 主分析对象（先做一个可发表版本） ======
# # Idents(treg_clean) <- "seurat_clusters"
# # treg_main <- subset(treg_clean, idents = c("0","3","4"))  # 4 先保留，后面用 marker 再确认
# # saveRDS(treg_main, "./temp/treg_main_for_paper.rds")
# # # 生成 Treg 主分析对象（先做一个可发表版本）结束
# #
# # # ====== 用“审稿人友好”的 marker 复核 cluster 4 ======
# # DefaultAssay(treg_clean) <- "RNA"
# # Idents(treg_clean) <- "seurat_clusters"
# #
# # markers_check <- base::intersect(
# #   c(
# #     # T/Treg core
# #     "Trac","Cd3e","Foxp3","Il2ra","Ikzf2","Ctla4",
# #     # activated/eTreg
# #     "Icos","Tnfrsf18","Tnfrsf4","Tnfrsf9",
# #     # naive/conv T
# #     "Tcf7","Lef1","Il7r",
# #     # cycling
# #     "Mki67","Top2a"
# #   ),
# #   rownames(treg_clean)
# # )
# #
# # p_dot <- DotPlot(treg_clean, features = markers_check, group.by = "seurat_clusters") + RotatedAxis() +
# #   ggtitle("Check cluster 4: Treg vs activated T markers")
# # ggsave("./temp/Check_cluster4_markers_DotPlot.pdf", p_dot, width = 12, height = 6)
# #
# # p_vln <- VlnPlot(treg_clean, features = c("Foxp3","Il2ra","Ikzf2","Ctla4","Icos","Tcf7","Il7r"),
# #                  group.by = "seurat_clusters", pt.size = 0)
# # ggsave("./temp/Check_cluster4_markers_VlnPlot.pdf", p_vln, width = 12, height = 6)
# #
# # p_fp <- FeaturePlot(treg_clean, features = c("Foxp3","Il2ra","Ikzf2","Ctla4","Icos","Tcf7","Il7r"),
# #                     reduction = "umap")
# # ggsave("./temp/Check_cluster4_markers_FeaturePlot.pdf", p_fp, width = 12, height = 8)
# # # 复核 cluster 4结束
# #
# # # 对 cluster 4 做一次细胞级联合门控并重新算“高置信 Treg”比例
# # # 放在：你已完成 treg_clean 的聚类/UMAP + 上述检查图之后
# # DefaultAssay(treg_clean) <- "RNA"
# # Idents(treg_clean) <- "seurat_clusters"
# #
# # cells4 <- WhichCells(treg_clean, idents = "4")
# # df4 <- FetchData(treg_clean, vars = base::intersect(c("Foxp3","Il2ra","Ikzf2","Il7r"), rownames(treg_clean)), cells = cells4)
# #
# # # 高置信 Treg（>=2 of 3 core）
# # core_mat <- as.matrix(df4[, base::intersect(c("Foxp3","Il2ra","Ikzf2"), colnames(df4)), drop = FALSE] > 0)
# # flag_Treg_strict <- rowSums(core_mat) >= 2
# #
# # # 可疑的 conventional/naive-like（Il7r 高但 Foxp3 低）
# # flag_conv_like <- ("Il7r" %in% colnames(df4)) & (df4$Il7r > 0) & (("Foxp3" %in% colnames(df4)) & (df4$Foxp3 == 0))
# #
# # cat("Cluster4 n =", nrow(df4), "\n")
# # cat("Cluster4 strict-Treg fraction (>=2/3):", mean(flag_Treg_strict), "\n")
# # cat("Cluster4 Il7r+ & Foxp3- fraction:", mean(flag_conv_like), "\n")
# # # 对 cluster 4 做一次细胞级联合门控并重新算“高置信 Treg”比例结束
# # #结果判读：如果 strict-Treg fraction 很高（比如 >0.6），cluster 4 可以非常放心地作为 eTreg/activated Treg 写进主结果。
# # #		如果 Il7r+Foxp3- 比例很高（比如 >0.3），则建议对 cluster 4 做一次细胞级剔除再重聚类（我也可以给你最简剔除代码）
# #
# # # ======================================================================
# # # Step 1: Treg 子群“发表级命名”（rTreg / eTreg / cycling / checkpoint-high）
# # # 目的：
# # #   1) 在 treg_main（高置信Treg主对象）上输出：UMAP + DotPlot + FeaturePlot + VlnPlot
# # #   2) 计算并保存 module score（亚群程序分数）
# # #   3) 输出“每个cluster的程序分数统计表”，便于命名
# # #   4) 保存中间对象，后续改代码可直接 readRDS 快速复跑
# # # ======================================================================
# #
# # cat("==== Step 1: Treg subclass naming START ====\n")
# #
# # # ---------- 0. 路径与基础设置 ----------
# # dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# #
# # # 强制用 base 的集合函数，避免 conflicted 报错
# # # （如果你已经 conflicts_prefer 过，也没问题）
# # base_intersect <- base::intersect
# # base_setdiff   <- base::setdiff
# # base_unname    <- base::unname
# #
# # # ---------- 1. 从 treg_clean 生成 treg_main（论文主分析对象）并缓存 ----------
# # # 【若你前面已经生成过 treg_main，可直接跳过本段：读取缓存即可】
# # treg_main_rds <- "./temp/treg_main_for_paper.rds"
# #
# # if (!file.exists(treg_main_rds)) {
# #   cat(">> No cached treg_main. Creating and saving...\n")
# #   Idents(treg_clean) <- "seurat_clusters"
# #   treg_main <- subset(treg_clean, idents = c("0","3","4"))
# #   saveRDS(treg_main, treg_main_rds)
# # } else {
# #   cat(">> Found cached treg_main. Loading...\n")
# #   treg_main <- readRDS(treg_main_rds)
# # }
# #
# # cat("treg_main cells:", ncol(treg_main), "\n")
# # cat("treg_main clusters:\n")
# # print(table(Idents(treg_main)))
# #
# # # ---------- 2. （可选）cluster4 轻度清理后的版本：如果你做过并保存了，就优先读取 ----------
# # # 你若在前面已经做了 cluster4 清理并保存为 treg_clean2_...，在这里可切换使用
# # # 默认不自动切换，以免你混用对象口径；如你确实想用清理后重聚类版本，把 use_clean2 <- TRUE
# # use_clean2 <- FALSE
# # treg_clean2_rds <- "./temp/treg_clean2_removed_cluster4_Il7rpos_Foxp3neg_cells.rds"
# #
# # if (use_clean2 && file.exists(treg_clean2_rds)) {
# #   cat(">> Using treg_clean2 (cluster4 cleaned) as the base.\n")
# #   treg_clean2 <- readRDS(treg_clean2_rds)
# #
# #   # 注意：清理后建议重新定义主分析对象（仍以高置信Treg为核心）
# #   Idents(treg_clean2) <- "seurat_clusters"
# #   # 这里不写死簇编号（因为清理后簇编号可能变），建议你先看 table(Idents)
# #   # 暂时给一个最保守做法：仍然取你确认的簇（若编号变化，你需改这里）
# #   treg_main <- subset(treg_clean2, idents = c("0","3","4"))
# #   saveRDS(treg_main, "./temp/treg_main_for_paper_from_clean2.rds")
# #   cat("treg_main (from clean2) cells:", ncol(treg_main), "\n")
# # }
# #
# # # ---------- 3. Step1 图：基础 UMAP（cluster） ----------
# # p_umap_cluster <- DimPlot(treg_main, reduction = "umap", label = TRUE) +
# #   ggtitle("Treg main - clusters")
# # ggsave("./temp/Step1_TregMain_UMAP_clusters.pdf", p_umap_cluster, width = 10, height = 8)
# #
# # # ---------- 4. Step1 marker 面板（rTreg/eTreg/checkpoint/cycling） ----------
# # DefaultAssay(treg_main) <- "RNA"
# # Idents(treg_main) <- "seurat_clusters"
# #
# # markers_subclass <- c(
# #   # Treg core / identity
# #   "Trac","Cd3e","Foxp3","Il2ra","Ikzf2","Ctla4",
# #   # activated/eTreg
# #   "Icos","Tnfrsf18","Tnfrsf4","Tnfrsf9",
# #   # checkpoint-high（Treg也会高，写作要用“checkpoint-high Treg”口径）
# #   "Pdcd1","Lag3","Tigit",
# #   # resting/naive-like（在Treg里通常较低，但可区分rTreg-like）
# #   "Tcf7","Lef1","Il7r",
# #   # cycling
# #   "Mki67","Top2a"
# # )
# # markers_subclass <- base_intersect(markers_subclass, rownames(treg_main))
# #
# # p_dot <- DotPlot(treg_main, features = markers_subclass) + RotatedAxis() +
# #   ggtitle("Treg subclass marker panel (cluster-level)")
# # ggsave("./temp/Step1_TregMain_subclass_markers_DotPlot.pdf", p_dot, width = 14, height = 6)
# #
# # # FeaturePlot（空间共定位，用于主图或补图）
# # fp_genes <- base_intersect(c("Foxp3","Il2ra","Icos","Ctla4","Tcf7","Mki67","Pdcd1"), rownames(treg_main))
# # p_fp <- FeaturePlot(treg_main, features = fp_genes, reduction = "umap")
# # ggsave("./temp/Step1_TregMain_subclass_FeaturePlot.pdf", p_fp, width = 12, height = 8)
# #
# # # VlnPlot（看簇内分布，判断是否混杂）
# # vln_genes <- base_intersect(c("Foxp3","Il2ra","Ikzf2","Ctla4","Icos","Tcf7","Il7r","Mki67","Pdcd1"), rownames(treg_main))
# # p_vln <- VlnPlot(treg_main, features = vln_genes, group.by = "seurat_clusters", pt.size = 0)
# # ggsave("./temp/Step1_TregMain_subclass_VlnPlot.pdf", p_vln, width = 14, height = 7)
# #
# # # ---------- 5. Step1 模块分数（module scores）并缓存对象 ----------
# # # 目的：用程序分数辅助命名：rTreg/eTreg/checkpoint/cycling
# # treg_main_scored_rds <- "./temp/treg_main_scored_step1.rds"
# #
# # if (!file.exists(treg_main_scored_rds)) {
# #   cat(">> No cached scored treg_main. Computing module scores...\n")
# #
# #   # 选基因集：尽量“审稿人不太好挑刺”，但也要避免过多基因缺失
# #   genes_rTreg <- c("Tcf7","Lef1","Il7r","Ccr7","Ltb")
# #   genes_eTreg <- c("Icos","Ctla4","Tnfrsf18","Tnfrsf4","Tnfrsf9","Ikzf2","Batf")
# #   genes_checkpoint <- c("Pdcd1","Lag3","Tigit","Havcr2")
# #   genes_cycling <- c("Mki67","Top2a","Stmn1","Hmgb2")
# #
# #   # 与数据基因取交集
# #   score_list <- list(
# #     rTreg = base_intersect(genes_rTreg, rownames(treg_main)),
# #     eTreg = base_intersect(genes_eTreg, rownames(treg_main)),
# #     checkpoint = base_intersect(genes_checkpoint, rownames(treg_main)),
# #     cycling = base_intersect(genes_cycling, rownames(treg_main))
# #   )
# #
# #   # 打印每个分数实际使用了多少基因（便于你判断是否太少）
# #   cat("Module score gene counts:\n")
# #   print(sapply(score_list, length))
# #
# #   # 逐个生成分数字段（列名可控，避免 AddModuleScore 自动命名带来的混乱）
# #   # 这里使用一个安全循环：每次 AddModuleScore 只传一个列表
# #   for (nm in names(score_list)) {
# #     feats <- score_list[[nm]]
# #     if (length(feats) < 3) {
# #       warning(paste0("Skipping ", nm, " (too few genes in dataset: ", length(feats), ")"))
# #       next
# #     }
# #     treg_main <- AddModuleScore(treg_main, features = list(feats), name = paste0(nm, "_Score"))
# #   }
# #
# #   # 保存
# #   saveRDS(treg_main, treg_main_scored_rds)
# # } else {
# #   cat(">> Found cached scored treg_main. Loading...\n")
# #   treg_main <- readRDS(treg_main_scored_rds)
# # }
# #
# # # 确认分数字段是否存在
# # score_cols <- grep("rTreg_Score|eTreg_Score|checkpoint_Score|cycling_Score", colnames(treg_main@meta.data), value = TRUE)
# # cat("Detected module score columns:\n")
# # print(score_cols)
# #
# # # ---------- 6. 输出“按cluster汇总的程序分数表”（用于命名） ----------
# # # 注意：AddModuleScore 会生成例如 rTreg_Score1 这样的列名
# # # 我们动态识别列名
# # get_score_col <- function(prefix) {
# #   # 匹配 prefix_Score 开头且以数字结尾的列名（例如 rTreg_Score1）
# #   hit <- grep(paste0("^", prefix, "_Score\\d+$"), colnames(treg_main@meta.data), value = TRUE)
# #   if (length(hit) == 0) return(NA_character_)
# #   hit[1]
# # }
# #
# # col_rTreg <- get_score_col("rTreg")
# # col_eTreg <- get_score_col("eTreg")
# # col_chk  <- get_score_col("checkpoint")
# # col_cyc  <- get_score_col("cycling")
# #
# # score_table <- data.frame(
# #   cluster = as.character(Idents(treg_main)),
# #   rTreg = if (!is.na(col_rTreg)) treg_main@meta.data[[col_rTreg]] else NA,
# #   eTreg = if (!is.na(col_eTreg)) treg_main@meta.data[[col_eTreg]] else NA,
# #   checkpoint = if (!is.na(col_chk)) treg_main@meta.data[[col_chk]] else NA,
# #   cycling = if (!is.na(col_cyc)) treg_main@meta.data[[col_cyc]] else NA
# # )
# #
# # # 汇总：每簇均值
# # score_summary <- aggregate(. ~ cluster, data = score_table, FUN = mean)
# # score_summary$n <- as.integer(table(score_table$cluster)[score_summary$cluster])
# # score_summary <- score_summary[order(score_summary$cluster), ]
# #
# # write.csv(score_summary, "./temp/Step1_TregMain_module_score_summary_by_cluster.csv", row.names = FALSE)
# #
# # # ---------- 7. 输出程序分数的可视化（对审稿人友好） ----------
# # # 7.1 DotPlot：程序分数（其实是meta数据，不是基因表达，所以用 VlnPlot 更直观）
# # # 这里给最稳的：VlnPlot 直接看分数分布
# # score_cols2 <- base_intersect(c(col_rTreg, col_eTreg, col_chk, col_cyc), colnames(treg_main@meta.data))
# # score_cols2 <- score_cols2[!is.na(score_cols2)]
# #
# # if (length(score_cols2) > 0) {
# #   p_score_vln <- VlnPlot(treg_main, features = score_cols2, group.by = "seurat_clusters", pt.size = 0) +
# #     ggtitle("Treg module scores by cluster")
# #   ggsave("./temp/Step1_TregMain_module_scores_VlnPlot.pdf", p_score_vln, width = 14, height = 7)
# # }
# #
# # # 7.2 FeaturePlot：把分数投到UMAP上（很直观）
# # # FeaturePlot 可以画meta列（Seurat支持）
# # if (!is.na(col_eTreg)) {
# #   p_eTreg <- FeaturePlot(treg_main, features = col_eTreg, reduction = "umap") + ggtitle("eTreg module score")
# #   ggsave("./temp/Step1_TregMain_eTregScore_FeaturePlot.pdf", p_eTreg, width = 10, height = 8)
# # }
# # if (!is.na(col_rTreg)) {
# #   p_rTreg <- FeaturePlot(treg_main, features = col_rTreg, reduction = "umap") + ggtitle("rTreg module score")
# #   ggsave("./temp/Step1_TregMain_rTregScore_FeaturePlot.pdf", p_rTreg, width = 10, height = 8)
# # }
# # if (!is.na(col_chk)) {
# #   p_chk <- FeaturePlot(treg_main, features = col_chk, reduction = "umap") + ggtitle("Checkpoint module score")
# #   ggsave("./temp/Step1_TregMain_CheckpointScore_FeaturePlot.pdf", p_chk, width = 10, height = 8)
# # }
# # if (!is.na(col_cyc)) {
# #   p_cyc <- FeaturePlot(treg_main, features = col_cyc, reduction = "umap") + ggtitle("Cycling module score")
# #   ggsave("./temp/Step1_TregMain_CyclingScore_FeaturePlot.pdf", p_cyc, width = 10, height = 8)
# # }
# #
# # cat("==== Step 1: Treg subclass naming DONE ====\n")
# #
# # #
# # # ======================================================================
# # # Step 1 后续：从“亚群命名完成”走到“组间描述性比较 + pseudo-bulk logFC 排名表”
# # # 重要边界（写 Methods/Results 时必须明确）：
# # #   - 每组只有 1 个 scRNA 测序库（虽然是 3 mice pooled）
# # #   - 因此组间比较仅做“描述性趋势”，不输出 p 值、不做显著性检验
# # # ======================================================================
# #
# # cat("==== Step 1.1+: Post-Step1 analyses START ====\n")
# # dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# #
# # # 强制使用 base 的集合函数，避免 conflicted 报错
# # base_intersect <- base::intersect
# # base_setdiff   <- base::setdiff
# # base_unname    <- base::unname
# #
# # DefaultAssay(treg_main) <- "RNA"
# # Idents(treg_main) <- "seurat_clusters"
# #
# # # ----------------------------------------------------------------------
# # # Step 1.1：把 cluster -> “可发表标签”（rTreg/eTreg/cycling/checkpoint-high）
# # # 目的：
# # #   - 明确图中每个簇的可发表命名（用于主文/图注）
# # # 产出：
# # #   - ./temp/Step1.1_TregMain_UMAP_subclass_labels.pdf
# # #   - ./temp/Step1.1_TregMain_cluster_to_label.csv（记录你这版命名口径，防止后续混乱）
# # # ----------------------------------------------------------------------
# #
# # cat("---- Step 1.1: Assign publishable subclass labels ----\n")
# #
# # # 你目前的命名（按你已完成的 Step1 结论）：
# # #   cluster 0 = eTreg
# # #   cluster 3 = cycling
# # #   cluster 4 = rTreg / transitional
# # # checkpoint-high 通常不是独立 cluster，而是 eTreg 或某一簇内的“程序状态”
# # # 因此这里给两层命名：
# # #   (A) cluster-level 主标签：rTreg / eTreg / cycling
# # #   (B) checkpoint-high：用 checkpoint module score 作为“状态标签”附加（下方 Step1.2）
# # cluster_to_label <- c(
# #   "0" = "eTreg-like",
# #   "3" = "cycling Treg",
# #   "4" = "rTreg-ish / transitional"
# # )
# #
# # treg_main$Treg_subclass_cluster <- base_unname(cluster_to_label[as.character(Idents(treg_main))])
# # treg_main$Treg_subclass_cluster[is.na(treg_main$Treg_subclass_cluster)] <- "Unassigned"
# #
# # # 保存口径（非常建议：以后写论文就引用这张表，避免改来改去）
# # df_map <- data.frame(
# #   cluster = names(cluster_to_label),
# #   label = as.character(cluster_to_label),
# #   stringsAsFactors = FALSE
# # )
# # write.csv(df_map, "./temp/Step1.1_TregMain_cluster_to_label.csv", row.names = FALSE)
# #
# # # 输出：按“可发表标签”上色的 UMAP
# # p_umap_label <- DimPlot(treg_main, reduction = "umap", group.by = "Treg_subclass_cluster", label = TRUE) +
# #   ggtitle("Treg subclasses (cluster-level publishable labels)")
# # ggsave("./temp/Step1.1_TregMain_UMAP_subclass_labels.pdf", p_umap_label, width = 12, height = 8)
# #
# # # ----------------------------------------------------------------------
# # # Step 1.2：定义 checkpoint-high “状态标签”（基于 module score 的阈值）
# # # 目的：
# # #   - checkpoint-high 往往不是独立 cluster，而是某些细胞的状态（Pdcd1/Lag3/Tigit 高）
# # #   - 用 module score 给出一个“审稿人友好”的状态标注，供 Results 描述或补图展示
# # # 产出：
# # #   - ./temp/Step1.2_TregMain_UMAP_checkpoint_state.pdf
# # # ----------------------------------------------------------------------
# #
# # cat("---- Step 1.2: Define checkpoint-high state (score-based) ----\n")
# #
# # # 复用你 Step1 里写的 get_score_col 逻辑（如果当前环境没有，就定义一次）
# # get_score_col2 <- function(obj, prefix) {
# #   hit <- grep(paste0("^", prefix, "_Score\\d+$"), colnames(obj@meta.data), value = TRUE)
# #   if (length(hit) == 0) return(NA_character_)
# #   hit[1]
# # }
# # col_chk <- get_score_col2(treg_main, "checkpoint")
# #
# # # 经验阈值：用分数的上分位数定义 checkpoint-high（不做统计，只做状态标注）
# # # 你可把 0.80 调成 0.75/0.85（更宽松/更严格）
# # if (!is.na(col_chk)) {
# #   thr_q <- 0.80
# #   thr <- as.numeric(stats::quantile(treg_main@meta.data[[col_chk]], probs = thr_q, na.rm = TRUE))
# #   treg_main$Checkpoint_state <- ifelse(treg_main@meta.data[[col_chk]] >= thr, "checkpoint-high", "checkpoint-low")
# #
# #   p_chk_state <- DimPlot(treg_main, reduction = "umap", group.by = "Checkpoint_state") +
# #     ggtitle(paste0("Checkpoint state (top ", thr_q*100, "% by checkpoint score; descriptive)"))
# #   ggsave("./temp/Step1.2_TregMain_UMAP_checkpoint_state.pdf", p_chk_state, width = 12, height = 8)
# # } else {
# #   warning("No checkpoint module score column found. Skip checkpoint-high state labeling.")
# # }
# #
# # # 保存带标签的对象（后面所有组间分析建议从这里读取）
# # saveRDS(treg_main, "./temp/treg_main_labeled_step1.2.rds")
# #
# # # ----------------------------------------------------------------------
# # # Step 1.5：识别处理组字段（group）
# # # 目的：
# # #   - 你的组间比较需要一个“处理组列”
# # # 说明：
# # #   - 若你已有 treg_main$group，则直接使用
# # #   - 否则自动尝试：condition / treatment / group / orig.ident
# # #   - 你也可以手动指定：group_col <- "你的列名"
# # # ----------------------------------------------------------------------
# #
# # cat("---- Step 1.5: Detect group column ----\n")
# #
# # group_col <- NA_character_
# # candidate_cols <- c("group", "condition", "treatment", "Group", "Condition", "Treatment", "orig.ident")
# # for (cc in candidate_cols) {
# #   if (cc %in% colnames(treg_main@meta.data)) { group_col <- cc; break }
# # }
# # if (is.na(group_col)) {
# #   stop("treg_main@meta.data 里找不到处理组字段。请把处理组写入 meta，例如：treg_main$group <- ...")
# # } else {
# #   cat("Using group column:", group_col, "\n")
# #   treg_main$group_use <- as.character(treg_main@meta.data[[group_col]])
# # }
# #
# # # ----------------------------------------------------------------------
# # # Step 1.6：组间比较（描述性）：composition（亚群比例）+（可选）下采样敏感性分析
# # # 目的（适合写 Results 的“现象段”）：
# # #   - 描述不同处理组在 Treg 亚群组成上的趋势差异（不做p值）
# # #   - 用“均衡细胞数下采样”检查趋势是否仅由细胞数差异导致
# # # 产出：
# # #   - ./temp/Step1.6_group_by_subclass_counts.csv
# # #   - ./temp/Step1.6_group_by_subclass_proportions.csv
# # #   - ./temp/Step1.6_group_composition_barplot.pdf
# # #   - （可选）./temp/Step1.6_group_composition_downsampled_barplot.pdf
# # # ----------------------------------------------------------------------
# #
# # cat("---- Step 1.6: Group-level composition (descriptive) ----\n")
# #
# # df_comp <- data.frame(
# #   cell = colnames(treg_main),
# #   group = treg_main$group_use,
# #   subclass = treg_main$Treg_subclass_cluster,
# #   stringsAsFactors = FALSE
# # )
# #
# # tab_n <- table(df_comp$group, df_comp$subclass)
# # tab_prop <- prop.table(tab_n, margin = 1)
# #
# # write.csv(as.data.frame.matrix(tab_n), "./temp/Step1.6_group_by_subclass_counts.csv")
# # write.csv(as.data.frame.matrix(tab_prop), "./temp/Step1.6_group_by_subclass_proportions.csv")
# #
# # # 画：堆叠比例柱状图
# # if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# # library(ggplot2)
# #
# # df_prop <- as.data.frame(tab_prop)
# # colnames(df_prop) <- c("group","subclass","prop")
# #
# # p_comp <- ggplot(df_prop, aes(x = group, y = prop, fill = subclass)) +
# #   geom_col(position = "fill") +
# #   ylab("Proportion within group (descriptive)") +
# #   xlab("Group (1 pooled library/group)") +
# #   ggtitle("Treg subclass composition by group (descriptive; no p-values)") +
# #   theme_classic()
# #
# # ggsave("./temp/Step1.6_group_composition_barplot.pdf", p_comp, width = 10, height = 6)
# #
# # # （可选）下采样敏感性分析：均衡每组细胞数再算比例，避免“细胞数不均衡”伪差异
# # do_downsample_sensitivity <- TRUE
# # if (do_downsample_sensitivity) {
# #   set.seed(123)
# #   n_by_group <- table(df_comp$group)
# #   min_n <- min(n_by_group)
# #   B <- 50
# #
# #   groups <- names(n_by_group)
# #   subclasses <- sort(unique(df_comp$subclass))
# #
# #   prop_mat <- matrix(0, nrow = B, ncol = length(groups)*length(subclasses))
# #   colnames(prop_mat) <- as.vector(outer(groups, subclasses, paste, sep="__"))
# #
# #   for (b in seq_len(B)) {
# #     df_ds <- do.call(rbind, lapply(groups, function(g) {
# #       df_g <- df_comp[df_comp$group == g, , drop = FALSE]
# #       df_g[sample(seq_len(nrow(df_g)), min_n), , drop = FALSE]
# #     }))
# #     tab_ds <- table(df_ds$group, df_ds$subclass)
# #     prop_ds <- prop.table(tab_ds, margin = 1)
# #
# #     for (g in groups) for (s in subclasses) {
# #       key <- paste(g, s, sep="__")
# #       prop_mat[b, key] <- ifelse(s %in% colnames(prop_ds), prop_ds[g, s], 0)
# #     }
# #   }
# #
# #   prop_mean <- colMeans(prop_mat)
# #   df_ds_mean <- do.call(rbind, lapply(groups, function(g) {
# #     data.frame(group=g, subclass=subclasses,
# #                prop=prop_mean[paste(g, subclasses, sep="__")],
# #                stringsAsFactors = FALSE)
# #   }))
# #   write.csv(df_ds_mean, "./temp/Step1.6_group_composition_downsampled_proportions.csv", row.names = FALSE)
# #
# #   p_comp_ds <- ggplot(df_ds_mean, aes(x = group, y = prop, fill = subclass)) +
# #     geom_col(position = "fill") +
# #     ylab(paste0("Proportion within group (downsampled to n=", min_n, ", mean over B=", B, ")")) +
# #     xlab("Group") +
# #     ggtitle("Composition sensitivity check (descriptive)") +
# #     theme_classic()
# #
# #   ggsave("./temp/Step1.6_group_composition_downsampled_barplot.pdf", p_comp_ds, width = 10, height = 6)
# # }
# #
# # saveRDS(list(tab_n=tab_n, tab_prop=tab_prop), "./temp/Step1.6_group_composition_cache.rds")
# #
# # # ----------------------------------------------------------------------
# # # Step 1.7 (FINAL, MERGED): cycling Treg publishable pseudo-bulk logFC (no p-values)
# # # ----------------------------------------------------------------------
# # # ============================================================
# # # Step 1.7 (PUBLISHABLE, CYCLING-ONLY):
# # #   - ensure subclass labels exist (patch from CSV if needed)
# # #   - Treg core gate
# # #   - extract cycling Treg
# # #   - [PATCH v2] publishable cleaning FIX: avoid B/NK any-hit false-kill
# # #   - pseudo-bulk logFC ranking (no p-values) with gene filters
# # #   - single-cell validation plots (FeaturePlot/VlnPlot)
# # #   - export per-group cell IDs (for Methods/Supplement)
# # # # ============================================================主分析cyclingTreg细胞开始
# # # cat("==== Step 1.7 开始【主分析cyclingTreg细胞】 \n")
# # # # ============================================================
# # # cat("==== Step 1.7 (CYCLING-ONLY PUBLISHABLE): START ====\n")
# # # suppressPackageStartupMessages({
# # #   library(Seurat)
# # #   library(Matrix)
# # #   library(ggplot2)
# # # })
# # #
# # # dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
# # #
# # # # ---------- 0) 防 conflicted：显式用 base ----------
# # # base_intersect <- base::intersect
# # # base_setdiff   <- base::setdiff
# # # base_unname    <- base::unname
# # #
# # # # ---------- 1) 读取对象（优先你已修复并带 subclass labels 的版本） ----------
# # # treg_rds1 <- "./temp/treg_main_scored_step1_with_subclass_labels.rds"
# # # treg_rds2 <- "./temp/treg_main_scored_step1.rds"
# # #
# # # if (file.exists(treg_rds1)) {
# # #   treg_main <- readRDS(treg_rds1)
# # #   cat(">> Loaded:", treg_rds1, "\n")
# # # } else if (file.exists(treg_rds2)) {
# # #   treg_main <- readRDS(treg_rds2)
# # #   cat(">> Loaded:", treg_rds2, "\n")
# # # } else {
# # #   stop("找不到 treg_main RDS：\n  ", treg_rds1, "\n  ", treg_rds2,
# # #        "\n请先完成 Step1 并保存 treg_main。")
# # # }
# # #
# # # DefaultAssay(treg_main) <- "RNA"
# # #
# # # # ---------- 2) 关键列名（按你现在的实际列名） ----------
# # # group_col    <- "group"
# # # subclass_col <- "Treg_subclass_cluster"  # 你 Step1.2/patch 用的这列
# # #
# # # if (!(group_col %in% colnames(treg_main@meta.data))) {
# # #   stop("meta.data 里找不到 group 列：", group_col,
# # #        "\n请改成你的真实列名（例如 orig.ident）。")
# # # }
# # #
# # # # ---------- 3) 如 subclass_col 不存在：用 cluster->label CSV 现场补回（Patch 合并在这里） ----------
# # # if (!(subclass_col %in% colnames(treg_main@meta.data))) {
# # #   cat(">> subclass_col not found. Running patch from CSV...\n")
# # #
# # #   map_file <- "./temp/Step1.1_TregMain_cluster_to_label.csv"
# # #   if (!file.exists(map_file)) {
# # #     stop("找不到 cluster->label 映射表：", map_file,
# # #          "\n请确认该 CSV 已生成并在 ./temp 下。")
# # #   }
# # #
# # #   mp <- read.csv(map_file, stringsAsFactors = FALSE)
# # #   if (!all(c("cluster","label") %in% colnames(mp))) {
# # #     stop("映射表必须包含列 cluster,label。你当前列名为：",
# # #          paste(colnames(mp), collapse = ", "))
# # #   }
# # #
# # #   Idents(treg_main) <- "seurat_clusters"
# # #   cl <- as.character(Idents(treg_main))
# # #
# # #   map_vec <- setNames(mp$label, as.character(mp$cluster))
# # #   treg_main@meta.data[[subclass_col]] <- base_unname(map_vec[cl])
# # #
# # #   cat("Subclass label distribution (after patch):\n")
# # #   print(sort(table(treg_main@meta.data[[subclass_col]], useNA = "ifany"), decreasing = TRUE))
# # #
# # #   saveRDS(treg_main, treg_rds1)
# # #   cat(">> Saved patched object to:", treg_rds1, "\n")
# # # }
# # #
# # # # ---------- 4) 入口：Treg core gate ----------
# # # need_core <- c("Foxp3","Il2ra","Trac","Cd3e")
# # # need_core <- base_intersect(need_core, rownames(treg_main))
# # # if (length(need_core) < 3) {
# # #   stop("关键 Treg core 基因缺失过多（Foxp3/Il2ra/Trac/Cd3e）。请检查 gene symbol 命名体系。")
# # # }
# # #
# # # treg_core <- subset(treg_main, subset = Foxp3 > 0 & Il2ra > 0 & Trac > 0 & Cd3e > 0)
# # # cat("Treg core cells:", ncol(treg_core), "\n")
# # # cat("Treg core per group:\n")
# # # print(sort(table(as.character(treg_core@meta.data[[group_col]])), decreasing = TRUE))
# # #
# # # # ---------- 5) 取 cycling Treg（最稳：用 cells 向量） ----------
# # # meta_core <- treg_core@meta.data
# # # cells_cycling <- rownames(meta_core)[as.character(meta_core[[subclass_col]]) == "cycling Treg"]
# # #
# # # cat("cycling cells (core, before clean):", length(cells_cycling), "\n")
# # # if (length(cells_cycling) == 0) {
# # #   cat("Available subclass labels:\n")
# # #   print(sort(table(as.character(meta_core[[subclass_col]])), decreasing = TRUE))
# # #   stop("未找到 cycling Treg。请检查 subclass label 字符串是否完全一致（空格/大小写）。")
# # # }
# # #
# # # cycling <- subset(treg_core, cells = cells_cycling)
# # # cat("cycling object cells:", ncol(cycling), "\n")
# # # cat("cycling cells per group:\n")
# # # print(sort(table(as.character(cycling@meta.data[[group_col]])), decreasing = TRUE))
# # #
# # # # ---------- 6) cycling publishable clean（[PATCH v2] 修复 B/NK any-hit 误杀） ----------
# # # # 说明：
# # # #   你之前用 B any / NK any 会导致 cycling Treg 被整批误杀（尤其 s40kD）
# # # #   这里改为：
# # # #     - B STRONG：>=2 且命中至少一个 B-core
# # # #     - NK STRONG：>=2 且命中至少一个 NK-core（并删除 Trdc/Trgc1）
# # # #     - melanoma / non-immune：改为 STRONG(>=2)，避免 any-hit 误杀
# # # thr_core <- 0.25
# # # thr_bad  <- 0.25
# # #
# # # DefaultAssay(cycling) <- "RNA"
# # #
# # # # --- core sets ---
# # # genes_Tcore    <- c("Trac","Cd3e","Cd3d")
# # # genes_Tregcore <- c("Foxp3","Il2ra","Ikzf2","Ctla4")
# # #
# # # # --- myeloid strong (>=2) ---
# # # genes_myeloid  <- c("Tyrobp","Lyz2","Csf1r","LST1","Itgam","S100a8","S100a9","Aif1","Fcgr3","Lgals3")
# # #
# # # # [PATCH v2] B: core + aux（避免 Igkc/Cd74 环境 RNA 误杀）
# # # genes_B_core <- c("Ms4a1","Cd79a","Cd79b","Jchain")
# # # genes_B_aux  <- c("Cd74","Ighm","Igkc")  # aux 参与计数，但必须 corehit 才算 strong
# # # genes_B_all  <- unique(c(genes_B_core, genes_B_aux))
# # #
# # # # [PATCH v2] NK/CTL: 删除 Trdc/Trgc1（γδT），用 NK STRONG（>=2 + corehit）
# # # genes_NK_core <- c("Ncr1","Klrb1c","Klrk1","Nkg7")
# # # genes_NK_aux  <- c("Prf1","Gzmb","Gzma","Cst7")
# # # genes_NK_all  <- unique(c(genes_NK_core, genes_NK_aux))
# # #
# # # # RBC any（保持）
# # # genes_RBC <- c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2")
# # #
# # # # [PATCH v2] non-immune：STRONG(>=2)
# # # genes_epi   <- c("Epcam","Tacstd2","Krt8","Krt18","Krt19","Krt14","Krt17")
# # # genes_endo  <- c("Pecam1","Kdr","Vwf","Eng","Ramp2")
# # # genes_fibro <- c("Col1a1","Col1a2","Dcn","Lum","Col3a1","Pdgfra")
# # # genes_peri  <- c("Acta2","Rgs5","Pdgfrb","Des","Cspg4","Tagln")
# # # genes_nonimmune_all <- unique(c(genes_epi, genes_endo, genes_fibro, genes_peri))
# # #
# # # # [PATCH v2] melanoma：STRONG(>=2)
# # # genes_melanoma <- c("Sox10","Pmel","Dct","Mlana","Tyrp1","Slc45a2","Trpm1","Slc24a5","Mitf","Tyr")
# # #
# # # panel_all <- unique(c(
# # #   genes_Tcore, genes_Tregcore,
# # #   genes_myeloid, genes_B_all, genes_NK_all, genes_RBC,
# # #   genes_nonimmune_all, genes_melanoma
# # # ))
# # # panel_present <- base_intersect(panel_all, rownames(cycling))
# # # cat("Panel requested:", length(panel_all), " | present:", length(panel_present), "\n")
# # #
# # # if (length(panel_present) < 10) {
# # #   stop("panel_present 太少（<10）。说明你对象里 gene 命名可能不是 symbol（例如 Ensembl）。请先统一命名体系。")
# # # }
# # #
# # # df <- FetchData(cycling, vars = panel_present)
# # # grp_raw <- as.character(cycling@meta.data[[group_col]])
# # # names(grp_raw) <- rownames(cycling@meta.data)
# # #
# # # # helper：any-hit 与 hit-count
# # # hit_any <- function(gset, thr = thr_bad) {
# # #   g2 <- base_intersect(gset, colnames(df))
# # #   if (length(g2) == 0) return(rep(FALSE, nrow(df)))
# # #   out <- rep(FALSE, nrow(df))
# # #   for (g in g2) out <- out | (df[[g]] > thr)
# # #   out
# # # }
# # # hit_n <- function(gset, thr = thr_bad) {
# # #   g2 <- base_intersect(gset, colnames(df))
# # #   if (length(g2) == 0) return(rep(0L, nrow(df)))
# # #   k <- rep(0L, nrow(df))
# # #   for (g in g2) k <- k + as.integer(df[[g]] > thr)
# # #   k
# # # }
# # #
# # # # core flags
# # # flag_T <- hit_any(genes_Tcore, thr = thr_core)
# # # flag_Foxp3 <- if ("Foxp3" %in% colnames(df)) (df$Foxp3 > thr_core) else rep(FALSE, nrow(df))
# # # flag_Il2ra <- if ("Il2ra" %in% colnames(df)) (df$Il2ra > thr_core) else rep(FALSE, nrow(df))
# # # flag_Treg_strict <- flag_Foxp3 & flag_Il2ra
# # #
# # # # myeloid STRONG >=2
# # # flag_my_strong <- hit_n(genes_myeloid, thr = thr_bad) >= 2
# # #
# # # # [PATCH v2] B STRONG: >=2 AND corehit
# # # B_hits <- hit_n(genes_B_all, thr = thr_bad)
# # # flag_B_corehit <- hit_any(genes_B_core, thr = thr_bad)
# # # flag_B_strong <- (B_hits >= 2) & flag_B_corehit
# # #
# # # # [PATCH v2] NK STRONG: >=2 AND corehit
# # # NK_hits <- hit_n(genes_NK_all, thr = thr_bad)
# # # flag_NK_corehit <- hit_any(genes_NK_core, thr = thr_bad)
# # # flag_NK_strong <- (NK_hits >= 2) & flag_NK_corehit
# # #
# # # # RBC any
# # # flag_RBC <- hit_any(genes_RBC, thr = thr_bad)
# # #
# # # # [PATCH v2] non-immune STRONG >=2
# # # flag_nonimmune_strong <- hit_n(genes_nonimmune_all, thr = thr_bad) >= 2
# # #
# # # # [PATCH v2] melanoma STRONG >=2
# # # flag_mel_strong <- hit_n(genes_melanoma, thr = thr_bad) >= 2
# # #
# # # # 诊断打印（非常关键：看是否仍对 s40kD “一刀切”）
# # # rate_by_group <- function(flag, nm) {
# # #   cells <- rownames(df)
# # #   g <- grp_raw[cells]
# # #   cat("\nRule hit-rate:", nm, "\n")
# # #   print(tapply(flag, g, mean))
# # # }
# # #
# # # cat("\nQC hit-rates by group (cycling, before clean):\n")
# # # cat("cycling per group:\n"); print(table(grp_raw[rownames(df)]))
# # # rate_by_group(flag_Treg_strict, "Treg strict (Foxp3&Il2ra)")
# # # rate_by_group(flag_my_strong, "Myeloid STRONG (>=2)")
# # # rate_by_group(flag_B_strong,  "B STRONG (>=2 + corehit)")
# # # rate_by_group(flag_NK_strong, "NK STRONG (>=2 + corehit)")
# # # rate_by_group(flag_RBC, "RBC any")
# # # rate_by_group(flag_nonimmune_strong, "Non-immune STRONG (>=2)")
# # # rate_by_group(flag_mel_strong, "Melanoma STRONG (>=2)")
# # #
# # # # [PATCH v2] keep：只剔除 strong 污染（不再用 B/NK any）
# # # keep <- flag_T & flag_Treg_strict &
# # #   (!flag_my_strong) &
# # #   (!flag_B_strong) &
# # #   (!flag_NK_strong) &
# # #   (!flag_RBC) &
# # #   (!flag_nonimmune_strong) &
# # #   (!flag_mel_strong)
# # #
# # # cat("\nPublishable keep fraction by group:\n")
# # # print(tapply(keep, grp_raw[rownames(df)], mean))
# # # cat("Publishable kept cell counts by group:\n")
# # # print(tapply(keep, grp_raw[rownames(df)], sum))
# # #
# # # keep_cells <- rownames(df)[keep]
# # # cycling_clean <- subset(cycling, cells = keep_cells)
# # #
# # # cat("\ncycling_clean cells:", ncol(cycling_clean), "\n")
# # # cat("cycling_clean cells per group:\n")
# # # print(sort(table(as.character(cycling_clean@meta.data[[group_col]])), decreasing = TRUE))
# # #
# # # cycling_clean_rds <- "./temp/cycling_treg_clean_for_pseudobulk_publishable.rds"
# # # saveRDS(cycling_clean, cycling_clean_rds)
# # # cat(">> Saved cycling_clean RDS:", cycling_clean_rds, "\n")
# # #
# # # # ---------- 7) 导出每组 cell IDs ----------
# # # outdir <- "./temp/Step1.7_cycling_publishable"
# # # dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# # #
# # # grp <- as.character(cycling_clean@meta.data[[group_col]])
# # # names(grp) <- colnames(cycling_clean)
# # #
# # # write.csv(
# # #   data.frame(cell = names(grp), group = grp, stringsAsFactors = FALSE),
# # #   file = file.path(outdir, "cycling_clean_cellIDs_all.csv"),
# # #   row.names = FALSE
# # # )
# # #
# # # for (g in sort(unique(grp))) {
# # #   cells_g <- names(grp)[grp == g]
# # #   write.csv(
# # #     data.frame(cell = cells_g, stringsAsFactors = FALSE),
# # #     file = file.path(outdir, paste0("cycling_clean_cellIDs_", g, ".csv")),
# # #     row.names = FALSE
# # #   )
# # # }
# # #
# # # # ---------- 8) 单细胞层面回查关键基因 ----------
# # # val_genes <- c("Il2ra","Ctla4","Il2rb","Foxp3","Mki67","Top2a")
# # # val_genes <- base_intersect(val_genes, rownames(cycling_clean))
# # #
# # # if (length(val_genes) == 0) {
# # #   warning("验证基因在对象里一个都找不到：请确认基因命名体系。跳过 Feature/Vln。")
# # # } else {
# # #   p_fp <- FeaturePlot(cycling_clean, features = val_genes, reduction = "umap", order = TRUE)
# # #   ggsave(file.path(outdir, "cycling_clean_validation_FeaturePlot.pdf"), p_fp, width = 12, height = 8)
# # #
# # #   p_vln <- VlnPlot(cycling_clean, features = val_genes, group.by = group_col, pt.size = 0)
# # #   ggsave(file.path(outdir, "cycling_clean_validation_VlnPlot_by_group.pdf"), p_vln, width = 12, height = 6)
# # # }
# # #
# # # # ---------- 9) pseudo-bulk logFC（不输出 p 值；描述性排名） ----------
# # # counts <- LayerData(cycling_clean, assay = "RNA", layer = "counts")
# # # if (is.null(counts)) stop("找不到 RNA counts layer（LayerData返回NULL）。")
# # #
# # # tab <- table(grp)
# # # cat("\nCells per group (cycling_clean):\n"); print(tab)
# # #
# # # # 至少两个组即可跑（如果 s40kD 仍太少，至少保证 control1 vs s1_hIL）
# # # min_cells_candidates <- c(20, 10, 5, 2, 1)
# # # ok_groups <- character(0)
# # # used_min_cells <- NA_integer_
# # #
# # # for (m in min_cells_candidates) {
# # #   cand <- names(tab)[as.integer(tab) >= m]
# # #   if (length(cand) >= 2) {
# # #     ok_groups <- cand
# # #     used_min_cells <- m
# # #     break
# # #   }
# # # }
# # #
# # # if (length(ok_groups) < 2) {
# # #   stop("满足最小细胞数阈值的组别不足2个：\n",
# # #        paste(names(tab), as.integer(tab), sep="=", collapse=", "),
# # #        "\n请降低过滤强度或检查是否把某组清空了。")
# # # }
# # #
# # # cat("Groups kept (min cells >=", used_min_cells, "): ", paste(ok_groups, collapse = ", "), "\n", sep="")
# # #
# # # pb_counts <- sapply(ok_groups, function(g) {
# # #   cells <- names(grp)[grp == g]
# # #   Matrix::rowSums(counts[, cells, drop = FALSE])
# # # })
# # # rownames(pb_counts) <- rownames(counts)
# # #
# # # libsize <- colSums(pb_counts)
# # # cat("Pseudo-bulk library sizes:\n"); print(libsize)
# # #
# # # cpm <- sweep(pb_counts, 2, libsize, FUN = "/") * 1e6
# # # logcpm <- log2(cpm + 1)
# # #
# # # write.csv(pb_counts, file.path(outdir, "pseudobulk_counts_gene_x_group.csv"))
# # # write.csv(logcpm,   file.path(outdir, "pseudobulk_logCPM_gene_x_group.csv"))
# # #
# # # # gene 过滤：避免极低 counts 冲榜
# # # min_max_counts <- 10
# # # min_sum_counts <- 20
# # # min_aveexpr    <- 1.0
# # #
# # # rank_one_pair <- function(g1, g2) {
# # #   c1 <- pb_counts[, g1]; c2 <- pb_counts[, g2]
# # #   AveExpr <- rowMeans(logcpm[, c(g1, g2), drop = FALSE])
# # #
# # #   keep_gene <- (pmax(c1, c2) >= min_max_counts) & ((c1 + c2) >= min_sum_counts) & (AveExpr >= min_aveexpr)
# # #
# # #   df_rank <- data.frame(
# # #     gene = rownames(logcpm),
# # #     logFC = as.numeric(logcpm[, g1] - logcpm[, g2]),
# # #     AveExpr = as.numeric(AveExpr),
# # #     c1 = as.integer(c1),
# # #     c2 = as.integer(c2),
# # #     stringsAsFactors = FALSE
# # #   )
# # #   df_rank <- df_rank[keep_gene, , drop = FALSE]
# # #   df_rank <- df_rank[order(df_rank$logFC, decreasing = TRUE), , drop = FALSE]
# # #
# # #   tag <- paste0(g1, "_vs_", g2)
# # #   write.csv(df_rank, file.path(outdir, paste0("logFC_rank_", tag, ".csv")), row.names = FALSE)
# # #   write.csv(head(df_rank, 200), file.path(outdir, paste0("logFC_rank_", tag, "_TOP200.csv")), row.names = FALSE)
# # #   write.csv(tail(df_rank, 200), file.path(outdir, paste0("logFC_rank_", tag, "_BOTTOM200.csv")), row.names = FALSE)
# # #
# # #   if (nrow(df_rank) > 0) {
# # #     topn <- min(20, nrow(df_rank))
# # #     df_up <- head(df_rank, topn)
# # #     df_dn <- tail(df_rank, topn)
# # #     df_dn <- df_dn[order(df_dn$logFC, decreasing = FALSE), ]
# # #
# # #     df_up$dir <- "Up"; df_dn$dir <- "Down"
# # #     df_plot <- rbind(df_up, df_dn)
# # #     df_plot$gene <- factor(df_plot$gene, levels = df_plot$gene)
# # #
# # #     p_bar <- ggplot(df_plot, aes(x = gene, y = logFC)) +
# # #       geom_col() + coord_flip() + facet_wrap(~dir, scales = "free_y") +
# # #       theme_bw() + ggtitle(paste0("cycling_clean | ", tag, " | Top logFC genes (no p-values)"))
# # #     ggsave(file.path(outdir, paste0("TopGenes_Barplot_", tag, ".pdf")), p_bar, width = 10, height = 8)
# # #
# # #     p_ma <- ggplot(df_rank, aes(x = AveExpr, y = logFC)) +
# # #       geom_point(size = 0.5, alpha = 0.5) +
# # #       geom_hline(yintercept = 0, linetype = 2) +
# # #       theme_bw() + ggtitle(paste0("cycling_clean | ", tag, " | MA plot (logCPM)"))
# # #     ggsave(file.path(outdir, paste0("MAplot_", tag, ".pdf")), p_ma, width = 8, height = 6)
# # #   }
# # #
# # #   cat("Done:", tag, " genes kept =", nrow(df_rank), "\n")
# # # }
# # #
# # # pairs <- combn(ok_groups, 2, simplify = FALSE)
# # # for (pp in pairs) rank_one_pair(pp[1], pp[2])
# # #
# # # cat("\nOutputs saved under: ", outdir, "\n", sep="")
# # # cat("==== Step 1.7 (CYCLING-ONLY PUBLISHABLE): DONE ====\n")
# # #
# # # # ===================== Step 1.7 (POSTCHECK): cycling_clean single-cell validation + cellID export =====================
# # # cat("==== Step 1.7 POSTCHECK: cycling_clean validation + cellID export START ====\n")
# # #
# # # dir.create("./temp/Step1.7_POSTCHECK", showWarnings = FALSE, recursive = TRUE)
# # # out_post <- "./temp/Step1.7_POSTCHECK"
# # #
# # # # ---- 0) 读取 cycling_clean（优先内存，否则读缓存） ----
# # # cycling_clean_rds <- "./temp/cycling_treg_clean_for_pseudobulk_publishable.rds"
# # #
# # # if (exists("cycling_clean")) {
# # #   cat(">> cycling_clean found in memory. cells =", ncol(cycling_clean), "\n")
# # # } else if (file.exists(cycling_clean_rds)) {
# # #   cycling_clean <- readRDS(cycling_clean_rds)
# # #   cat(">> cycling_clean loaded from:", cycling_clean_rds, " cells =", ncol(cycling_clean), "\n")
# # # } else {
# # #   stop("找不到 cycling_clean，也找不到缓存：", cycling_clean_rds,
# # #        "\n请先运行 Step1.7 FINAL (MERGED) 生成 cycling_clean。")
# # # }
# # #
# # # DefaultAssay(cycling_clean) <- "RNA"
# # #
# # # # ---- 1) 基础列名检查（只改这里即可） ----
# # # group_col <- "group"
# # # if (!(group_col %in% colnames(cycling_clean@meta.data))) {
# # #   stop("meta.data 里找不到 group_col = ", group_col,
# # #        "\n请用 colnames(cycling_clean@meta.data) 查看真实列名，然后改 group_col。")
# # # }
# # #
# # # cat("Cells per group in cycling_clean:\n")
# # # print(table(cycling_clean@meta.data[[group_col]]))
# # #
# # # # ---- 2) 关键基因面板：在对象里实际存在的基因 ----
# # # base_intersect <- base::intersect
# # #
# # # genes_check <- c("Il2ra","Ctla4","Il2rb","Foxp3","Mki67","Top2a")
# # # genes_present <- base_intersect(genes_check, rownames(cycling_clean))
# # #
# # # cat("Requested genes:", paste(genes_check, collapse = ", "), "\n")
# # # cat("Present genes:", paste(genes_present, collapse = ", "), "\n")
# # #
# # # if (length(genes_present) == 0) {
# # #   cat("Example rownames(cycling_clean):\n"); print(head(rownames(cycling_clean), 20))
# # #   stop("上述基因在对象中全部缺失：很可能是基因命名体系不一致（Ensembl vs Symbol）。")
# # # }
# # #
# # # # ---- 3) UMAP 保障：如果 cycling_clean 没有 umap，则从 RNA 重新计算（只针对该子集，开销很小） ----
# # # # 目的：FeaturePlot 需要一个二维嵌入（umap/tsne），不想依赖上游对象一定携带
# # # if (!("umap" %in% names(cycling_clean@reductions))) {
# # #   cat(">> cycling_clean has no UMAP reduction. Recomputing UMAP for visualization...\n")
# # #   suppressPackageStartupMessages(library(Seurat))
# # #
# # #   # 这里用最稳的“轻量可视化管线”：Normalize -> HVG -> Scale -> PCA -> UMAP
# # #   cycling_clean <- NormalizeData(cycling_clean, verbose = FALSE)
# # #   cycling_clean <- FindVariableFeatures(cycling_clean, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# # #   cycling_clean <- ScaleData(cycling_clean, features = VariableFeatures(cycling_clean), verbose = FALSE)
# # #   cycling_clean <- RunPCA(cycling_clean, npcs = 30, verbose = FALSE)
# # #   cycling_clean <- RunUMAP(cycling_clean, dims = 1:20, verbose = FALSE)
# # #
# # #   saveRDS(cycling_clean, file.path(out_post, "cycling_clean_with_umap.rds"))
# # #   cat("Saved cycling_clean_with_umap.rds for fast reuse.\n")
# # # } else {
# # #   cat(">> cycling_clean already has UMAP. Using existing embedding.\n")
# # # }
# # #
# # # # ---- 4) FeaturePlot：关键基因回到单细胞层面（UMAP） ----
# # # suppressPackageStartupMessages(library(ggplot2))
# # # suppressPackageStartupMessages(library(Seurat))
# # #
# # # p_fp <- FeaturePlot(
# # #   cycling_clean,
# # #   features = genes_present,
# # #   reduction = "umap",
# # #   order = TRUE
# # # ) + ggtitle("cycling_clean: single-cell FeaturePlot (key genes)")
# # #
# # # ggsave(
# # #   filename = file.path(out_post, "cycling_clean_keygenes_FeaturePlot.pdf"),
# # #   plot = p_fp, width = 14, height = 10
# # # )
# # #
# # # # （可选）同时输出分组的 UMAP（检查是否某组缺失）
# # # p_umap_group <- DimPlot(
# # #   cycling_clean, reduction = "umap",
# # #   group.by = group_col, label = FALSE
# # # ) + ggtitle("cycling_clean: UMAP colored by group")
# # #
# # # ggsave(
# # #   filename = file.path(out_post, "cycling_clean_UMAP_by_group.pdf"),
# # #   plot = p_umap_group, width = 10, height = 8
# # # )
# # #
# # # # ---- 5) VlnPlot：按 group 的表达分布（重点：样本量很小时，看到是不是“极少数细胞拉高均值”） ----
# # # # 注意：pt.size 不设为0，便于你看到“每组到底几个点”
# # # p_vln <- VlnPlot(
# # #   cycling_clean,
# # #   features = genes_present,
# # #   group.by = group_col,
# # #   pt.size = 0.8
# # # ) + ggtitle("cycling_clean: VlnPlot by group (key genes)")
# # #
# # # ggsave(
# # #   filename = file.path(out_post, "cycling_clean_keygenes_VlnPlot_by_group.pdf"),
# # #   plot = p_vln, width = 14, height = 8
# # # )
# # #
# # # # ---- 6) 导出每组细胞 ID 列表（CSV + TXT），以及一个汇总表 ----
# # # meta <- cycling_clean@meta.data
# # # meta$cell_id <- rownames(meta)
# # # meta$group   <- as.character(meta[[group_col]])
# # #
# # # # 6.1 汇总：每组细胞数
# # # count_df <- as.data.frame(table(meta$group), stringsAsFactors = FALSE)
# # # colnames(count_df) <- c("group", "n_cells")
# # # write.csv(count_df, file.path(out_post, "cycling_clean_cells_per_group_summary.csv"), row.names = FALSE)
# # #
# # # # 6.2 每组 cell list：CSV（两列：group, cell_id）
# # # write.csv(meta[, c("group","cell_id")], file.path(out_post, "cycling_clean_cellIDs_long.csv"), row.names = FALSE)
# # #
# # # # 6.3 每组单独一个 TXT（写方法/补充材料最方便）
# # # groups_here <- sort(unique(meta$group))
# # # for (g in groups_here) {
# # #   ids <- meta$cell_id[meta$group == g]
# # #   fn  <- file.path(out_post, paste0("cycling_clean_cellIDs_", g, ".txt"))
# # #   writeLines(ids, con = fn)
# # # }
# # #
# # # cat("POSTCHECK outputs saved under: ", out_post, "\n", sep = "")
# # # cat("==== Step 1.7 POSTCHECK: cycling_clean validation + cellID export DONE ====\n")
# # #
# # # cat("==== Step 1.7 【主分析cyclingTreg细胞】结束 \n")
# # # ============================================================主分析cyclingTreg细胞结束
# #
# # ============================================================step1.8_主分析eTreg细胞开始
# # ============================================================
# # Step 1.7 (PUBLISHABLE, eTreg-like):
# #   - load treg_main_labeled_step1.2.rds (with Treg_subclass_cluster)
# #   - Treg core gate
# #   - extract eTreg-like
# #   - publishable cleaning (FIX: B/NK false-kill using STRONG rules)
# #   - pseudo-bulk logFC ranking (no p-values) with gene filters
# #   - single-cell validation plots (FeaturePlot/VlnPlot)
# #   - export per-group cell IDs
# # ============================================================
#
# cat("==== Step 1.7 (eTreg-like PUBLISHABLE): START ====\n")
#
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(Matrix)
#   library(ggplot2)
# })
#
# dir.create("./temp", showWarnings = FALSE, recursive = TRUE)
#
# # ---------- 0) 显式用 base ----------
# base_intersect <- base::intersect
# base_setdiff   <- base::setdiff
# base_unname    <- base::unname
#
# # ---------- 1) 读取带 publishable label 的对象 ----------
# treg_rds <- "./temp/treg_main_labeled_step1.2.rds"
# if (!file.exists(treg_rds)) stop("找不到：", treg_rds, "\n请先完成 Step1.2 并 saveRDS。")
#
# treg_main <- readRDS(treg_rds)
# DefaultAssay(treg_main) <- "RNA"
# cat(">> Loaded:", treg_rds, "\n")
#
# # ---------- 2) 列名（按你当前项目口径） ----------
# group_col    <- "group"
# subclass_col <- "Treg_subclass_cluster"
#
# if (!(group_col %in% colnames(treg_main@meta.data))) {
#   stop("meta.data 里找不到 group 列：", group_col,
#        "\n请改成你的真实列名（例如 orig.ident）。")
# }
# if (!(subclass_col %in% colnames(treg_main@meta.data))) {
#   stop("meta.data 里找不到 subclass 列：", subclass_col,
#        "\n请确认你已运行 Step1.1/Step1.2 并保存 treg_main_labeled_step1.2.rds。")
# }
#
# # ---------- 3) Treg core gate（保证仍在 Treg 主谱系内） ----------
# need_core <- c("Foxp3","Il2ra","Trac","Cd3e")
# need_core <- base_intersect(need_core, rownames(treg_main))
# if (length(need_core) < 3) {
#   stop("关键 Treg core 基因缺失过多（Foxp3/Il2ra/Trac/Cd3e）。请检查 gene symbol 命名体系。")
# }
# treg_core <- subset(treg_main, subset = Foxp3 > 0 & Il2ra > 0 & Trac > 0 & Cd3e > 0)
# cat("Treg core cells:", ncol(treg_core), "\n")
# cat("Treg core per group:\n"); print(table(as.character(treg_core@meta.data[[group_col]])))
#
# # ---------- 4) 提取 eTreg-like（用 cells 向量，最稳） ----------
# meta_core <- treg_core@meta.data
# cells_etreg <- rownames(meta_core)[as.character(meta_core[[subclass_col]]) == "eTreg-like"]
#
# cat("eTreg-like cells (core, before clean):", length(cells_etreg), "\n")
# if (length(cells_etreg) == 0) {
#   cat("Available subclass labels:\n")
#   print(sort(table(as.character(meta_core[[subclass_col]])), decreasing = TRUE))
#   stop("未找到 eTreg-like。请检查 subclass label 字符串是否一致（空格/大小写）。")
# }
#
# etreg <- subset(treg_core, cells = cells_etreg)
# cat("eTreg-like object cells:", ncol(etreg), "\n")
# cat("eTreg-like cells per group:\n")
# print(sort(table(as.character(etreg@meta.data[[group_col]])), decreasing = TRUE))
#
# # ---------- 5) publishable clean（FIX B/NK false-kill：用 STRONG 规则，不用 any-hit） ----------
# thr_core <- 0.25  # log-normalized data threshold
# thr_bad  <- 0.25
#
# # marker panels（尽量保持与你 cycling 的口径一致，便于补充材料统一）
# genes_Tcore    <- c("Trac","Cd3e","Cd3d")
# genes_Tregcore <- c("Foxp3","Il2ra","Ikzf2","Ctla4")
#
# genes_myeloid  <- c("Tyrobp","Lyz2","Csf1r","LST1","Itgam","S100a8","S100a9","Aif1","Fcgr3","Lgals3")
# genes_B        <- c("Ms4a1","Cd79a","Cd79b","Cd74","Igkc","Jchain","Ighm")
#
# # 关键 FIX：NK 列表去掉 Trdc/Trgc1（容易把“像T的细胞”误判成 NK），只保留更典型 NK/CTL 轴
# genes_NK       <- c("Ncr1","Klrb1c","Klrk1","Prf1","Gzmb")
#
# genes_RBC      <- c("Hbb-bs","Hbb-bt","Hba-a1","Hba-a2","Alas2")
#
# genes_epi   <- c("Epcam","Tacstd2","Krt8","Krt18","Krt19","Krt14","Krt17")
# genes_endo  <- c("Pecam1","Kdr","Vwf","Eng","Ramp2")
# genes_fibro <- c("Col1a1","Col1a2","Dcn","Lum","Col3a1","Pdgfra")
# genes_peri  <- c("Acta2","Rgs5","Pdgfrb","Des","Cspg4","Tagln")
#
# genes_melanoma <- c("Sox10","Pmel","Dct","Mlana","Tyrp1","Slc45a2","Trpm1","Slc24a5","Mitf","Tyr")
#
# panel_all <- unique(c(
#   genes_Tcore, genes_Tregcore,
#   genes_myeloid, genes_B, genes_NK, genes_RBC,
#   genes_epi, genes_endo, genes_fibro, genes_peri,
#   genes_melanoma
# ))
# panel_present <- base_intersect(panel_all, rownames(etreg))
# cat("Panel requested:", length(panel_all), " | present:", length(panel_present), "\n")
#
# if (length(panel_present) < 10) {
#   stop("panel_present 太少（<10）。说明 gene 命名可能不是 symbol（例如 Ensembl）。请先统一命名体系。")
# }
#
# df <- FetchData(etreg, vars = panel_present)
#
# hit_any <- function(gset, thr = thr_bad) {
#   g2 <- base_intersect(gset, colnames(df))
#   if (length(g2) == 0) return(rep(FALSE, nrow(df)))
#   out <- rep(FALSE, nrow(df))
#   for (g in g2) out <- out | (df[[g]] > thr)
#   out
# }
#
# hit_ge2 <- function(gset, thr = thr_bad) {
#   g2 <- base_intersect(gset, colnames(df))
#   if (length(g2) == 0) return(rep(FALSE, nrow(df)))
#   k <- rep(0L, nrow(df))
#   for (g in g2) k <- k + (df[[g]] > thr)
#   (k >= 2)
# }
#
# # core / Treg identity flags
# flag_T <- hit_any(genes_Tcore, thr = thr_core)
#
# flag_Foxp3 <- if ("Foxp3" %in% colnames(df)) (df$Foxp3 > thr_core) else rep(FALSE, nrow(df))
# flag_Il2ra <- if ("Il2ra" %in% colnames(df)) (df$Il2ra > thr_core) else rep(FALSE, nrow(df))
# flag_Treg_strict <- flag_Foxp3 & flag_Il2ra
#
# # contaminations (STRONG rules)
# flag_my_strong <- hit_ge2(genes_myeloid, thr = thr_bad)
#
# # B/NK false-kill FIX：只用 STRONG（>=2），不做 any-hit
# flag_B_strong  <- hit_ge2(genes_B,  thr = thr_bad)
# flag_NK_strong <- hit_ge2(genes_NK, thr = thr_bad)
#
# flag_RBC_any <- hit_any(genes_RBC, thr = thr_bad)
#
# # 非免疫：用 STRONG（>=2）更稳，避免因单个背景基因误杀
# nonimmune_gset <- unique(c(genes_epi, genes_endo, genes_fibro, genes_peri))
# flag_nonimmune_strong <- hit_ge2(nonimmune_gset, thr = thr_bad)
#
# # 黑色素细胞/肿瘤：用 STRONG（>=2）避免单基因背景触发
# flag_mel_strong <- hit_ge2(genes_melanoma, thr = thr_bad)
#
# # ---------- 5.1) 输出 QC hit-rate（按组，便于定位是哪条规则在“杀”某组） ----------
# grp_all <- as.character(etreg@meta.data[[group_col]])
# names(grp_all) <- rownames(etreg@meta.data)
# grp <- grp_all[rownames(df)]
#
# report_rate <- function(flag, name) {
#   tab <- tapply(flag, grp, mean)
#   cat("\nRule hit-rate: ", name, "\n", sep = "")
#   print(tab)
# }
#
# cat("\nQC hit-rates by group (eTreg-like, before clean):\n")
# cat("eTreg-like per group:\n"); print(table(grp))
#
# report_rate(flag_Treg_strict, "Treg strict (Foxp3&Il2ra)")
# report_rate(flag_my_strong, "Myeloid STRONG (>=2)")
# report_rate(flag_B_strong,  "B STRONG (>=2)")
# report_rate(flag_NK_strong, "NK STRONG (>=2; no Trdc/Trgc1)")
# report_rate(flag_RBC_any,   "RBC any")
# report_rate(flag_nonimmune_strong, "Non-immune STRONG (>=2)")
# report_rate(flag_mel_strong, "Melanoma STRONG (>=2)")
#
# # ---------- 5.2) publishable keep ----------
# keep <- flag_T & flag_Treg_strict &
#   (!flag_my_strong) &
#   (!flag_B_strong) & (!flag_NK_strong) &
#   (!flag_RBC_any) &
#   (!flag_nonimmune_strong) &
#   (!flag_mel_strong)
#
# keep_cells <- rownames(df)[keep]
# etreg_clean <- subset(etreg, cells = keep_cells)
#
# cat("\neTreg_clean cells:", ncol(etreg_clean), "\n")
# cat("eTreg_clean cells per group:\n")
# print(sort(table(as.character(etreg_clean@meta.data[[group_col]])), decreasing = TRUE))
#
# etreg_clean_rds <- "./temp/eTreg_like_clean_for_pseudobulk_publishable.rds"
# saveRDS(etreg_clean, etreg_clean_rds)
# cat(">> Saved eTreg_clean RDS:", etreg_clean_rds, "\n")
#
# # ---------- 6) 输出目录 ----------
# outdir <- "./temp/Step1.7_eTreg_publishable"
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
#
# # ---------- 7) 导出每组 cell IDs（Methods/补充材料/排查） ----------
# grp2 <- as.character(etreg_clean@meta.data[[group_col]])
# names(grp2) <- colnames(etreg_clean)
#
# write.csv(
#   data.frame(cell = names(grp2), group = grp2, stringsAsFactors = FALSE),
#   file = file.path(outdir, "eTreg_clean_cellIDs_all.csv"),
#   row.names = FALSE
# )
# for (g in sort(unique(grp2))) {
#   cells_g <- names(grp2)[grp2 == g]
#   write.csv(
#     data.frame(cell = cells_g, stringsAsFactors = FALSE),
#     file = file.path(outdir, paste0("eTreg_clean_cellIDs_", g, ".csv")),
#     row.names = FALSE
#   )
# }
#
# # ---------- 8) 单细胞层面验证图：不用 Seurat::VlnPlot，手写 ggplot 版（规避 S4SXP bug） ----------
# val_genes <- c(
#   "Foxp3","Il2ra","Ctla4","Il2rb",
#   "Icos","Tnfrsf4","Tnfrsf18","Batf","Maf",
#   "Pdcd1","Tigit","Lag3",
#   "Mki67","Top2a"
# )
#
# val_genes <- base::intersect(val_genes, rownames(etreg_clean))
# cat("Validation genes present:", paste(val_genes, collapse = ", "), "\n")
#
# if (length(val_genes) > 0) {
#
#   # 1) 抓取表达 + 分组
#   vars_need <- c(group_col, val_genes)
#   df_v <- Seurat::FetchData(etreg_clean, vars = vars_need)
#
#   # 2) 确保分组是字符/因子，并按你想要的顺序（可选）
#   df_v[[group_col]] <- as.character(df_v[[group_col]])
#   # 如你固定三组顺序（可按需改）
#   df_v[[group_col]] <- factor(df_v[[group_col]], levels = c("control1","s1_hIL","s40kD"))
#
#   # 3) 转成长表：用 utils::stack（显式命名空间，避免 stack 不可见）
# df_expr <- utils::stack(df_v[, val_genes, drop = FALSE])
#
# df_long <- data.frame(
#   group = rep(df_v[[group_col]], times = length(val_genes)),
#   gene  = as.character(df_expr$ind),
#   expr  = as.numeric(df_expr$values),
#   stringsAsFactors = FALSE
# )
#   # 修正 group：每个基因都对应一遍所有细胞（stack 的顺序是按列堆叠）
#   df_long$group <- rep(df_v[[group_col]], times = length(val_genes))
#
#   # 4) 画图：violin + 可选箱线（更像 Seurat 风格）
#   p_vln_manual <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = expr)) +
#     ggplot2::geom_violin(trim = TRUE, scale = "width") +
#     ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.2) +
#     ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 5) +
#     ggplot2::theme_bw() +
#     ggplot2::labs(
#       title = "eTreg_clean | Key genes by group (manual violin; Seurat VlnPlot bypass)",
#       x = "Group", y = "Expression (log-normalized data)"
#     ) +
#     ggplot2::theme(
#       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
#       strip.text  = ggplot2::element_text(size = 9)
#     )
#
#   # 5) 保存（这里用 ggsave 一般不会触发你之前那个 bug；如果仍担心，也可改 pdf+print）
#   ggplot2::ggsave(
#     filename = file.path(outdir, "eTreg_clean_validation_VlnPlot_by_group.pdf"),
#     plot     = p_vln_manual,
#     width    = 14, height = 8
#   )
#
# } else {
#   warning("验证基因在对象里一个都找不到（可能是命名体系不一致）。跳过绘图。")
# }
# # ---------- 9) pseudo-bulk logFC（不输出 p 值；描述性排名） ----------
# counts <- LayerData(etreg_clean, assay = "RNA", layer = "counts")
# if (is.null(counts)) stop("找不到 RNA counts layer（LayerData返回NULL）。")
#
# tab <- table(grp2)
# cat("\nCells per group (eTreg_clean):\n"); print(tab)
#
# # ---------- 固定阈值纳入：s40kD 只要 n>=20 就纳入；其他组不再被“更高阈值”踢掉 ----------
# wanted_groups <- c("control1", "s1_hIL", "s40kD")
# tab <- tab[names(tab) %in% wanted_groups]  # 只看三组
# cat("\nCells per group (wanted groups only):\n"); print(tab)
#
# # 1) 先把实际存在的组都纳入（不做 50/30/20 的阶梯淘汰）
# ok_groups <- names(tab)
#
# # 2) s40kD 的“硬阈值”：<20 就剔除；>=20 强制保留
# min_cells_s40kD <- 20
# if ("s40kD" %in% ok_groups) {
#   if (as.integer(tab["s40kD"]) < min_cells_s40kD) {
#     ok_groups <- setdiff(ok_groups, "s40kD")
#     cat("NOTE: s40kD dropped because n < ", min_cells_s40kD,
#         " (n=", as.integer(tab["s40kD"]), ").\n", sep = "")
#   } else {
#     cat("NOTE: s40kD kept (n >= ", min_cells_s40kD,
#         "; n=", as.integer(tab["s40kD"]), ").\n", sep = "")
#   }
# }
#
# # 3) 最低要求：至少两组才能做 pairwise logFC
# if (length(ok_groups) < 2) {
#   stop("可用组别不足2个，无法做 pseudo-bulk logFC。\n当前 tab = ",
#        paste(names(tab), as.integer(tab), sep="=", collapse=", "))
# }
#
# # 4) 提醒：如果某组很小（例如 control1=16），结果只能作为“描述性/补充材料”
# small_n <- ok_groups[as.integer(tab[ok_groups]) < 20]
# if (length(small_n) > 0) {
#   cat("WARNING: Some groups have <20 cells: ",
#       paste0(small_n, "=", as.integer(tab[small_n]), collapse=", "),
#       "\nThese comparisons should be treated as descriptive/supplementary.\n", sep = "")
# }
#
# cat("Groups kept (fixed-inclusion logic): ", paste(ok_groups, collapse = ", "), "\n", sep = "")
#
# pb_counts <- sapply(ok_groups, function(g) {
#   cells <- names(grp2)[grp2 == g]
#   Matrix::rowSums(counts[, cells, drop = FALSE])
# })
# rownames(pb_counts) <- rownames(counts)
#
# libsize <- colSums(pb_counts)
# cat("Pseudo-bulk library sizes:\n"); print(libsize)
#
# cpm <- sweep(pb_counts, 2, libsize, FUN = "/") * 1e6
# logcpm <- log2(cpm + 1)
#
# write.csv(pb_counts, file.path(outdir, "pseudobulk_counts_gene_x_group.csv"))
# write.csv(logcpm,   file.path(outdir, "pseudobulk_logCPM_gene_x_group.csv"))
#
# # ---------- 9.1) pseudo-bulk counts ----------
# pb_counts <- sapply(ok_groups, function(g) {
#   cells <- names(grp2)[grp2 == g]
#   Matrix::rowSums(counts[, cells, drop = FALSE])
# })
# rownames(pb_counts) <- rownames(counts)
#
# libsize <- colSums(pb_counts)
# cat("\nPseudo-bulk library sizes:\n"); print(libsize)
#
# # ---------- 9.2) logCPM ----------
# cpm <- sweep(pb_counts, 2, libsize, FUN = "/") * 1e6
# logcpm <- log2(cpm + 1)
#
# write.csv(pb_counts, file.path(outdir, "pseudobulk_counts_gene_x_group.csv"))
# write.csv(logcpm,   file.path(outdir, "pseudobulk_logCPM_gene_x_group.csv"))
#
# # ---------- 9.3) gene 过滤：避免极低 counts 冲榜（可按稀疏程度调） ----------
# min_max_counts <- 10
# min_sum_counts <- 20
# min_aveexpr    <- 1.0
#
# # 小样本组（例如 s40kD 只有 20）更容易受低counts影响：
# # 如果 ok_groups 里包含 s40kD 且它刚好在阈值边缘，建议稍微加严 min_sum_counts/min_aveexpr
# if ("s40kD" %in% ok_groups && as.integer(tab["s40kD"]) <= 25) {
#   cat("\n[NOTE] s40kD is small (<=25 cells). Tightening gene filters for robustness.\n")
#   min_sum_counts <- 30
#   min_aveexpr    <- 1.2
# }
#
# rank_one_pair <- function(g1, g2) {
#   c1 <- pb_counts[, g1]; c2 <- pb_counts[, g2]
#   AveExpr <- rowMeans(logcpm[, c(g1, g2), drop = FALSE])
#
#   keep_gene <- (pmax(c1, c2) >= min_max_counts) &
#     ((c1 + c2) >= min_sum_counts) &
#     (AveExpr >= min_aveexpr)
#
#   df_rank <- data.frame(
#     gene = rownames(logcpm),
#     logFC = as.numeric(logcpm[, g1] - logcpm[, g2]),
#     AveExpr = as.numeric(AveExpr),
#     c1 = as.integer(c1),
#     c2 = as.integer(c2),
#     stringsAsFactors = FALSE
#   )
#   df_rank <- df_rank[keep_gene, , drop = FALSE]
#   df_rank <- df_rank[order(df_rank$logFC, decreasing = TRUE), , drop = FALSE]
#
#   tag <- paste0(g1, "_vs_", g2)
#   write.csv(df_rank, file.path(outdir, paste0("logFC_rank_", tag, ".csv")), row.names = FALSE)
#   write.csv(head(df_rank, 200), file.path(outdir, paste0("logFC_rank_", tag, "_TOP200.csv")), row.names = FALSE)
#   write.csv(tail(df_rank, 200), file.path(outdir, paste0("logFC_rank_", tag, "_BOTTOM200.csv")), row.names = FALSE)
#
#   # Results-friendly plots
#   if (nrow(df_rank) > 0) {
#     topn <- min(20, nrow(df_rank))
#     df_up <- head(df_rank, topn)
#     df_dn <- tail(df_rank, topn)
#     df_dn <- df_dn[order(df_dn$logFC, decreasing = FALSE), ]
#
#     df_up$dir <- "Up"; df_dn$dir <- "Down"
#     df_plot <- rbind(df_up, df_dn)
#     df_plot$gene <- factor(df_plot$gene, levels = df_plot$gene)
#
#     p_bar <- ggplot(df_plot, aes(x = gene, y = logFC)) +
#       geom_col() + coord_flip() + facet_wrap(~dir, scales = "free_y") +
#       theme_bw() + ggtitle(paste0("eTreg_clean | ", tag, " | Top logFC genes (no p-values)"))
#     ggsave(file.path(outdir, paste0("TopGenes_Barplot_", tag, ".pdf")), p_bar, width = 10, height = 8)
#
#     p_ma <- ggplot(df_rank, aes(x = AveExpr, y = logFC)) +
#       geom_point(size = 0.5, alpha = 0.5) +
#       geom_hline(yintercept = 0, linetype = 2) +
#       theme_bw() + ggtitle(paste0("eTreg_clean | ", tag, " | MA plot (logCPM)"))
#     ggsave(file.path(outdir, paste0("MAplot_", tag, ".pdf")), p_ma, width = 8, height = 6)
#   }
#
#   cat("Done:", tag, " genes kept =", nrow(df_rank), "\n")
# }
#
# pairs <- combn(ok_groups, 2, simplify = FALSE)
# cat("\nPairwise comparisons to run:\n")
# print(pairs)
#
# for (pp in pairs) rank_one_pair(pp[1], pp[2])
#
# cat("\nOutputs saved under: ", outdir, "\n", sep="")
# cat("==== Step 1.7 (eTreg-like PUBLISHABLE): DONE ====\n")
# #  # # ---------- 下一步了回答“比例减少”拆成 A/B/C 三类机制（trafficking/retention、proliferation、stress/apoptosis + IL2/STAT5 axis）做描述性证据链
# # # ============================================================
# # # ============================================================
# # ============================================================
# # [MOD] Step 1.7-CLOSURE-EXT (Minimal closure; no full rerun)
# #   目标：在不重跑全流程的前提下，快速回答：
# #     (1) s40kD 为什么 Treg 比例减少？
# #     (2) s40kD vs s1_hIL 对 Treg 作用为何不同？
# #
# #   最短闭环输出：
# #     A) coreTreg / treg_main 分母比例（descriptive）
# #     A2) coreTreg / 全样本所有细胞（全局分母）比例（descriptive）   <-- [ADD]
# #     B) eTreg-like / cycling / rTreg-ish 结构是否重排（within coreTreg 的构成）
# #     C) 机制拆解（3类模块 + IL-2 轴模块）
# #
# #   说明：本段代码依赖 treg_main_labeled_step1.2.rds；
# #        若要算 A2（全局分母），会自动尝试读取 ./temp/combined_* 全细胞对象
# # ============================================================
#
# cat("\n==== [MOD] Step 1.7-CLOSURE-EXT (Minimal closure; no full rerun) START ====\n")
#
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(Matrix)
#   library(ggplot2)
# })
#
# # ---------- [MOD-FIX] 0) 显式 base/utils，规避 conflicted & write.csv not found ----------
# base_intersect <- base::intersect
# base_setdiff   <- base::setdiff
# base_unname    <- base::unname
# base_table     <- base::table
# base_names     <- base::names
# base_match     <- base::match
#
# outdir_closure <- "./temp/Step1.7_CLOSURE_EXT"
# dir.create(outdir_closure, showWarnings = FALSE, recursive = TRUE)
# #
# # ---------- 1) 读取对象（Step1.2 已保存） ----------
# treg_rds <- "./temp/treg_main_labeled_step1.2.rds"
# if (!file.exists(treg_rds)) stop("找不到：", treg_rds, "。请先完成 Step1.2 并保存该 RDS。")
# treg_main <- readRDS(treg_rds)
# DefaultAssay(treg_main) <- "RNA"
# cat(">> Loaded:", treg_rds, "\n")
#
# group_col    <- "group"
# subclass_col <- "Treg_subclass_cluster"
# if (!(group_col %in% colnames(treg_main@meta.data))) stop("treg_main@meta.data 缺少 group 列：", group_col)
# if (!(subclass_col %in% colnames(treg_main@meta.data))) stop("treg_main@meta.data 缺少 subclass 列：", subclass_col)
#
# wanted_groups <- c("control1","s1_hIL","s40kD")
#
# # ---------- 2) Node-0：treg_main 规模 ----------
# cat("\n[Node-0] treg_main cells per group:\n")
# print(base_table(as.character(treg_main@meta.data[[group_col]])))
#
# # ---------- 3) Node-1：treg_core（与你 Step1.7 一致） ----------
# need_core <- base_intersect(c("Foxp3","Il2ra","Trac","Cd3e"), rownames(treg_main))
# if (length(need_core) < 3) stop("关键 Treg core 基因缺失过多（Foxp3/Il2ra/Trac/Cd3e）。请检查基因命名体系。")
#
# treg_core <- subset(treg_main, subset = Foxp3 > 0 & Il2ra > 0 & Trac > 0 & Cd3e > 0)
# cat("\n[Node-1] treg_core cells:", ncol(treg_core), "\n")
# cat("[Node-1] treg_core per group:\n")
# print(base_table(as.character(treg_core@meta.data[[group_col]])))
#
# # ---------- [A] Minimal table：coreTreg fraction within treg_main（你之前已做过的A表） ----------
# tab_main <- base_table(as.character(treg_main@meta.data[[group_col]]))
# tab_core <- base_table(as.character(treg_core@meta.data[[group_col]]))
# wanted <- base_intersect(base_names(tab_main), base_names(tab_core))
# wanted <- wanted[wanted %in% wanted_groups]
#
# df_core_within_treg_main <- data.frame(
#   group = wanted,
#   n_treg_main = as.integer(tab_main[wanted]),
#   n_treg_core = as.integer(tab_core[wanted]),
#   frac_core_within_treg_main = as.numeric(tab_core[wanted]) / as.numeric(tab_main[wanted]),
#   stringsAsFactors = FALSE
# )
# df_core_within_treg_main <- df_core_within_treg_main[base_match(wanted_groups, df_core_within_treg_main$group), , drop = FALSE]
# cat("\n[Node-1+] coreTreg fraction within treg_main:\n")
# print(df_core_within_treg_main)
#
# out_A <- file.path(outdir_closure, "coreTreg_fraction_within_treg_main.csv")
# utils::write.csv(df_core_within_treg_main, out_A, row.names = FALSE)
# cat("[Node-1+] Wrote:", out_A, "\n")
#
# # ---------- [ADD A2] Minimal table：coreTreg / 全样本所有细胞（全局分母） ----------
# # 说明：需要 combined 全对象；若内存没有，则自动尝试读取 ./temp/combined_* 缓存
# cat("\n[Node-1++] coreTreg fraction within ALL CELLS (global denominator):\n")
#
# detect_group_col <- function(md) {
#   cand <- c("group", "orig.ident", "sample", "Sample", "condition", "treatment")
#   hit <- cand[cand %in% colnames(md)]
#   if (length(hit) == 0) return(NA_character_)
#   hit[1]
# }
#
# if (!exists("combined")) {
#   cand_rds <- c(
#     "./temp/combined_after_marker_6.rds",
#     "./temp/combined_umap_5.rds",
#     "./temp/combined_normalization_4.rds",
#     "./temp/combined_raw_data_1.rds"
#   )
#   cand_rds <- cand_rds[file.exists(cand_rds)]
#
#   if (length(cand_rds) == 0) {
#     cat("[WARN] 未找到 ./temp/combined_* RDS，因此无法输出 A2（coreTreg/全局分母）表。\n")
#     cat("       你若需要 A2，请确保保存了 combined_after_marker_6.rds 或 combined_umap_5.rds。\n")
#   } else {
#     combined <- readRDS(cand_rds[1])
#     cat(">> Loaded combined from:", cand_rds[1], "\n")
#   }
# } else {
#   cat(">> Using combined from memory.\n")
# }
#
# if (exists("combined")) {
#   group_col_all <- detect_group_col(combined@meta.data)
#   if (is.na(group_col_all)) {
#     stop("combined@meta.data 中找不到 group/orig.ident/sample 等列，无法计算全局分母。")
#   }
#
#   tab_all <- base_table(as.character(combined@meta.data[[group_col_all]]))
#   tab_all <- tab_all[base_names(tab_all) %in% wanted_groups]
#
#   tab_core2 <- base_table(as.character(treg_core@meta.data[[group_col]]))
#   tab_core2 <- tab_core2[base_names(tab_core2) %in% wanted_groups]
#
#   wanted2 <- base_intersect(base_names(tab_all), base_names(tab_core2))
#   df_core_within_allcells <- data.frame(
#     group = wanted2,
#     n_all_cells = as.integer(tab_all[wanted2]),
#     n_coreTreg  = as.integer(tab_core2[wanted2]),
#     frac_coreTreg_within_allcells = as.numeric(tab_core2[wanted2]) / as.numeric(tab_all[wanted2]),
#     stringsAsFactors = FALSE
#   )
#   df_core_within_allcells <- df_core_within_allcells[base_match(wanted_groups, df_core_within_allcells$group), , drop = FALSE]
#
#   print(df_core_within_allcells)
#
#   out_A2 <- file.path(outdir_closure, "coreTreg_fraction_within_allcells.csv")
#   utils::write.csv(df_core_within_allcells, out_A2, row.names = FALSE)
#   cat("[Node-1++] Wrote:", out_A2, "\n")
# }
#
# # ---------- [B] Node-2：coreTreg 内部 subclasses 构成（counts + proportions） ----------
# meta_core <- treg_core@meta.data
# meta_core$grp <- as.character(meta_core[[group_col]])
# meta_core$sub <- as.character(meta_core[[subclass_col]])
#
# keep_idx <- meta_core$grp %in% wanted_groups & !is.na(meta_core$sub)
# meta_core2 <- meta_core[keep_idx, , drop = FALSE]
#
# tab_sub_count <- base_table(meta_core2$sub, meta_core2$grp)
# cat("\n[Node-2] Treg subclass counts (within treg_core):\n")
# print(tab_sub_count)
#
# tab_sub_prop <- sweep(tab_sub_count, 2, colSums(tab_sub_count), FUN = "/")
# cat("\n[Node-2] Treg subclass proportions (within treg_core):\n")
# print(round(tab_sub_prop, 4))
#
# utils::write.csv(as.data.frame.matrix(tab_sub_count), file.path(outdir_closure, "treg_core_subclass_counts.csv"))
# utils::write.csv(as.data.frame.matrix(tab_sub_prop),  file.path(outdir_closure, "treg_core_subclass_proportions.csv"))
#
# # ---------- [C] Node-3：机制拆解模块 + IL-2 轴模块（AddModuleScore + VlnPlot） ----------
# score_gene_sets <- list(
#   Prolif_CC = c("Mki67","Top2a","Mcm5","Mcm6","Mcm7","Pcna","Tubb5","Cenpf","Cdk1","Birc5"),
#   Apoptosis = c("Bax","Bak1","Bcl2","Bcl2l1","Mcl1","Casp3","Casp8","Fas","Tradd","Tnfrsf1a"),
#   Migration = c("Ccr7","Ccr4","Ccr8","Cxcr3","Cxcr4","Sell","Itga4","Itgb1","Icam1","S1pr1"),
#   IL2_axis  = c("Il2ra","Il2rb","Il2rg","Stat5a","Stat5b","Socs1","Socs3","Cish")
# )
#
# score_gene_sets_present <- lapply(score_gene_sets, function(gs) base_intersect(gs, rownames(treg_core)))
# lens <- sapply(score_gene_sets_present, length)
# score_gene_sets_present <- score_gene_sets_present[lens >= 3]
#
# if (length(score_gene_sets_present) == 0) {
#   stop("模块基因集在对象中匹配过少（全部<3）。请检查 gene symbol 命名体系。")
# }
#
# cat("\n[Node-3] Module genes present per set:\n")
# print(sapply(score_gene_sets_present, length))
#
# treg_core_scored <- AddModuleScore(
#   object   = treg_core,
#   features = base_unname(score_gene_sets_present),
#   name     = paste0(names(score_gene_sets_present), "_MS"),
#   assay    = "RNA",
#   search   = FALSE
# )
#
# # ---- 自动匹配 AddModuleScore 生成的列名（避免 *_MS1 不存在）----
# score_prefixes <- paste0(names(score_gene_sets_present), "_MS")
# find_ms_col <- function(meta_cols, prefix) {
#   hit <- grep(paste0("^", prefix, "\\d+$"), meta_cols, value = TRUE)
#   if (length(hit) == 0) return(NA_character_)
#   hit[length(hit)]
# }
# meta_cols <- colnames(treg_core_scored@meta.data)
# score_cols <- vapply(score_prefixes, function(pf) find_ms_col(meta_cols, pf), FUN.VALUE = character(1))
#
# if (any(is.na(score_cols))) {
#   miss <- score_prefixes[is.na(score_cols)]
#   cat("\n[ERROR] AddModuleScore 后未找到以下模块 score 列：\n"); print(miss)
#   cat("\n[DEBUG] meta.data 中匹配 *_MS\\d+$ 的列：\n"); print(grep("_MS\\d+$", meta_cols, value = TRUE))
#   stop("AddModuleScore 模块分数列名匹配失败。")
# }
#
# cat("\n[OK] Matched module score columns:\n")
# print(data.frame(set = names(score_gene_sets_present), prefix = score_prefixes, col = score_cols, row.names = NULL))
#
# # ---- 聚焦：eTreg-like + cycling + rTreg（用来解释“比例减少”拆解）----
# sub_focus <- c("eTreg-like","cycling Treg","rTreg-ish / transitional")
# meta_s <- treg_core_scored@meta.data
#
# cells_focus <- rownames(meta_s)[
#   as.character(meta_s[[subclass_col]]) %in% sub_focus &
#   as.character(meta_s[[group_col]]) %in% wanted_groups
# ]
# treg_focus <- subset(treg_core_scored, cells = cells_focus)
#
# cat("\n[Node-4] treg_focus cells:", ncol(treg_focus), "\n")
# cat("[Node-4] treg_focus per group:\n");    print(base_table(as.character(treg_focus@meta.data[[group_col]])))
# cat("[Node-4] treg_focus per subclass:\n"); print(base_table(as.character(treg_focus@meta.data[[subclass_col]])))
#
# # ============================================================
# # [EXT++/FUNC] coreTreg functional state scoring (single-cell; stratified by subclass)
# #   Goal: Compare functional states across groups within each coreTreg subclass
# #         (cycling Treg / eTreg-like / rTreg-ish)
# #   Method: score = per-cell mean(log-normalized expr) over gene set
# #   Output (PDF only): ./temp/Step1.7_CLOSURE_EXT/SC_coreTreg_FunctionScores/
# # ============================================================
#
# cat("\n[EXT++/FUNC] coreTreg functional scoring (subclass-stratified; PDF only) ...\n")
#
# out_func <- file.path(outdir_closure, "SC_coreTreg_FunctionScores")
# dir.create(out_func, showWarnings = FALSE, recursive = TRUE)
#
# # --- Factor order (keep stable for figures) ---
# subclass_levels <- c("cycling Treg", "eTreg-like", "rTreg-ish / transitional")
# group_levels <- wanted_groups
#
# # --- Define functional gene sets (Treg-suppressive vs activation/inflammation axes) ---
# # Notes:
# # 1) "Suppressive" emphasizes canonical suppressive/checkpoint machinery.
# # 2) "EffectorTreg" emphasizes eTreg activation/co-stim receptors / TFs.
# # 3) "Inflamm_NFkB" and "IFN_Response" capture inflammatory activation that can modulate Treg stability/function.
# # 4) These sets are intentionally compact and interpretable for figures.
# func_gene_sets <- list(
#   Treg_Suppressive = c("Foxp3","Il2ra","Ctla4","Tigit","Lag3","Pdcd1","Icos","Entpd1","Nt5e","Il10","Tgfb1","Ebi3","Prdm1","Gzmb","Areg"),
#   Treg_EffectorTreg = c("Foxp3","Ikzf2","Icos","Ctla4","Tnfrsf4","Tnfrsf18","Batf","Maf","Tigit"),
#   Inflamm_NFkB = c("Nfkb1","Nfkb2","Rel","Rela","Relb","Tnf","Tnfaip3","Nfkbia","Nfkbie","Il1b","Il6","Cxcl10","Ccl2","Ccl3","Ccl4"),
#   IFN_Response = c("Stat1","Irf1","Irf7","Isg15","Ifit1","Ifit2","Ifit3","Oas1a","Oas1b","Mx1","Cxcl9","Cxcl10"),
#   IL2_STAT5_Downstream = c("Il2ra","Il2rb","Il2rg","Stat5a","Stat5b","Cish","Socs1","Socs3","Bcl2","Mcl1")
# )
#
# # --- Match genes ---
# func_gene_sets_present <- lapply(func_gene_sets, function(gs) base::intersect(gs, rownames(treg_focus)))
# lens <- vapply(func_gene_sets_present, length, integer(1))
# func_gene_sets_present <- func_gene_sets_present[lens >= 5]  # keep sets with enough overlap
#
# if (length(func_gene_sets_present) == 0) {
#   cat("[EXT++/FUNC] No functional gene sets passed overlap threshold (>=5). Skip.\n")
# } else {
#
#   cat("[EXT++/FUNC] Functional sets kept (>=5 matched genes):\n")
#   print(data.frame(
#     set = names(func_gene_sets_present),
#     n_genes = vapply(func_gene_sets_present, length, integer(1)),
#     stringsAsFactors = FALSE
#   ))
#
#   # --- Pull expression matrix once (log-normalized data) ---
#   expr_mat <- Seurat::GetAssayData(treg_focus, assay = "RNA", slot = "data")
#
#   # --- Compute per-cell scores (mean expression across genes) ---
#   score_mat <- sapply(names(func_gene_sets_present), function(nm) {
#     gs <- func_gene_sets_present[[nm]]
#     if (length(gs) == 1) {
#       as.numeric(expr_mat[gs, ])
#     } else {
#       Matrix::colMeans(expr_mat[gs, , drop = FALSE])
#     }
#   })
#   score_mat <- as.data.frame(score_mat, check.names = FALSE)
#   score_mat$cell <- colnames(treg_focus)
#
#   md <- treg_focus@meta.data
#   md$cell <- rownames(md)
#
#   # Ensure columns exist
#   if (!(group_col %in% colnames(md))) stop("[EXT++/FUNC] meta.data missing group_col: ", group_col)
#   if (!(subclass_col %in% colnames(md))) stop("[EXT++/FUNC] meta.data missing subclass_col: ", subclass_col)
#
#   df_wide <- merge(md[, c("cell", group_col, subclass_col), drop = FALSE], score_mat, by = "cell", all.x = TRUE)
#   colnames(df_wide)[colnames(df_wide) == group_col] <- "group"
#   colnames(df_wide)[colnames(df_wide) == subclass_col] <- "subclass"
#
#   df_wide$group <- factor(as.character(df_wide$group), levels = group_levels)
#   df_wide$subclass <- factor(as.character(df_wide$subclass), levels = subclass_levels)
#
#   # Keep only wanted groups & subclasses
#   df_wide <- df_wide[!is.na(df_wide$group) & !is.na(df_wide$subclass), , drop = FALSE]
#   df_wide <- df_wide[df_wide$group %in% group_levels & df_wide$subclass %in% subclass_levels, , drop = FALSE]
#
#   utils::write.csv(df_wide, file.path(out_func, "coreTreg_function_scores_wide.csv"), row.names = FALSE)
#
#   # --- Long format using utils::stack (avoid base::stack error) ---
#   score_cols <- names(func_gene_sets_present)
#   st <- utils::stack(df_wide[, score_cols, drop = FALSE])  # values, ind
#
#   df_long <- data.frame(
#     cell = rep(df_wide$cell, times = length(score_cols)),
#     group = rep(df_wide$group, times = length(score_cols)),
#     subclass = rep(df_wide$subclass, times = length(score_cols)),
#     signature = as.character(st$ind),
#     score = as.numeric(st$values),
#     stringsAsFactors = FALSE
#   )
#   df_long$group <- factor(df_long$group, levels = group_levels)
#   df_long$subclass <- factor(df_long$subclass, levels = subclass_levels)
#   df_long$signature <- factor(df_long$signature, levels = score_cols)
#
#   utils::write.csv(df_long, file.path(out_func, "coreTreg_function_scores_long.csv"), row.names = FALSE)
#
#   # --- (A) Violin+box by group, facet by subclass; one PDF per signature ---
#   for (sig in score_cols) {
#     d <- df_long[df_long$signature == sig, , drop = FALSE]
#     if (nrow(d) == 0) next
#
#     p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = score)) +
#       ggplot2::geom_violin(trim = TRUE, scale = "width") +
#       ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.2) +
#       ggplot2::facet_wrap(~subclass, nrow = 1, scales = "free_y") +
#       ggplot2::theme_bw() +
#       ggplot2::labs(
#         title = paste0("coreTreg functional score | ", sig, " (by group; stratified by subclass)"),
#         x = "Group", y = "Mean log-normalized expression (gene-set score)"
#       )
#
#     ggplot2::ggsave(
#       filename = file.path(out_func, paste0("coreTreg_", sig, "_violin_bySubclass.pdf")),
#       plot = p, width = 13, height = 4.8
#     )
#   }
#
#   # --- (B) Wilcoxon tests at cell level within each subclass (descriptive; unbalanced n warning applies) ---
#   pair_list <- list(
#     c("control1", "s1_hIL"),
#     c("control1", "s40kD"),
#     c("s1_hIL", "s40kD")
#   )
#
#   safe_wilcox <- function(x, g, g1, g2) {
#     x1 <- x[g == g1]; x2 <- x[g == g2]
#     if (length(x1) < 5 || length(x2) < 5) return(c(p = NA_real_, W = NA_real_))
#     tt <- stats::wilcox.test(x1, x2, exact = FALSE)
#     c(p = as.numeric(tt$p.value), W = as.numeric(tt$statistic))
#   }
#
#   res_list <- list()
#   idx <- 1L
#   for (sb in subclass_levels) {
#     for (sig in score_cols) {
#       d <- df_long[df_long$subclass == sb & df_long$signature == sig, , drop = FALSE]
#       if (nrow(d) == 0) next
#       for (pp in pair_list) {
#         g1 <- pp[1]; g2 <- pp[2]
#         if (!(g1 %in% d$group) || !(g2 %in% d$group)) next
#         ww <- safe_wilcox(d$score, d$group, g1, g2)
#         res_list[[idx]] <- data.frame(
#           subclass = sb,
#           signature = sig,
#           contrast = paste0(g1, "_vs_", g2),
#           n_g1 = sum(d$group == g1),
#           n_g2 = sum(d$group == g2),
#           W = ww["W"],
#           p = ww["p"],
#           stringsAsFactors = FALSE
#         )
#         idx <- idx + 1L
#       }
#     }
#   }
#
#   if (length(res_list) > 0) {
#     res <- do.call(rbind, res_list)
#     res$padj <- stats::p.adjust(res$p, method = "BH")
#     utils::write.csv(res, file.path(out_func, "coreTreg_function_wilcox_bySubclass_celllevel.csv"), row.names = FALSE)
#   } else {
#     cat("[EXT++/FUNC] No Wilcoxon results produced (insufficient cells per group within subclasses).\n")
#   }
#
#   cat("[EXT++/FUNC] Outputs (PDF+CSV): ", out_func, "\n", sep = "")
# }
#
#
# # ---- eTreg-like 内部：模块分数 + 关键基因（解释 s40kD vs s1_hIL）----
# meta_f <- treg_focus@meta.data
# cells_et <- rownames(meta_f)[as.character(meta_f[[subclass_col]]) == "eTreg-like"]
# et_focus <- subset(treg_focus, cells = cells_et)
#
# cat("\n[Node-5] eTreg-like focus cells:", ncol(et_focus), "\n")
# cat("[Node-5] eTreg-like per group:\n"); print(base_table(as.character(et_focus@meta.data[[group_col]])))
#
# # 先创建基础图
# # ---------- [REPLACE] 7.2 Plot WITHOUT VlnPlot (base boxplot + jitter; robust) ----------
# # 目的：绕开 Seurat::VlnPlot 的 tidy-eval/patchwork 路径，避免 'S4SXP' 内部错误
# # 输出文件名保持不变，便于你后续流程无缝衔接
#
# plot_boxjitter_pdf <- function(df, group_col, feat_cols, out_pdf,
#                                groups_order = NULL,
#                                nrow = 2, ncol = 2, point_cex = 0.25) {
#   if (!group_col %in% colnames(df)) stop("group_col not in df: ", group_col)
#   feat_cols <- feat_cols[feat_cols %in% colnames(df)]
#   if (length(feat_cols) == 0) stop("No valid feature columns to plot.")
#
#   grp <- as.character(df[[group_col]])
#   if (!is.null(groups_order)) {
#     keep <- grp %in% groups_order
#     df <- df[keep, , drop = FALSE]
#     grp <- grp[keep]
#     grp <- factor(grp, levels = groups_order)
#   } else {
#     grp <- factor(grp)
#   }
#
#   grDevices::pdf(out_pdf, width = 14, height = 7, onefile = TRUE)
#   on.exit(grDevices::dev.off(), add = TRUE)
#
#   per_page <- nrow * ncol
#   oldpar <- graphics::par(no.readonly = TRUE)
#   on.exit(graphics::par(oldpar), add = TRUE)
#
#   graphics::par(mfrow = c(nrow, ncol), mar = c(6, 4, 3, 1), oma = c(0, 0, 2, 0))
#
#   for (i in seq_along(feat_cols)) {
#     feat <- feat_cols[i]
#     y <- df[[feat]]
#     y <- as.numeric(y)
#
#     graphics::boxplot(y ~ grp, outline = FALSE, las = 2,
#                       main = feat, ylab = "Value", xlab = "")
#     graphics::stripchart(y ~ grp, vertical = TRUE, method = "jitter",
#                          pch = 16, cex = point_cex, add = TRUE)
#
#     if (i %% per_page == 0 && i < length(feat_cols)) {
#       # 新页：重置面板
#       graphics::par(mfrow = c(nrow, ncol), mar = c(6, 4, 3, 1), oma = c(0, 0, 2, 0))
#     }
#   }
# }
# # 7.2 模块分数图：只画 AddModuleScore 的 4 个模块列（避免 score_cols 被后续覆盖）
# score_cols_ms <- c("Prolif_CC_MS1", "Apoptosis_MS2", "Migration_MS3", "IL2_axis_MS4")
#
# # 再做一次稳健过滤（防止未来某次 AddModuleScore 列号变化）
# md_cols <- colnames(et_focus@meta.data)
# score_cols_ms <- score_cols_ms[score_cols_ms %in% md_cols]
# if (length(score_cols_ms) < 1) {
#   stop("No module score columns found. Available MS-like columns: ",
#        paste(grep("_MS\\d+$", md_cols, value = TRUE), collapse = ", "))
# }
#
# df_ms <- et_focus@meta.data[, c(group_col, score_cols_ms), drop = FALSE]
#
# out_ms <- file.path(outdir_closure, "eTreg_like_ModuleScores_VlnPlot_by_group.pdf")
# plot_boxjitter_pdf(
#   df = df_ms,
#   group_col = group_col,
#   feat_cols = score_cols_ms,
#   out_pdf = out_ms,
#   groups_order = wanted_groups,
#   nrow = 2, ncol = 2
# )
# cat("[OK] Wrote:", out_ms, "\n")
#
# # 7.3 关键基因图：用 FetchData 拉表达（避免 Seurat 绘图层）
# key_genes <- c("Foxp3","Il2ra","Il2rb","Ctla4","Mki67","Top2a","Tigit","Lag3","Pdcd1")
# key_genes_present <- base_intersect(key_genes, rownames(et_focus))
#
# if (length(key_genes_present) > 0) {
#   df_key <- Seurat::FetchData(et_focus, vars = c(group_col, key_genes_present))
#   out_key <- file.path(outdir_closure, "eTreg_like_KeyGenes_VlnPlot_by_group.pdf")
#
#   # 基因数可能>4，自动多页（2x3 每页6个更合适）
#   plot_boxjitter_pdf(
#     df = df_key,
#     group_col = group_col,
#     feat_cols = key_genes_present,
#     out_pdf = out_key,
#     groups_order = wanted_groups,
#     nrow = 2, ncol = 3
#   )
#   cat("[OK] Wrote:", out_key, "\n")
# } else {
#   warning("Key genes not found in object. Skip eTreg key-gene plot.")
# }
# # ---------- [Node-6] Write interpretation template ----------
# txt <- file.path(outdir_closure, "closure_interpretation_template.txt")
#
# cat(
#   "Minimal closure readout (to write Supplement Results/Discussion):\n",
#   "\nA) coreTreg fraction within treg_main:\n",
#   "- coreTreg_fraction_within_treg_main.csv\n",
#   "\nA2) coreTreg fraction within ALL CELLS (global denominator):\n",
#   "- coreTreg_fraction_within_allcells.csv (if combined_* available)\n",
#   "\nB) Within-coreTreg redistribution:\n",
#   "- treg_core_subclass_counts.csv / treg_core_subclass_proportions.csv\n",
#   "\nC) Mechanism decomposition (3 classes + IL-2 axis):\n",
#   "- eTreg_like_ModuleScores_VlnPlot_by_group.pdf\n",
#   "- eTreg_like_KeyGenes_VlnPlot_by_group.pdf\n",
#   "\nInterpretation logic:\n",
#   "- If A2 shows coreTreg/allcells lower in s40kD, it supports true abundance/composition reduction.\n",
#   "- If B shows eTreg-like fraction within coreTreg drops in s40kD, it supports internal subset redistribution.\n",
#   "- If Prolif_CC lower and/or Apoptosis higher and/or Migration lower in s40kD, classify '比例减少' into 3 mechanisms.\n",
#   "- Compare IL2_axis + Il2ra/Il2rb between s40kD and s1_hIL to explain divergent IL-2 signaling engagement.\n",
#   file = txt
# )
#
# cat("\n[Node-6] Wrote interpretation template:", txt, "\n")
# cat(
#   "Minimal closure readout (to write Supplement Results/Discussion):\n",
#   "\nA) coreTreg fraction within treg_main:\n",
#   "- coreTreg_fraction_within_treg_main.csv\n",
#   "\nA2) coreTreg fraction within ALL CELLS (global denominator):\n",
#   "- coreTreg_fraction_within_allcells.csv (if combined_* available)\n",
#   "\nB) Within-coreTreg redistribution:\n",
#   "- treg_core_subclass_counts.csv / treg_core_subclass_proportions.csv\n",
#   "\nC) Mechanism decomposition (3 classes + IL-2 axis):\n",
#   "- eTreg_like_ModuleScores_VlnPlot_by_group.pdf\n",
#   "- eTreg_like_KeyGenes_VlnPlot_by_group.pdf\n",
#   "\nInterpretation logic:\n",
#   "- If A2 shows coreTreg/allcells lower in s40kD, it supports true abundance/composition reduction.\n",
#   "- If B shows eTreg-like fraction within coreTreg drops in s40kD, it supports internal subset redistribution.\n",
#   "- If Prolif_CC lower and/or Apoptosis higher and/or Migration lower in s40kD, classify '比例减少' into 3 mechanisms.\n",
#   "- Compare IL2_axis + Il2ra/Il2rb between s40kD and s1_hIL to explain divergent IL-2 signaling engagement.\n",
#   file = txt
# )
# cat("\n[Node-6] Wrote interpretation template:", txt, "\n")
#
# cat("\n==== [MOD] Step 1.7-CLOSURE-EXT DONE. Outputs in: ", outdir_closure, " ====\n", sep = "")
#
# # #============================================================
# # [ADD] Step 1.7-CLOSURE-EXT+ (TIGHTENING PACK; minimal extra work)
# #   目的：在不增加大工程的前提下，夯实“比例减少”与“效应差异”的证据链
# #   输出：
# #     1) eTreg-like cell-cycle phase composition by group (G1/S/G2M)
# #     2) Bootstrap CI for coreTreg subclass proportions by group
# #     3) Optional pre-ranked GSEA (if fgsea+msigdbr installed)
# # ============================================================
#
# cat("\n==== [ADD] Step 1.7-CLOSURE-EXT+ (TIGHTENING PACK) START ====\n")
#
# # ---------- [ADD-FIX] base:: 显式引用，避免 conflicted ----------
# base_intersect <- base::intersect
# base_unname    <- base::unname
# base_table     <- base::table
# base_names     <- base::names
#
# # ---------- [ADD] 输出目录 ----------
# outdir_tight <- "./temp/Step1.7_CLOSURE_EXT_TIGHT"
# dir.create(outdir_tight, showWarnings = FALSE, recursive = TRUE)
#
# # ---------- [ADD] 依赖你在 CLOSURE-EXT 里已经创建的对象 ----------
# # 需要这些对象在内存中存在；若你是分段运行脚本，可自动重载
# if (!exists("treg_main")) {
#   treg_rds <- "./temp/treg_main_labeled_step1.2.rds"
#   if (!file.exists(treg_rds)) stop("找不到：", treg_rds, "。请先完成 Step1.2 并保存该 RDS。")
#   treg_main <- readRDS(treg_rds)
#   DefaultAssay(treg_main) <- "RNA"
# }
# if (!exists("treg_core")) {
#   # 与你 CLOSURE-EXT 一致的 core gate（显式 base::intersect）
#   need_core <- base_intersect(c("Foxp3","Il2ra","Trac","Cd3e"), rownames(treg_main))
#   if (length(need_core) < 3) stop("关键 Treg core 基因缺失过多（Foxp3/Il2ra/Trac/Cd3e）。请检查基因命名体系。")
#   treg_core <- subset(treg_main, subset = Foxp3 > 0 & Il2ra > 0 & Trac > 0 & Cd3e > 0)
# }
# if (!exists("group_col"))    group_col <- "group"
# if (!exists("subclass_col")) subclass_col <- "Treg_subclass_cluster"
# wanted_groups <- c("control1","s1_hIL","s40kD")
#
# # ============================================================
# # (1) eTreg-like 内：Cell-cycle phase proportions (G1/S/G2M)
# # ============================================================
# cat("\n[ADD-1] eTreg-like cell-cycle phase composition by group ...\n")
#
# # 取 eTreg-like（在 coreTreg 内）
# meta_core <- treg_core@meta.data
# cells_et <- rownames(meta_core)[as.character(meta_core[[subclass_col]]) == "eTreg-like" &
#                                   as.character(meta_core[[group_col]]) %in% wanted_groups]
# et_focus_cc <- subset(treg_core, cells = cells_et)
#
# cat("[ADD-1] eTreg-like cells:", ncol(et_focus_cc), "\n")
# cat("[ADD-1] eTreg-like per group:\n")
# print(base_table(as.character(et_focus_cc@meta.data[[group_col]])))
#
# # CellCycleScoring 需要 S/G2M 基因表：使用 Seurat 内置 cc.genes.updated.2019（若存在）
# # 并用 CaseMatch 适配你的基因命名（鼠基因常为首字母大写：Mcm5 等）
# cc_ok <- TRUE
# if (exists("cc.genes.updated.2019")) {
#   cc <- cc.genes.updated.2019
# } else if ("Seurat" %in% loadedNamespaces() && exists("cc.genes.updated.2019", where = asNamespace("Seurat"))) {
#   cc <- get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
# } else {
#   cc_ok <- FALSE
#   warning("Seurat 内置 cc.genes.updated.2019 未找到。跳过 CellCycleScoring；你仍可用 Prolif module score 作为替代。")
# }
#
# if (cc_ok) {
#   s.genes   <- Seurat::CaseMatch(cc$s.genes, rownames(et_focus_cc))
#   g2m.genes <- Seurat::CaseMatch(cc$g2m.genes, rownames(et_focus_cc))
#
#   if (length(s.genes) < 10 || length(g2m.genes) < 10) {
#     warning("Cell-cycle gene match 太少（S 或 G2M < 10）。可能是命名体系不匹配。跳过 CellCycleScoring。")
#   } else {
#     # 不重跑大流程：CellCycleScoring 用 normalize data；若不存在，轻量 NormalizeData
#     if (!"RNA" %in% names(et_focus_cc@assays)) stop("对象缺少 RNA assay。")
#     # 如果对象 data slot 为空，做一次轻量 NormalizeData
#     if (ncol(Seurat::GetAssayData(et_focus_cc, assay = "RNA", slot = "data")) == 0) {
#       et_focus_cc <- NormalizeData(et_focus_cc, verbose = FALSE)
#     }
#
#     et_focus_cc <- Seurat::CellCycleScoring(
#       et_focus_cc,
#       s.features   = s.genes,
#       g2m.features = g2m.genes,
#       set.ident    = FALSE
#     )
#
#     # 统计 Phase 占比（按组）
#     df_phase <- et_focus_cc@meta.data[, c(group_col, "Phase"), drop = FALSE]
#     df_phase <- df_phase[as.character(df_phase[[group_col]]) %in% wanted_groups, , drop = FALSE]
#
#     tab_phase <- base_table(as.character(df_phase$Phase), as.character(df_phase[[group_col]]))
#     prop_phase <- sweep(tab_phase, 2, colSums(tab_phase), FUN = "/")
#
#     # 输出表  [FIX] Export tables (write.csv is in utils, not base) ----------
# utils::write.csv(
#   as.data.frame.matrix(tab_phase),
#   file = file.path(outdir_tight, "eTreg_like_cellcycle_phase_counts_by_group.csv"),
#   row.names = TRUE
# )
#
# utils::write.csv(
#   as.data.frame.matrix(prop_phase),
#   file = file.path(outdir_tight, "eTreg_like_cellcycle_phase_prop_by_group.csv"),
#   row.names = TRUE
# )
#
#     # 画图（堆叠比例）
#     dfp <- as.data.frame(as.table(prop_phase), stringsAsFactors = FALSE)
#     colnames(dfp) <- c("Phase","group","prop")
#
#     p_phase <- ggplot(dfp, aes(x = group, y = prop, fill = Phase)) +
#       geom_col() +
#       theme_bw() +
#       labs(title = "eTreg-like: Cell-cycle phase composition by group (descriptive)",
#            x = "Group", y = "Proportion within eTreg-like")
#
#     ggsave(file.path(outdir_tight, "eTreg_like_cellcycle_phase_prop_by_group.png"),
#            p_phase, width = 10, height = 5)
#
#     cat("[ADD-1] Wrote cell-cycle outputs under:", outdir_tight, "\n")
#   }
# }
#
# # ============================================================
# # (2) coreTreg 内：subclass proportions bootstrap 95% CI
# # ============================================================
# cat("\n[ADD-2] Bootstrap CI for coreTreg subclass proportions ...\n")
#
# meta2 <- treg_core@meta.data
# meta2$grp <- as.character(meta2[[group_col]])
# meta2$sub <- as.character(meta2[[subclass_col]])
# meta2 <- meta2[meta2$grp %in% wanted_groups & !is.na(meta2$sub), , drop = FALSE]
#
# sub_levels <- sort(unique(meta2$sub))
# grp_levels <- wanted_groups[wanted_groups %in% sort(unique(meta2$grp))]
#
# # bootstrap 函数：组内重采样，计算各 sub 的比例
# bootstrap_props_one_group <- function(subvec, B = 500L) {
#   n <- length(subvec)
#   out <- matrix(NA_real_, nrow = B, ncol = length(sub_levels))
#   colnames(out) <- sub_levels
#   for (b in seq_len(B)) {
#     s <- sample(subvec, size = n, replace = TRUE)
#     tab <- base::table(factor(s, levels = sub_levels))
#     out[b, ] <- as.numeric(tab) / sum(tab)
#   }
#   out
# }
#
# B <- 500L  # 500 次足够“最短闭环”；想更稳可改 1000
# boot_list <- list()
# for (g in grp_levels) {
#   subvec <- meta2$sub[meta2$grp == g]
#   boot_list[[g]] <- bootstrap_props_one_group(subvec, B = B)
# }
#
# # 汇总：mean + 95% CI
# df_ci <- do.call(rbind, lapply(grp_levels, function(g) {
#   mat <- boot_list[[g]]
#   data.frame(
#     group = g,
#     subclass = colnames(mat),
#     mean = apply(mat, 2, mean, na.rm = TRUE),
#     lo95 = apply(mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
#     hi95 = apply(mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
#     stringsAsFactors = FALSE
#   )
# }))
#
# utils::write.csv(df_ci, file = file.path(outdir_tight, "coreTreg_subclass_prop_bootstrap95CI.csv"),
#                 row.names = FALSE)
#
# # 作图：误差条（每个 subclass 一个 facet，审稿友好）
# p_ci <- ggplot(df_ci, aes(x = group, y = mean)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.15) +
#   facet_wrap(~subclass, scales = "free_y") +
#   theme_bw() +
#   labs(title = paste0("coreTreg subclass proportion by group (bootstrap 95% CI; B=", B, ")"),
#        x = "Group", y = "Proportion within coreTreg")
#
# ggsave(file.path(outdir_tight, "coreTreg_subclass_prop_bootstrap95CI.png"),
#        p_ci, width = 12, height = 7)
#
# cat("[ADD-2] Wrote bootstrap CI outputs under:", outdir_tight, "\n")
#
# # ============================================================
# # (3) Optional: pre-ranked GSEA for eTreg-like pseudo-bulk ranks
#
# # ---------- [ADD-3 MIN] Always export standardized prerank lists ----------
# rank_dir <- "./temp/Step1.7_eTreg_publishable"
# out_prerank_dir <- outdir_tight
# dir.create(out_prerank_dir, showWarnings = FALSE, recursive = TRUE)
#
# targets <- c(
#   "control1_vs_s1_hIL",
#   "control1_vs_s40kD",
#   "s1_hIL_vs_s40kD"
# )
#
# for (nm in targets) {
#   in_file <- file.path(rank_dir, paste0("logFC_rank_", nm, ".csv"))
#   if (!file.exists(in_file)) {
#     cat("[ADD-3 MIN][WARN] missing:", in_file, "\n")
#     next
#   }
#   df <- utils::read.csv(in_file, stringsAsFactors = FALSE)
#   if (!all(c("gene","logFC") %in% colnames(df))) {
#     stop("Rank file missing gene/logFC columns: ", in_file)
#   }
#   df <- df[!is.na(df$gene) & !is.na(df$logFC), c("gene","logFC"), drop = FALSE]
#   df <- df[order(df$logFC, decreasing = TRUE), , drop = FALSE]
#
#   out_file <- file.path(out_prerank_dir, paste0("GSEA_prerank_logFC_rank_", nm, ".csv"))
#   utils::write.csv(df, out_file, row.names = FALSE)
#   cat("[ADD-3 MIN] wrote:", out_file, "\n")
# }
# # ============================================================
# cat("\n[ADD-3] Optional pre-ranked GSEA (skip gracefully if packages missing) ...\n")
#
# # 你 Step1.7_eTreg_publishable 已输出 logFC rank CSV：优先用这些作为 pre-ranked 输入
# rank_dir <- "./temp/Step1.7_eTreg_publishable"
# rank_files <- c(
#   "logFC_rank_control1_vs_s40kD.csv",
#   "logFC_rank_s1_hIL_vs_s40kD.csv",
#   "logFC_rank_control1_vs_s1_hIL.csv"
# )
# rank_files_full <- file.path(rank_dir, rank_files)
# rank_files_full <- rank_files_full[file.exists(rank_files_full)]
#
# if (length(rank_files_full) == 0) {
#   cat("[ADD-3] No rank files found under:", rank_dir, "\n",
#       "Skip GSEA. (You can rerun Step1.7_eTreg_publishable first.)\n")
# } else {
#   # 尝试使用 fgsea + msigdbr（若安装了）
#   has_fgsea  <- requireNamespace("fgsea", quietly = TRUE)
#   has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
#
#   if (!has_fgsea || !has_msigdbr) {
#     cat("[ADD-3] fgsea/msigdbr not installed. Will export ranked lists for later GSEA.\n")
#     # 仅导出最小 ranked list（gene -> stat），供你后续在服务器/本机装包后直接跑
#     for (f in rank_files_full) {
#       df <- utils::read.csv(f, stringsAsFactors = FALSE)
#       if (!all(c("gene","logFC") %in% colnames(df))) next
#       r <- df[, c("gene","logFC")]
#       out <- file.path(outdir_tight, paste0("GSEA_prerank_", basename(f)))
#       utils::write.csv(r, out, row.names = FALSE)
#       cat("[ADD-3] Wrote pre-rank list:", out, "\n")
#     }
#   } else {
#     cat("[ADD-3] fgsea + msigdbr available. Running pre-ranked GSEA (Hallmark + C2/CP) ...\n")
#
#     # 使用 msigdbr 的 Hallmark (H) 作为最稳的最小集合；物种默认 human
#     # 若你是 mouse symbol，后续可将 species = "Mus musculus"（若 msigdbr 支持你的版本）
#     msig_h <- msigdbr::msigdbr(category = "H")
#     pathways <- split(msig_h$gene_symbol, msig_h$gs_name)
#
#     for (f in rank_files_full) {
#       df <- utils::read.csv(f, stringsAsFactors = FALSE)
#       if (!all(c("gene","logFC") %in% colnames(df))) next
#
#       stats <- df$logFC
#       names(stats) <- df$gene
#       stats <- sort(stats, decreasing = TRUE)
#
#       res <- fgsea::fgsea(pathways = pathways, stats = stats, nperm = 5000)
#       res <- res[order(res$pval, -abs(res$NES)), ]
#
#       out <- file.path(outdir_tight, paste0("fgsea_Hallmark_", tools::file_path_sans_ext(basename(f)), ".csv"))
#       utils::write.csv(res, out, row.names = FALSE)
#       cat("[ADD-3] Wrote:", out, "\n")
#     }
#   }
# }
#
# cat("\n==== [ADD] Step 1.7-CLOSURE-EXT+ DONE. Outputs in: ", outdir_tight, " ====\n", sep = "")
#
# # ============================================================
# # [ADD] Step 1.7-CLOSURE-EXT++ : GSEA prerank via fgsea (publishable)
# #   Input: ./temp/Step1.7_CLOSURE_EXT_TIGHT/GSEA_prerank_logFC_rank_*.csv (gene, logFC)
# #   Output: ./temp/Step1.7_CLOSURE_EXT/GSEA_fGSEA/
# # ============================================================
#
# cat("\n==== [ADD] Step 1.7-CLOSURE-EXT++ (fgsea prerank; publishable) START ====\n")
# # ============================================================
# # ============================================================
# # [ADD] Step 1.7-CLOSURE-EXT++ : GSEA prerank via fgsea (publishable)
# #   Input:  ./temp/Step1.7_CLOSURE_EXT_TIGHT/GSEA_prerank_logFC_rank_*.csv  (gene, logFC)
# #           或 ./temp/Step1.7_CLOSURE_EXT_TIGHT/GSEA_prerank_logFC_rank_*.csv
# #   Output: ./temp/Step1.7_CLOSURE_EXT/GSEA_fGSEA/fgsea_Hallmark_*.csv
# # ============================================================
#
# cat("\n==== [ADD] Step 1.7-CLOSURE-EXT++ (fgsea prerank; publishable) START ====\n")
#
# outdir_fgsea <- "./temp/Step1.7_CLOSURE_EXT/GSEA_fGSEA"
# dir.create(outdir_fgsea, showWarnings = FALSE, recursive = TRUE)
#
# has_fgsea   <- requireNamespace("fgsea", quietly = TRUE)
# has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
#
# if (!has_fgsea || !has_msigdbr) {
#   cat("[EXT++] fgsea/msigdbr not installed. Skip fgsea run.\n")
#   cat("Outputs not generated under: ", outdir_fgsea, "\n", sep = "")
# } else {
#
#   # ---- 1) 读取 prerank 文件（你现在已有的三个）----
#   prerank_dir <- "./temp/Step1.7_CLOSURE_EXT_TIGHT"
#   rank_files <- list.files(prerank_dir, pattern = "^GSEA_prerank_logFC_rank_.*\\.csv$", full.names = TRUE)
#
#   # 如果你把 prerank 放在 Step1.7_eTreg_publishable，也兼容
#   if (length(rank_files) == 0) {
#     prerank_dir2 <- "./temp/Step1.7_eTreg_publishable"
#     rank_files2 <- list.files(prerank_dir2, pattern = "^logFC_rank_.*\\.csv$", full.names = TRUE)
#     # 把 logFC_rank_* 转成 prerank（gene, logFC）
#     if (length(rank_files2) > 0) {
#       for (f in rank_files2) {
#         df <- utils::read.csv(f, stringsAsFactors = FALSE)
#         if (!all(c("gene","logFC") %in% colnames(df))) next
#         out <- file.path(prerank_dir2, paste0("GSEA_prerank_", basename(f)))
#         utils::write.csv(df[, c("gene","logFC")], out, row.names = FALSE)
#       }
#       rank_files <- list.files(prerank_dir2, pattern = "^GSEA_prerank_logFC_rank_.*\\.csv$", full.names = TRUE)
#     }
#   }
#
#   if (length(rank_files) == 0) {
#     stop("[EXT++] No prerank CSV found. Expect: GSEA_prerank_logFC_rank_*.csv")
#   } else {
#     cat("[EXT++] Found prerank files:\n")
#     print(rank_files)
#   }
#
#   # ---- 2) 关键修正：用小鼠 Hallmark（否则极易 overlap=0 导致空结果）----
#   msig_h <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
#   pathways <- split(msig_h$gene_symbol, msig_h$gs_name)
#
#   # ---- 3) 循环跑 fgsea（带 overlap 诊断 + leadingEdge 写出修复）----
#   run_one <- function(f) {
#     df <- utils::read.csv(f, stringsAsFactors = FALSE)
#     if (!all(c("gene","logFC") %in% colnames(df))) {
#       cat("[SKIP] Missing gene/logFC columns in: ", f, "\n", sep = "")
#       return(invisible(NULL))
#     }
#
#     # stats：去 NA、去空 gene、合并重复 gene（取绝对值更大的那个 logFC）
#     df <- df[is.finite(df$logFC) & nzchar(df$gene), , drop = FALSE]
#     df <- df[order(abs(df$logFC), decreasing = TRUE), , drop = FALSE]
#     df <- df[!duplicated(df$gene), , drop = FALSE]
#
#     stats <- df$logFC
#     names(stats) <- df$gene
#     stats <- sort(stats, decreasing = TRUE)
#
#     # overlap 诊断：如果这里很小，fgsea 很可能空
#     ov <- vapply(pathways, function(p) sum(p %in% names(stats)), FUN.VALUE = integer(1))
#     cat("\n[EXT++] ", basename(f), " | stats genes=", length(stats),
#         " | pathway overlap summary: ",
#         "median=", stats::median(ov), ", max=", max(ov), "\n", sep = "")
#     # 如 overlap 太低，就把 minSize 下调
#     minSize_use <- if (max(ov) < 15) 5 else 10
#
#     res <- fgsea::fgsea(
#       pathways = pathways,
#       stats    = stats,
#       nperm    = 5000,
#       minSize  = minSize_use,
#       maxSize  = 500
#     )
#     res <- as.data.frame(res)
#     if (nrow(res) == 0) {
#       cat("[EXT++] fgsea returned 0 rows for: ", basename(f),
#           " (likely overlap too low even after minSize adjustment)\n", sep = "")
#     } else {
#       # 解决你之前的报错：leadingEdge 是 list 列，write.csv 会炸
#       if ("leadingEdge" %in% colnames(res)) {
#         res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), FUN.VALUE = character(1))
#       }
#       res <- res[order(res$pval, -abs(res$NES)), , drop = FALSE]
#     }
#
#     tag <- sub("^GSEA_prerank_logFC_rank_", "", tools::file_path_sans_ext(basename(f)))
#     out_all <- file.path(outdir_fgsea, paste0("fgsea_Hallmark_", tag, ".csv"))
#     out_top <- file.path(outdir_fgsea, paste0("fgsea_Hallmark_", tag, "_TOP30.csv"))
#
#     utils::write.csv(res, out_all, row.names = FALSE)
#     utils::write.csv(utils::head(res, 30), out_top, row.names = FALSE)
#
#     cat("[OK] Wrote: ", out_all, "\n", sep = "")
#     invisible(res)
#   }
#
#   for (f in rank_files) run_one(f)
# }
#
# # =========================
# # [VIS] fgsea Hallmark visualization (PDF only)
# # =========================
# suppressPackageStartupMessages({
#   library(ggplot2)
# })
#
# base_setdiff <- base::setdiff
# base_order   <- base::order
#
# fg_dir <- "./temp/Step1.7_CLOSURE_EXT/GSEA_fGSEA"
# out_vis <- file.path(fg_dir, "VIS_PDF")
# dir.create(out_vis, showWarnings = FALSE, recursive = TRUE)
#
# # 自动同时兼容你“带编号”的文件名：1_fgsea_... / fgsea_...
# fg_files <- list.files(fg_dir, pattern = "fgsea_Hallmark_.*\\.csv$", full.names = TRUE)
#
# safe_read_fg <- function(f) {
#   df <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
#   if (is.null(df) || nrow(df) == 0) return(NULL)
#   need <- c("pathway","NES","pval","padj","size")
#   miss <- base_setdiff(need, colnames(df))
#   if (length(miss) > 0) return(NULL)
#
#   df$pathway2 <- gsub("^HALLMARK_", "", df$pathway)
#   df$pathway2 <- gsub("_", " ", df$pathway2)
#   df$padj2 <- pmax(df$padj, 1e-300)
#   df$mlog10padj <- -log10(df$padj2)
#   df
# }
#
# dfs <- lapply(fg_files, safe_read_fg)
# dfs <- Filter(Negate(is.null), dfs)
#
# if (length(dfs) == 0) {
#   stop("fgsea csv 全是空的或列不全。请确认 fgsea 输出目录与文件内容。")
# }
#
# for (i in seq_along(dfs)) {
#   df <- dfs[[i]]
#   bn <- tools::file_path_sans_ext(basename(fg_files[[i]]))
#   # bn 形如 fgsea_Hallmark_control1_vs_s40kD 或 1_fgsea_Hallmark_...
#   bn <- sub("^\\d+_", "", bn)
#
#   # 取 padj 最小的 top20 作 dotplot
#   df1 <- df[base_order(df$padj, -abs(df$NES)), , drop = FALSE]
#   df_top <- utils::head(df1, 20)
#
#   # dotplot
#   df_top$pathway2 <- factor(df_top$pathway2, levels = rev(df_top$pathway2))
#   p_dot <- ggplot(df_top, aes(x = NES, y = pathway2, size = mlog10padj)) +
#     geom_point() +
#     theme_bw() +
#     labs(title = paste0(bn, " | Hallmark (fgsea)"),
#          x = "NES", y = NULL, size = "-log10(padj)")
#
#   ggsave(file.path(out_vis, paste0(bn, "_dotplot.pdf")), p_dot, width = 10, height = 6)
#
#   # barplot: top 10 up / top 10 down
#   df_up <- df[base_order(-df$NES), , drop = FALSE]; df_up <- utils::head(df_up, 10)
#   df_dn <- df[base_order(df$NES),  , drop = FALSE]; df_dn <- utils::head(df_dn, 10)
#   df_bar <- rbind(df_up, df_dn)
#   df_bar$pathway2 <- gsub("^HALLMARK_", "", df_bar$pathway)
#   df_bar$pathway2 <- gsub("_", " ", df_bar$pathway2)
#   df_bar$pathway2 <- factor(df_bar$pathway2, levels = df_bar$pathway2[base_order(df_bar$NES)])
#
#   p_bar <- ggplot(df_bar, aes(x = pathway2, y = NES)) +
#     geom_col() +
#     coord_flip() +
#     theme_bw() +
#     labs(title = paste0(bn, " | Top NES pathways"),
#          x = NULL, y = "NES")
#
#   ggsave(file.path(out_vis, paste0(bn, "_barplot.pdf")), p_bar, width = 10, height = 7)
# }
#
# cat("[OK] fgsea VIS PDFs saved in: ", out_vis, "\n", sep = "")
#
#
# cat("==== [ADD] Step 1.7-CLOSURE-EXT++ DONE ====\n")
#
# # ============================================================
# # [ADD] Step 1.7-CLOSURE-EXT++ : Single-cell Hallmark scoring + tests
# #   3.1 Hallmark: single-cell score + group test (recommended over pseudobulk fgsea)
# #   3.2 Key regulators: gene panel (violin/dot + simple logFC on per-cell mean)
# #   3.4 (Optional) IL-2 ligand source in ALL cells + receptor/response in eTreg-like
# #
# #   Outputs:
# #     ./temp/Step1.7_CLOSURE_EXT/SC_Hallmark/
# #     ./temp/Step1.7_CLOSURE_EXT/KeyGenes/
# #     ./temp/Step1.7_CLOSURE_EXT/IL2_LigandSource/
# # ============================================================
#
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(Matrix)
#   library(ggplot2)
# })
#
# # ---- explicit base/utils to avoid conflicted ----
# b_intersect <- base::intersect
# b_setdiff   <- base::setdiff
# b_match     <- base::match
# b_table     <- base::table
# b_names     <- base::names
# b_any       <- base::any
# b_all       <- base::all
#
# # ---- output roots (reuse your CLOSURE_EXT) ----
# out_root <- "./temp/Step1.7_CLOSURE_EXT"
# dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
#
# out_sc_hm  <- file.path(out_root, "SC_Hallmark")
# out_key    <- file.path(out_root, "KeyGenes")
# out_il2src <- file.path(out_root, "IL2_LigandSource")
# dir.create(out_sc_hm,  showWarnings = FALSE, recursive = TRUE)
# dir.create(out_key,    showWarnings = FALSE, recursive = TRUE)
# dir.create(out_il2src, showWarnings = FALSE, recursive = TRUE)
#
# # ---- reload minimal objects if not in memory ----
# treg_rds <- "./temp/treg_main_labeled_step1.2.rds"
# if (!exists("treg_main")) {
#   if (!file.exists(treg_rds)) stop("Missing RDS: ", treg_rds)
#   treg_main <- readRDS(treg_rds)
#   DefaultAssay(treg_main) <- "RNA"
# }
# group_col    <- if (exists("group_col")) group_col else "group"
# subclass_col <- if (exists("subclass_col")) subclass_col else "Treg_subclass_cluster"
# wanted_groups <- c("control1","s1_hIL","s40kD")
#
# if (!(group_col %in% colnames(treg_main@meta.data))) {
#   stop("treg_main@meta.data missing group_col = ", group_col)
# }
# if (!(subclass_col %in% colnames(treg_main@meta.data))) {
#   stop("treg_main@meta.data missing subclass_col = ", subclass_col)
# }
#
# # ---- core gate (consistent with your earlier steps) ----
# need_core <- b_intersect(c("Foxp3","Il2ra","Trac","Cd3e"), rownames(treg_main))
# if (length(need_core) < 3) stop("Core Treg genes missing too much (Foxp3/Il2ra/Trac/Cd3e).")
# treg_core <- subset(treg_main, subset = Foxp3 > 0 & Il2ra > 0 & Trac > 0 & Cd3e > 0)
#
# # ---- eTreg-like focus object ----
# meta_core <- treg_core@meta.data
# cells_et <- rownames(meta_core)[
#   as.character(meta_core[[subclass_col]]) == "eTreg-like" &
#     as.character(meta_core[[group_col]]) %in% wanted_groups
# ]
# if (length(cells_et) == 0) stop("No eTreg-like cells found in treg_core. Check label string.")
# et_focus <- subset(treg_core, cells = cells_et)
#
# cat("\n[EXT++/SC] eTreg-like cells:", ncol(et_focus), "\n")
# cat("[EXT++/SC] eTreg-like per group:\n")
# print(b_table(as.character(et_focus@meta.data[[group_col]])))
#
# # ============================================================
# # 3.1 Hallmark: single-cell scoring + group tests
# #   - use msigdbr Hallmark (mouse) -> Seurat::AddModuleScore (no AUCell/GSVA needed)
# #   - avoid Seurat::VlnPlot; use ggplot2 directly
# # ============================================================
#
# has_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
# if (!has_msigdbr) {
#   cat("\n[EXT++/SC] msigdbr NOT installed -> skip Hallmark single-cell scoring.\n")
#   cat("  Install: install.packages('msigdbr')\n")
# } else {
#   # --- hallmark sets (mouse) ---
#   hm_all <- msigdbr::msigdbr(species = "Mus musculus", category = "H")
#   pathways_all <- split(hm_all$gene_symbol, hm_all$gs_name)
#
#   # focus sets (edit if you want more)
#   hm_focus <- c(
#     "HALLMARK_IL2_STAT5_SIGNALING",
#     "HALLMARK_E2F_TARGETS",
#     "HALLMARK_G2M_CHECKPOINT",
#     "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#     "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#     "HALLMARK_APOPTOSIS",
#     "HALLMARK_INFLAMMATORY_RESPONSE"
#   )
#   hm_focus <- hm_focus[hm_focus %in% names(pathways_all)]
#   if (length(hm_focus) == 0) stop("No hallmark focus sets found in msigdbr output (mouse).")
#
#   pathways_focus <- pathways_all[hm_focus]
#   # intersect with object genes; keep sets with >=10 matched genes (adjustable)
#   pathways_focus2 <- lapply(pathways_focus, function(gs) b_intersect(gs, rownames(et_focus)))
#   lens <- vapply(pathways_focus2, length, integer(1))
#   pathways_focus2 <- pathways_focus2[lens >= 10]
#
#   if (length(pathways_focus2) == 0) {
#     cat("\n[EXT++/SC] Hallmark matched genes too few in et_focus. Check gene symbol system.\n")
#   } else {
#     cat("\n[EXT++/SC] Hallmark sets kept (>=10 matched genes):\n")
#     print(data.frame(pathway = names(pathways_focus2), n_genes = vapply(pathways_focus2, length, integer(1))))
#
#     # AddModuleScore
#     # NOTE: AddModuleScore will append columns like HM_MS1, HM_MS2...
#     et_sc <- Seurat::AddModuleScore(
#       object   = et_focus,
#       features = base::unname(pathways_focus2),
#       name     = "HM_MS",
#       assay    = "RNA",
#       search   = FALSE
#     )
#
#     # map module score cols to pathways (order preserved)
#     ms_cols <- grep("^HM_MS\\d+$", colnames(et_sc@meta.data), value = TRUE)
#     if (length(ms_cols) != length(pathways_focus2)) {
#       # fallback: take the last N matches
#       ms_cols <- tail(ms_cols, length(pathways_focus2))
#     }
#     map_ms <- data.frame(
#       pathway = names(pathways_focus2),
#       score_col = ms_cols,
#       stringsAsFactors = FALSE
#     )
#     utils::write.csv(map_ms, file.path(out_sc_hm, "Hallmark_scorecol_map.csv"), row.names = FALSE)
#
#     # fetch scores + group to plain df (avoid Seurat plots)
#     df_sc <- Seurat::FetchData(et_sc, vars = c(group_col, ms_cols))
#     df_sc$group <- as.character(df_sc[[group_col]])
#     df_sc <- df_sc[df_sc$group %in% wanted_groups, , drop = FALSE]
#
#     # make long table (no tidyr): utils::stack
#     df_long <- utils::stack(df_sc[, ms_cols, drop = FALSE])
#     df_long$group <- rep(df_sc$group, times = length(ms_cols))
#     df_long$score_col <- as.character(df_long$ind)
#     df_long$score <- as.numeric(df_long$values)
#     df_long$pathway <- map_ms$pathway[match(df_long$score_col, map_ms$score_col)]
#
#     df_long$pathway2 <- gsub("^HALLMARK_", "", df_long$pathway)
#     df_long$pathway2 <- gsub("_", " ", df_long$pathway2)
#
#     utils::write.csv(df_long, file.path(out_sc_hm, "Hallmark_singlecell_scores_long.csv"), row.names = FALSE)
#
#     # ---- cell-level Wilcoxon (descriptive; n huge -> p very small; keep effect size too) ----
#     pair_list <- list(c("control1","s1_hIL"), c("control1","s40kD"), c("s1_hIL","s40kD"))
#     res_list <- list()
#     for (pw in unique(df_long$pathway)) {
#       subd <- df_long[df_long$pathway == pw, , drop = FALSE]
#       for (pp in pair_list) {
#         g1 <- pp[1]; g2 <- pp[2]
#         x1 <- subd$score[subd$group == g1]
#         x2 <- subd$score[subd$group == g2]
#         if (length(x1) < 5 || length(x2) < 5) next
#         wt <- stats::wilcox.test(x1, x2, exact = FALSE)
#         # effect: median diff
#         eff <- stats::median(x1) - stats::median(x2)
#         res_list[[length(res_list) + 1]] <- data.frame(
#           pathway = pw,
#           contrast = paste0(g1, "_vs_", g2),
#           n1 = length(x1), n2 = length(x2),
#           med1 = stats::median(x1), med2 = stats::median(x2),
#           med_diff = eff,
#           pval = wt$p.value,
#           stringsAsFactors = FALSE
#         )
#       }
#     }
#     res_cell <- do.call(rbind, res_list)
#     if (!is.null(res_cell) && nrow(res_cell) > 0) {
#       res_cell$padj <- stats::p.adjust(res_cell$pval, method = "BH")
#       utils::write.csv(res_cell, file.path(out_sc_hm, "Hallmark_wilcox_celllevel.csv"), row.names = FALSE)
#     }
#
#     # ---- optional sample-level aggregation if replicate column exists ----
#     detect_rep_col <- function(md) {
#       cand <- c("sample","Sample","orig.ident","donor","mouse","replicate","rep","library")
#       hit <- cand[cand %in% colnames(md)]
#       if (length(hit) == 0) return(NA_character_)
#       hit[1]
#     }
#     rep_col <- detect_rep_col(et_sc@meta.data)
#     if (!is.na(rep_col)) {
#       df_sc2 <- Seurat::FetchData(et_sc, vars = c(group_col, rep_col, ms_cols))
#       df_sc2$group <- as.character(df_sc2[[group_col]])
#       df_sc2$rep   <- as.character(df_sc2[[rep_col]])
#       df_sc2 <- df_sc2[df_sc2$group %in% wanted_groups, , drop = FALSE]
#
#       # per (group,rep) mean score
#       agg_list <- list()
#       for (mc in ms_cols) {
#         tmp <- aggregate(df_sc2[[mc]], by = list(group=df_sc2$group, rep=df_sc2$rep), FUN = mean)
#         colnames(tmp)[3] <- "mean_score"
#         tmp$score_col <- mc
#         tmp$pathway <- map_ms$pathway[match(mc, map_ms$score_col)]
#         agg_list[[length(agg_list)+1]] <- tmp
#       }
#       df_rep <- do.call(rbind, agg_list)
#       if (!is.null(df_rep) && nrow(df_rep) > 0) {
#         utils::write.csv(df_rep, file.path(out_sc_hm, "Hallmark_sampleMean_scores.csv"), row.names = FALSE)
#
#         # sample-level tests only if each group has >=2 reps
#         tab_rep <- table(df_rep$group, df_rep$rep)
#         nrep_by_group <- rowSums(tab_rep > 0)
#         if (all(nrep_by_group[wanted_groups] >= 2, na.rm = TRUE)) {
#           res_rep_list <- list()
#           for (pw in unique(df_rep$pathway)) {
#             subd <- df_rep[df_rep$pathway == pw, , drop = FALSE]
#             for (pp in pair_list) {
#               g1 <- pp[1]; g2 <- pp[2]
#               x1 <- subd$mean_score[subd$group == g1]
#               x2 <- subd$mean_score[subd$group == g2]
#               if (length(x1) < 2 || length(x2) < 2) next
#               wt <- stats::wilcox.test(x1, x2, exact = FALSE)
#               res_rep_list[[length(res_rep_list)+1]] <- data.frame(
#                 pathway = pw,
#                 contrast = paste0(g1, "_vs_", g2),
#                 nrep1 = length(x1), nrep2 = length(x2),
#                 med1 = stats::median(x1), med2 = stats::median(x2),
#                 med_diff = stats::median(x1) - stats::median(x2),
#                 pval = wt$p.value,
#                 stringsAsFactors = FALSE
#               )
#             }
#           }
#           res_rep <- do.call(rbind, res_rep_list)
#           res_rep$padj <- stats::p.adjust(res_rep$pval, method = "BH")
#           utils::write.csv(res_rep, file.path(out_sc_hm, "Hallmark_wilcox_samplelevel.csv"), row.names = FALSE)
#         } else {
#           cat("\n[EXT++/SC] Replicate column detected (", rep_col, "), but groups lack >=2 reps; sample-level tests skipped.\n", sep="")
#         }
#       }
#     }
#
#     # ---- Visualization (no Seurat plots): Violin + box for each pathway ----
#     for (pw in unique(df_long$pathway)) {
#       subd <- df_long[df_long$pathway == pw, , drop = FALSE]
#       subd$group <- factor(subd$group, levels = wanted_groups)
#
#       p <- ggplot(subd, aes(x = group, y = score)) +
#         geom_violin(trim = TRUE, scale = "width") +
#         geom_boxplot(width = 0.15, outlier.size = 0.2) +
#         theme_bw() +
#         labs(
#           title = paste0(gsub("^HALLMARK_", "", pw), " (eTreg-like single-cell score)"),
#           x = "Group", y = "Module score (AddModuleScore)"
#         )
#       ggsave(file.path(out_sc_hm, paste0("Hallmark_", pw, "_violin.pdf")), p, width = 9, height = 5)
#     }
#
#     # ---- Summary volcano-like (effect vs -log10 padj) using cell-level results ----
#     if (exists("res_cell") && !is.null(res_cell) && nrow(res_cell) > 0) {
#       res_cell$pathway2 <- gsub("^HALLMARK_", "", res_cell$pathway)
#       res_cell$mlog10padj <- -log10(pmax(res_cell$padj, 1e-300))
#       p2 <- ggplot(res_cell, aes(x = med_diff, y = mlog10padj, label = pathway2)) +
#         geom_point() +
#         theme_bw() +
#         facet_wrap(~contrast, scales = "free") +
#         labs(title = "Hallmark single-cell score: effect (median diff) vs significance", x = "Median(score_g1) - Median(score_g2)", y = "-log10(BH padj)")
#       ggsave(file.path(out_sc_hm, "Hallmark_effect_vs_significance.pdf"), p2, width = 12, height = 7)
#     }
#
#     cat("\n[EXT++/SC] Hallmark single-cell scoring outputs: ", out_sc_hm, "\n", sep="")
#   }
# }
# # ============================================================
# # [EXT++/IFN] IFN core ISG genes + correlation vs Prolif/IL2 (eTreg-like; PDF only)
# #   1) IFN core genes (IFN-I ISG + IFNγ genes): dotplot-like + per-gene violin
# #   2) Correlation: IFN score vs Prolif_CC / IL2_axis module scores (same eTreg-like cells)
# #   Output: ./temp/Step1.7_CLOSURE_EXT/SC_IFN/
# # ============================================================
#
# cat("\n[EXT++/IFN] IFN core genes + IFN score correlation analysis ...\n")
#
# # ---- prerequisites ----
# if (!exists("et_focus")) stop("[EXT++/IFN] et_focus not found. Make sure EXT++ already created eTreg-like focus object.")
# if (!exists("group_col")) group_col <- "group"
#
# out_ifn <- "./temp/Step1.7_CLOSURE_EXT/SC_IFN"
# dir.create(out_ifn, showWarnings = FALSE, recursive = TRUE)
#
# # ---- safe base refs (avoid conflicted) ----
# b_intersect <- base::intersect
# b_setdiff   <- base::setdiff
#
# # ---- 1) Define IFN core genes (mouse symbols as you specified) ----
# # IFN-I ISG core
# IFN_I_genes <- c("Isg15","Ifit1","Ifit3","Oas1a","Oas2","Irf7","Stat1")
# # IFNγ / response-oriented core
# IFNG_genes  <- c("Cxcl9","Cxcl10","Irf1","Stat1","H2-Ab1")
#
# # keep only genes present
# IFN_I_genes_present <- b_intersect(IFN_I_genes, rownames(et_focus))
# IFNG_genes_present  <- b_intersect(IFNG_genes,  rownames(et_focus))
#
# cat("[EXT++/IFN] Matched genes: IFN-I=", length(IFN_I_genes_present),
#     " | IFNγ=", length(IFNG_genes_present), "\n", sep="")
#
# if (length(IFN_I_genes_present) < 3 && length(IFNG_genes_present) < 3) {
#   cat("[EXT++/IFN] Too few IFN genes matched (<3). Skip IFN block.\n")
# } else {
#
#   # ---- 1A) Expression long table + per-gene violin (NO Seurat plotting) ----
#   IFN_all_genes <- unique(c(IFN_I_genes_present, IFNG_genes_present))
#   df_g <- Seurat::FetchData(et_focus, vars = c(group_col, IFN_all_genes))
#   df_g$group <- as.character(df_g[[group_col]])
#
#   # keep target groups if available
#   if (exists("wanted_groups")) {
#     df_g <- df_g[df_g$group %in% wanted_groups, , drop = FALSE]
#     df_g$group <- factor(df_g$group, levels = wanted_groups)
#   }
#
#   # long format via utils::stack
#   st <- utils::stack(df_g[, IFN_all_genes, drop = FALSE])   # values, ind
#   df_long_g <- data.frame(
#     group = rep(df_g$group, times = length(IFN_all_genes)),
#     gene  = as.character(st$ind),
#     expr  = as.numeric(st$values),
#     stringsAsFactors = FALSE
#   )
#
#   # annotate IFN class
#   df_long_g$IFN_class <- "Shared/Other"
#   df_long_g$IFN_class[df_long_g$gene %in% IFN_I_genes_present] <- "IFN-I (ISG)"
#   df_long_g$IFN_class[df_long_g$gene %in% IFNG_genes_present]  <- "IFN-γ"
#
#   # order genes by biology (NOT alphabetic)
#   gene_order <- c(
#     b_intersect(IFN_I_genes, IFN_all_genes),
#     b_intersect(b_setdiff(IFNG_genes, "Stat1"), IFN_all_genes)
#   )
#   # ensure unique and keep present
#   gene_order <- unique(gene_order)
#   df_long_g$gene <- factor(df_long_g$gene, levels = gene_order)
#
#   utils::write.csv(df_long_g, file.path(out_ifn, "IFN_coreGenes_singlecell_expr_long.csv"), row.names = FALSE)
#
#   # per-gene violin PDFs
#   for (gn in gene_order) {
#     subd <- df_long_g[df_long_g$gene == gn, , drop = FALSE]
#     if (nrow(subd) == 0) next
#     p <- ggplot(subd, aes(x = group, y = expr)) +
#       geom_violin(trim = TRUE, scale = "width") +
#       geom_boxplot(width = 0.15, outlier.size = 0.2) +
#       theme_bw() +
#       labs(title = paste0("eTreg-like | IFN core gene: ", gn),
#            x = "Group", y = "Expression (FetchData: RNA)")
#     ggsave(file.path(out_ifn, paste0("IFN_gene_", gn, "_violin.pdf")), p, width = 8, height = 5)
#   }
#
#   # dotplot-like summary: mean expr + pct>0 (one combined figure)
#   df_long_g$expr_pos <- df_long_g$expr > 0
#   agg_mean <- aggregate(expr ~ group + gene + IFN_class, data = df_long_g, FUN = mean)
#   agg_pct  <- aggregate(expr_pos ~ group + gene + IFN_class, data = df_long_g, FUN = mean)
#   colnames(agg_pct)[4] <- "pct_pos"
#   df_dot <- merge(agg_mean, agg_pct, by = c("group","gene","IFN_class"), all = TRUE)
#
#   # color gradient (low-mid-high) as you prefer
#   low_col <- "#F2F2F2"
#   mid_col <- "#A5C6DC"
#   high_col <- "#F4C09B"
#   midpoint_val <- stats::median(df_dot$expr, na.rm = TRUE)
#
#   p_dot <- ggplot(df_dot, aes(x = group, y = gene)) +
#     geom_point(aes(size = pct_pos, color = expr)) +
#     scale_color_gradient2(low = low_col, mid = mid_col, high = high_col, midpoint = midpoint_val) +
#     theme_bw() +
#     labs(title = "eTreg-like | IFN core genes (mean expr & %>0)",
#          x = "Group", y = "Gene", color = "Mean expr", size = "%>0")
#
#   ggsave(file.path(out_ifn, "IFN_coreGenes_dotplot_like.pdf"), p_dot, width = 10, height = 7)
#   utils::write.csv(df_dot, file.path(out_ifn, "IFN_coreGenes_dot_summary_mean_pct.csv"), row.names = FALSE)
#
#   # ---- 2) IFN score vs Prolif/IL2 correlation (same eTreg-like single cells) ----
#   # We compute IFN scores by AddModuleScore (robust & self-contained),
#   # then correlate with your internal module scores (Prolif_CC_MS1, IL2_axis_MS4) if present.
#
#   IFN_sets <- list(
#     IFN_I_CORE = IFN_I_genes_present,
#     IFNG_CORE  = IFNG_genes_present
#   )
#   # require >=5 genes for stability if possible, else allow >=3
#   lens <- vapply(IFN_sets, length, integer(1))
#   keep_sets <- names(IFN_sets)[lens >= 3]
#   IFN_sets <- IFN_sets[keep_sets]
#
#   if (length(IFN_sets) == 0) {
#     cat("[EXT++/IFN] No IFN set has >=3 genes matched. Skip correlation block.\n")
#   } else {
#
#     # AddModuleScore returns columns like "IFN_I_CORE_MS1" etc
#     prefixes <- paste0(names(IFN_sets), "_MS")
#     et_ifn <- Seurat::AddModuleScore(
#       object   = et_focus,
#       features = base::unname(IFN_sets),
#       name     = prefixes,
#       assay    = "RNA",
#       search   = FALSE
#     )
#
#     meta_cols <- colnames(et_ifn@meta.data)
#
#     find_ms_col <- function(meta_cols, prefix) {
#       hit <- grep(paste0("^", prefix, "\\d+$"), meta_cols, value = TRUE)
#       if (length(hit) == 0) return(NA_character_)
#       hit[length(hit)]
#     }
#     ifn_score_cols <- vapply(prefixes, function(pf) find_ms_col(meta_cols, pf), character(1))
#
#     if (any(is.na(ifn_score_cols))) {
#       cat("\n[EXT++/IFN][DEBUG] meta.data cols ending with _MS#:\n")
#       print(grep("_MS\\d+$", meta_cols, value = TRUE))
#       stop("[EXT++/IFN] Failed to match IFN AddModuleScore columns for: ",
#            paste(prefixes[is.na(ifn_score_cols)], collapse = ", "))
#     }
#
#     df_map <- data.frame(
#       set = names(IFN_sets),
#       prefix = prefixes,
#       score_col = ifn_score_cols,
#       n_genes = vapply(IFN_sets, length, integer(1)),
#       genes = vapply(IFN_sets, function(x) paste(x, collapse = ";"), character(1)),
#       stringsAsFactors = FALSE
#     )
#     utils::write.csv(df_map, file.path(out_ifn, "IFN_signature_scorecol_map.csv"), row.names = FALSE)
#
#     # gather wide table for correlation
#     md <- et_ifn@meta.data
#     df_sc <- data.frame(
#       cell = rownames(md),
#       group = as.character(md[[group_col]]),
#       stringsAsFactors = FALSE
#     )
#     if (exists("wanted_groups")) df_sc <- df_sc[df_sc$group %in% wanted_groups, , drop = FALSE]
#
#     # attach IFN scores
#     for (i in seq_along(names(IFN_sets))) {
#       nm <- names(IFN_sets)[i]
#       df_sc[[nm]] <- md[df_sc$cell, ifn_score_cols[i]]
#     }
#
#     # attach Prolif/IL2 scores if present (your internal module scores)
#     # prefer Prolif_CC_MS1 & IL2_axis_MS4 if they exist; otherwise try fuzzy match.
#     get_col_or_na <- function(md, candidates) {
#       hit <- candidates[candidates %in% colnames(md)]
#       if (length(hit) == 0) return(NA_character_)
#       hit[1]
#     }
#
#     prolif_col <- get_col_or_na(md, c("Prolif_CC_MS1"))
#     il2_col    <- get_col_or_na(md, c("IL2_axis_MS4"))
#
#     if (is.na(prolif_col) || is.na(il2_col)) {
#       cat("[EXT++/IFN] Prolif/IL2 module score columns not found (need Prolif_CC_MS1 and IL2_axis_MS4). Skip correlation plots.\n")
#     } else {
#       df_sc$Prolif_CC <- md[df_sc$cell, prolif_col]
#       df_sc$IL2_axis  <- md[df_sc$cell, il2_col]
#
#       # correlations (Spearman): overall + per-group
#       cor_list <- list()
#       pairs <- list(
#         c("IFN_I_CORE", "Prolif_CC"),
#         c("IFN_I_CORE", "IL2_axis"),
#         c("IFNG_CORE",  "Prolif_CC"),
#         c("IFNG_CORE",  "IL2_axis")
#       )
#
#       do_cor <- function(dsub, x, y) {
#         ok <- is.finite(dsub[[x]]) & is.finite(dsub[[y]])
#         if (sum(ok) < 20) return(NULL)
#         ct <- stats::cor.test(dsub[[x]][ok], dsub[[y]][ok], method = "spearman")
#         data.frame(
#           x = x, y = y,
#           n = sum(ok),
#           rho = unname(ct$estimate),
#           pval = ct$p.value,
#           stringsAsFactors = FALSE
#         )
#       }
#
#       # overall
#       for (pp in pairs) {
#         tmp <- do_cor(df_sc, pp[1], pp[2])
#         if (!is.null(tmp)) {
#           tmp$scope <- "ALL"
#           cor_list[[length(cor_list)+1]] <- tmp
#         }
#       }
#       # by group
#       for (g in unique(df_sc$group)) {
#         dsub <- df_sc[df_sc$group == g, , drop = FALSE]
#         for (pp in pairs) {
#           tmp <- do_cor(dsub, pp[1], pp[2])
#           if (!is.null(tmp)) {
#             tmp$scope <- g
#             cor_list[[length(cor_list)+1]] <- tmp
#           }
#         }
#       }
#
#       df_cor <- do.call(rbind, cor_list)
#       if (!is.null(df_cor) && nrow(df_cor) > 0) {
#         df_cor$padj <- stats::p.adjust(df_cor$pval, method = "BH")
#         utils::write.csv(df_cor, file.path(out_ifn, "IFN_vs_Prolif_IL2_spearman.csv"), row.names = FALSE)
#       }
#
#       # correlation scatter PDFs (facet by group and pair)
#       # build long for plotting
#       plot_long <- list()
#       for (pp in pairs) {
#         x <- pp[1]; y <- pp[2]
#         tmp <- data.frame(
#           group = df_sc$group,
#           x_name = x,
#           y_name = y,
#           x = df_sc[[x]],
#           y = df_sc[[y]],
#           stringsAsFactors = FALSE
#         )
#         plot_long[[length(plot_long)+1]] <- tmp
#       }
#       df_plot <- do.call(rbind, plot_long)
#       if (exists("wanted_groups")) df_plot$group <- factor(df_plot$group, levels = wanted_groups)
#
#       p_cor <- ggplot(df_plot, aes(x = x, y = y)) +
#         geom_point(alpha = 0.35, size = 0.8) +
#         geom_smooth(method = "lm", se = FALSE) +
#         theme_bw() +
#         facet_grid(y_name ~ group, scales = "free_y") +
#         labs(title = "eTreg-like | IFN score vs Prolif/IL2 module scores (single-cell)",
#              x = "IFN signature score", y = "Target module score")
#
#       ggsave(file.path(out_ifn, "IFN_vs_Prolif_IL2_correlation.pdf"), p_cor, width = 12, height = 8)
#
#       utils::write.csv(df_sc, file.path(out_ifn, "IFN_signature_scores_wide.csv"), row.names = FALSE)
#     }
#
#     cat("[EXT++/IFN] Outputs: ", out_ifn, "\n", sep="")
#   }
# }
#
# # ============================================================
# # [EXT++/STAT5] STAT5 downstream response signatures (single-cell)
# #   - Robust alternative to pSTAT5 at transcriptomic level
# #   - Avoid Seurat::VlnPlot (tidy-eval / S4SXP issue); use ggplot2 directly
# #   Output: ./temp/Step1.7_CLOSURE_EXT/SC_STAT5/
# # ============================================================
#
# cat("\n[EXT++/STAT5] Single-cell STAT5 downstream signature scoring ...\n")
#
# # ---- prerequisites: et_focus must exist (eTreg-like object), group_col must exist
# if (!exists("et_focus")) stop("[EXT++/STAT5] et_focus not found. Make sure you already created eTreg-like focus object in EXT++.")
# if (!exists("group_col")) group_col <- "group"
#
# outdir_stat5 <- "./temp/Step1.7_CLOSURE_EXT/SC_STAT5"
# dir.create(outdir_stat5, showWarnings = FALSE, recursive = TRUE)
#
# # ---- explicit base/utils/stats to avoid conflicted errors
# base_intersect <- base::intersect
# base_setdiff   <- base::setdiff
# base_table     <- base::table
# base_names     <- base::names
#
# # ---- 1) define STAT5 signatures (mouse symbols; adjust if your data uses different casing)
# STAT5_sets <- list(
#   STAT5_CORE = c(
#     "Il2ra","Il2rb","Il2rg","Stat5a","Stat5b",
#     "Bcl2","Mcl1","Foxp3","Ctla4"
#   ),
#   STAT5_FEEDBACK = c(
#     "Cish","Socs1","Socs2","Socs3","Pim1","Pim2","Bcl2","Il2ra"
#   ),
#   # "TREG_DOWN" here means canonical Treg effector/stability program that often tracks with IL2-STAT5 engagement;
#   # rename if you prefer; it is NOT a "downregulated geneset" literally.
#   STAT5_TREG_DOWN = c(
#     "Ikzf2","Tigit","Icos","Tnfrsf18","Tnfrsf4","Il2ra","Ctla4","Maf","Batf","Foxp3"
#   )
# )
#
# # ---- 2) match genes present; require >=5 genes to keep a set
# STAT5_sets_present <- lapply(STAT5_sets, function(gs) base_intersect(gs, rownames(et_focus)))
# lens <- vapply(STAT5_sets_present, length, integer(1))
# STAT5_sets_present <- STAT5_sets_present[lens >= 5]
#
# if (length(STAT5_sets_present) == 0) {
#   stop("[EXT++/STAT5] No STAT5 gene sets have >=5 matched genes. Check gene symbol naming (mouse vs human, Ensembl vs Symbol).")
# }
#
# df_map <- data.frame(
#   set = base_names(STAT5_sets_present),
#   n_genes_matched = vapply(STAT5_sets_present, length, integer(1)),
#   genes = vapply(STAT5_sets_present, function(x) paste(x, collapse = ";"), character(1)),
#   stringsAsFactors = FALSE
# )
# utils::write.csv(df_map, file = file.path(outdir_stat5, "STAT5_sets_matched.csv"), row.names = FALSE)
#
# # ---- 3) AddModuleScore (safe naming)
# # Important: AddModuleScore returns columns like "STAT5_CORE_MS1" etc.
# prefixes <- paste0(base_names(STAT5_sets_present), "_MS")
# et_stat5 <- Seurat::AddModuleScore(
#   object   = et_focus,
#   features = base::unname(STAT5_sets_present),
#   name     = prefixes,
#   assay    = "RNA",
#   search   = FALSE
# )
#
# # ---- 4) locate generated score columns robustly
# meta_cols <- colnames(et_stat5@meta.data)
#
# find_ms_col <- function(meta_cols, prefix) {
#   hit <- grep(paste0("^", prefix, "\\d+$"), meta_cols, value = TRUE)
#   if (length(hit) == 0) return(NA_character_)
#   hit[length(hit)]
# }
#
# score_cols <- vapply(prefixes, function(pf) find_ms_col(meta_cols, pf), character(1))
# if (any(is.na(score_cols))) {
#   cat("\n[EXT++/STAT5][DEBUG] meta.data cols ending with _MS#:\n")
#   print(grep("_MS\\d+$", meta_cols, value = TRUE))
#   stop("[EXT++/STAT5] Failed to match AddModuleScore output columns for: ",
#        paste(prefixes[is.na(score_cols)], collapse = ", "))
# }
#
# df_scorecol_map <- data.frame(
#   set = base_names(STAT5_sets_present),
#   prefix = prefixes,
#   score_col = score_cols,
#   stringsAsFactors = FALSE
# )
# utils::write.csv(df_scorecol_map, file = file.path(outdir_stat5, "STAT5_scorecol_map.csv"), row.names = FALSE)
#
# # ---- 5) build long table (no tidyr) for ggplot
# md <- et_stat5@meta.data
# grp <- as.character(md[[group_col]])
# names(grp) <- rownames(md)
#
# # keep only wanted groups if you have them defined, else keep all
# if (exists("wanted_groups")) {
#   keep_cells <- rownames(md)[grp %in% wanted_groups]
# } else {
#   keep_cells <- rownames(md)
# }
#
# # Fetch scores + group into a plain data.frame (avoid tidy-eval)
# df_wide <- data.frame(
#   cell  = keep_cells,
#   group = grp[keep_cells],
#   stringsAsFactors = FALSE
# )
# for (i in seq_along(score_cols)) {
#   cn <- base_names(STAT5_sets_present)[i]
#   sc <- score_cols[i]
#   df_wide[[cn]] <- as.numeric(md[keep_cells, sc])
# }
#
# # Convert to long format via utils::stack  (NOTE: stack() is in utils, not base)
# mat_scores <- df_wide[, base_names(STAT5_sets_present), drop = FALSE]
#
# if (ncol(mat_scores) == 0) {
#   stop("[EXT++/STAT5] No STAT5 signature score columns found in df_wide. ",
#        "Check STAT5_sets_present and the columns created by AddModuleScore.")
# }
#
# st <- utils::stack(mat_scores)  # returns: values, ind
#
# df_long <- data.frame(
#   cell      = rep(df_wide$cell,  times = ncol(mat_scores)),
#   group     = rep(df_wide$group, times = ncol(mat_scores)),
#   signature = as.character(st$ind),
#   score     = as.numeric(st$values),
#   stringsAsFactors = FALSE
# )
#
# utils::write.csv(df_long, file = file.path(outdir_stat5, "STAT5_signature_scores_long.csv"), row.names = FALSE)
#
# # ---- 6) cell-level Wilcoxon tests (descriptive; treat with caution if n is small)
# # contrasts: control1 vs s1_hIL, control1 vs s40kD, s1_hIL vs s40kD if present
# all_groups <- sort(unique(df_long$group))
# make_pairs <- function(gs) {
#   # ordered preference if those groups exist
#   pref <- list(c("control1","s1_hIL"), c("control1","s40kD"), c("s1_hIL","s40kD"))
#   out <- list()
#   for (pp in pref) if (all(pp %in% gs)) out[[length(out)+1]] <- pp
#   out
# }
# pairs <- make_pairs(all_groups)
#
# wilcox_one <- function(sig, g1, g2) {
#   x <- df_long$score[df_long$signature == sig & df_long$group == g1]
#   y <- df_long$score[df_long$signature == sig & df_long$group == g2]
#   if (length(x) < 5 || length(y) < 5) {
#     return(data.frame(signature=sig, group1=g1, group2=g2, n1=length(x), n2=length(y),
#                       p=NA_real_, stringsAsFactors=FALSE))
#   }
#   p <- tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
#   data.frame(signature=sig, group1=g1, group2=g2, n1=length(x), n2=length(y),
#              p=p, stringsAsFactors=FALSE)
# }
#
# df_wcx <- do.call(rbind, lapply(unique(df_long$signature), function(sig) {
#   if (length(pairs) == 0) return(NULL)
#   do.call(rbind, lapply(pairs, function(pp) wilcox_one(sig, pp[1], pp[2])))
# }))
#
# if (!is.null(df_wcx) && nrow(df_wcx) > 0) {
#   # BH adjust within each signature
#   df_wcx$padj <- NA_real_
#   for (sig in unique(df_wcx$signature)) {
#     idx <- which(df_wcx$signature == sig & !is.na(df_wcx$p))
#     if (length(idx) > 0) df_wcx$padj[idx] <- stats::p.adjust(df_wcx$p[idx], method = "BH")
#   }
#   utils::write.csv(df_wcx, file = file.path(outdir_stat5, "STAT5_signature_wilcox_celllevel.csv"), row.names = FALSE)
# }
#
# # ---- 7) PDF plot (violin + box); no png output
# p_stat5 <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = score)) +
#   ggplot2::geom_violin(trim = TRUE, scale = "width") +
#   ggplot2::geom_boxplot(width = 0.18, outlier.size = 0.2) +
#   ggplot2::facet_wrap(~signature, scales = "free_y") +
#   ggplot2::theme_bw() +
#   ggplot2::labs(
#     title = "eTreg-like | STAT5 downstream response signatures (single-cell scores)",
#     x = "Group",
#     y = "Signature score (AddModuleScore)"
#   )
#
# ggplot2::ggsave(
#   filename = file.path(outdir_stat5, "STAT5_signature_scores_violin_box.pdf"),
#   plot = p_stat5,
#   width = 12, height = 5
# )
#
# cat("[EXT++/STAT5] Outputs:\n",
#     "  - ", outdir_stat5, "/STAT5_sets_matched.csv\n",
#     "  - ", outdir_stat5, "/STAT5_scorecol_map.csv\n",
#     "  - ", outdir_stat5, "/STAT5_signature_scores_long.csv\n",
#     "  - ", outdir_stat5, "/STAT5_signature_wilcox_celllevel.csv\n",
#     "  - ", outdir_stat5, "/STAT5_signature_scores_violin_box.pdf\n",
#     sep = "")
#
# # ============================================================
# # [EXT++/IL2-STAT5-COR] Correlate IL2RA / STAT5-downstream with Cish/Socs feedback
# #   - per-cell Spearman correlation within eTreg-like (et_focus)
# #   - group-wise Fisher-z comparison (DESCRIPTIVE; single-cell pseudo-replication caveat)
# #   Output: ./temp/Step1.7_CLOSURE_EXT/SC_STAT5_Correlation/
# # ============================================================
#
# cat("\n[EXT++/IL2-STAT5-COR] Correlating IL2RA/STAT5-downstream with Cish/Socs feedback ...\n")
#
# out_cor <- file.path(outdir_closure, "SC_STAT5_Correlation")
# dir.create(out_cor, showWarnings = FALSE, recursive = TRUE)
#
# # ---- resolve gene symbols (mouse vs human) ----
# pick_gene <- function(cands, rn) {
#   for (g in cands) if (g %in% rn) return(g)
#   return(NA_character_)
# }
# g_Il2ra <- pick_gene(c("Il2ra","IL2RA"), rownames(et_focus))
# g_Cish  <- pick_gene(c("Cish","CISH"), rownames(et_focus))
# g_Socs1 <- pick_gene(c("Socs1","SOCS1"), rownames(et_focus))
# g_Socs3 <- pick_gene(c("Socs3","SOCS3"), rownames(et_focus))
#
# # STAT5 feedback score column (preferred) from EXT++/STAT5 block
# col_stat5_core     <- NA_character_
# col_stat5_feedback <- NA_character_
# col_stat5_tregdown <- NA_character_
# if (exists("score_map")) {
#   if ("STAT5_core"     %in% score_map$signature) col_stat5_core     <- score_map$score_col[match("STAT5_core",     score_map$signature)]
#   if ("STAT5_feedback" %in% score_map$signature) col_stat5_feedback <- score_map$score_col[match("STAT5_feedback", score_map$signature)]
#   if ("STAT5_Treg_down"%in% score_map$signature) col_stat5_tregdown <- score_map$score_col[match("STAT5_Treg_down",score_map$signature)]
# }
#
# # STAT5-downstream alternative: use IL2_axis_MS4 if present; else STAT5_core score
# col_il2axis <- if ("IL2_axis_MS4" %in% colnames(et_focus@meta.data)) "IL2_axis_MS4" else NA_character_
# col_stat5_down <- if (!is.na(col_il2axis)) col_il2axis else col_stat5_core
#
# # ---- assemble per-cell table ----
# vars_meta <- c(group_col, col_stat5_down, col_stat5_feedback, col_stat5_core, col_stat5_tregdown)
# vars_meta <- vars_meta[!is.na(vars_meta)]
# vars_meta <- unique(vars_meta)
#
# df_meta <- et_focus@meta.data[, vars_meta, drop = FALSE]
# df_meta$cell <- rownames(df_meta)
# df_meta$group <- as.character(df_meta[[group_col]])
#
# df_meta <- df_meta[df_meta$group %in% wanted_groups, , drop = FALSE]
# df_meta$group <- factor(df_meta$group, levels = wanted_groups)
#
# # gene expression pull
# genes_to_pull <- c(g_Il2ra, g_Cish, g_Socs1, g_Socs3)
# genes_to_pull <- genes_to_pull[!is.na(genes_to_pull)]
# if (length(genes_to_pull) > 0) {
#   df_g <- Seurat::FetchData(et_focus, vars = genes_to_pull)
#   df_g$cell <- rownames(df_g)
#   df_meta <- merge(df_meta, df_g, by = "cell", all.x = TRUE)
# }
#
# # define x/y for correlation
# df_meta$IL2RA_expr <- if (!is.na(g_Il2ra) && g_Il2ra %in% colnames(df_meta)) as.numeric(df_meta[[g_Il2ra]]) else NA_real_
#
# if (!is.na(col_stat5_feedback) && col_stat5_feedback %in% colnames(df_meta)) {
#   df_meta$STAT5_feedback_score <- as.numeric(df_meta[[col_stat5_feedback]])
# } else {
#   fb_genes <- c(g_Cish, g_Socs1, g_Socs3)
#   fb_genes <- fb_genes[!is.na(fb_genes)]
#   if (length(fb_genes) > 0) {
#     df_meta$STAT5_feedback_score <- rowMeans(df_meta[, fb_genes, drop = FALSE], na.rm = TRUE)
#   } else {
#     df_meta$STAT5_feedback_score <- NA_real_
#   }
# }
#
# df_meta$STAT5_down_score <- if (!is.na(col_stat5_down) && col_stat5_down %in% colnames(df_meta)) as.numeric(df_meta[[col_stat5_down]]) else NA_real_
#
# df_cor <- df_meta[is.finite(df_meta$STAT5_feedback_score) & (is.finite(df_meta$STAT5_down_score) | is.finite(df_meta$IL2RA_expr)), , drop = FALSE]
#
# utils::write.csv(df_cor[, c("cell","group","IL2RA_expr","STAT5_down_score","STAT5_feedback_score")],
#                  file.path(out_cor, "IL2RA_STAT5down_STAT5feedback_perCell.csv"),
#                  row.names = FALSE)
#
# spearman_by_group <- function(df, xcol, ycol) {
#   out <- list()
#   for (g in base::levels(df$group)) {
#     d <- df[df$group == g &
#               base::is.finite(df[[xcol]]) &
#               base::is.finite(df[[ycol]]), , drop = FALSE]
#     n <- base::nrow(d)
#     if (n < 10) next
#     ct <- base::suppressWarnings(stats::cor.test(d[[xcol]], d[[ycol]],
#                                                  method = "spearman", exact = FALSE))
#     out[[base::length(out) + 1]] <- base::data.frame(
#       group = g, x = xcol, y = ycol,
#       n = n, rho = base::unname(ct$estimate), pval = ct$p.value,
#       stringsAsFactors = FALSE
#     )
#   }
#   if (base::length(out) == 0) return(base::data.frame())
#   base::do.call(base::rbind, out)
# }
#
# res1 <- spearman_by_group(df_cor, "STAT5_down_score", "STAT5_feedback_score")
# res2 <- spearman_by_group(df_cor, "IL2RA_expr",       "STAT5_feedback_score")
# res_all <- rbind(res1, res2)
# utils::write.csv(res_all, file.path(out_cor, "Spearman_byGroup.csv"), row.names = FALSE)
#
# fisher_z_pairwise <- function(res_df) {
#   if (nrow(res_df) == 0) return(data.frame())
#   pairs <- list()
#   for (i in 1:(nrow(res_df)-1)) {
#     for (j in (i+1):nrow(res_df)) {
#       if (res_df$x[i] != res_df$x[j] || res_df$y[i] != res_df$y[j]) next
#       r1 <- res_df$rho[i]; r2 <- res_df$rho[j]
#       n1 <- res_df$n[i];   n2 <- res_df$n[j]
#       if (abs(r1) >= 0.999 || abs(r2) >= 0.999 || n1 < 10 || n2 < 10) next
#       z1 <- atanh(r1); z2 <- atanh(r2)
#       z  <- (z1 - z2) / sqrt(1/(n1-3) + 1/(n2-3))
#       p  <- 2 * stats::pnorm(-abs(z))
#       pairs[[length(pairs)+1]] <- data.frame(
#         x = res_df$x[i], y = res_df$y[i],
#         group1 = res_df$group[i], group2 = res_df$group[j],
#         rho1 = r1, n1 = n1, rho2 = r2, n2 = n2,
#         z = z, pval = p,
#         stringsAsFactors = FALSE
#       )
#     }
#   }
#   if (length(pairs) == 0) return(data.frame())
#   do.call(rbind, pairs)
# }
# pair_z <- fisher_z_pairwise(res_all)
# utils::write.csv(pair_z, file.path(out_cor, "FisherZ_pairwise_DESCRIPTIVE.csv"), row.names = FALSE)
#
# pal3 <- setNames(c("#BDBDBD", "#A5C6DC", "#F4C09B"), wanted_groups)
#
# plot_scatter_facet <- function(df, xcol, ycol, out_pdf, title) {
#   dd <- df[df$group %in% wanted_groups & is.finite(df[[xcol]]) & is.finite(df[[ycol]]), , drop=FALSE]
#   if (nrow(dd) < 50) return(invisible(NULL))
#   p <- ggplot(dd, aes_string(x = xcol, y = ycol)) +
#     geom_point(aes(color = group), size = 0.6, alpha = 0.5) +
#     geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
#     facet_wrap(~ group, scales = "free") +
#     scale_color_manual(values = pal3, drop = FALSE) +
#     theme_bw() +
#     labs(title = title, x = xcol, y = ycol)
#   ggsave(out_pdf, p, width = 10, height = 3.6)
# }
#
# plot_scatter_facet(df_cor, "STAT5_down_score", "STAT5_feedback_score",
#                    file.path(out_cor, "Scatter_STAT5down_vs_STAT5feedback.pdf"),
#                    "eTreg-like | STAT5-downstream vs Cish/Socs feedback (per-cell; descriptive)")
#
# plot_scatter_facet(df_cor, "IL2RA_expr", "STAT5_feedback_score",
#                    file.path(out_cor, "Scatter_IL2RAexpr_vs_STAT5feedback.pdf"),
#                    "eTreg-like | IL2RA expression vs Cish/Socs feedback (per-cell; descriptive)")
#
# cat("[EXT++/IL2-STAT5-COR] Outputs: ", out_cor, "\n", sep = "")
#
# # ============================================================
# # 3.2 Key regulators panel: Violin/Dot + per-cell mean logFC (descriptive)
# #   - avoid Seurat::VlnPlot/DotPlot; use FetchData + ggplot2
# # ============================================================
#
# key_gene_sets <- list(
#   IL2_axis = c("Il2ra","Il2rb","Il2rg","Stat5a","Stat5b","Socs1","Socs3","Cish"),
#   eTreg_core = c("Foxp3","Ikzf2","Ctla4","Icos","Tnfrsf4","Tnfrsf18","Batf","Maf","Tigit","Lag3","Pdcd1"),
#   CellCycle = c("Mki67","Top2a","Mcm5","Mcm6","Mcm7","Pcna","Cenpf","Cdk1","Birc5"),
#   Homing = c("Ccr7","Ccr4","Ccr8","Cxcr4","Sell","Itga4","S1pr1"),
#   Survival_anti = c("Bcl2","Mcl1","Bcl2l1"),
#   Death_pro = c("Bax","Bak1","Fas","Casp3","Casp8")
# )
#
# key_genes <- unique(unlist(key_gene_sets))
# key_genes <- b_intersect(key_genes, rownames(et_focus))
# if (length(key_genes) == 0) {
#   cat("\n[EXT++/KEY] No key genes found in et_focus (symbol mismatch). Skip key-gene panel.\n")
# } else {
#   df_g <- Seurat::FetchData(et_focus, vars = c(group_col, key_genes))
#   df_g$group <- as.character(df_g[[group_col]])
#   df_g <- df_g[df_g$group %in% wanted_groups, , drop = FALSE]
#
#   # long format (utils::stack)
#   df_expr <- utils::stack(df_g[, key_genes, drop = FALSE])
#   df_long_g <- data.frame(
#     group = rep(df_g$group, times = length(key_genes)),
#     gene  = as.character(df_expr$ind),
#     expr  = as.numeric(df_expr$values),
#     stringsAsFactors = FALSE
#   )
#   df_long_g$group <- factor(df_long_g$group, levels = wanted_groups)
#
#   utils::write.csv(df_long_g, file.path(out_key, "KeyGenes_singlecell_expr_long.csv"), row.names = FALSE)
#
#   # per gene violin
#   for (gn in unique(df_long_g$gene)) {
#     subd <- df_long_g[df_long_g$gene == gn, , drop = FALSE]
#     p <- ggplot(subd, aes(x = group, y = expr)) +
#       geom_violin(trim = TRUE, scale = "width") +
#       geom_boxplot(width = 0.15, outlier.size = 0.2) +
#       theme_bw() +
#       labs(title = paste0("eTreg-like | ", gn), x = "Group", y = "Expression (FetchData: RNA data)")
#     ggsave(file.path(out_key, paste0("KeyGene_", gn, "_violin.pdf")), p, width = 8, height = 5)
#   }
#
#   # DotPlot-like summary: mean expr + pct>0
# # ---- 1) gene order by function blocks ----
# gene_order <- c(
#   key_gene_sets$IL2_axis,
#   key_gene_sets$eTreg_core,
#   key_gene_sets$CellCycle,
#   key_gene_sets$Homing,
#   key_gene_sets$Survival_anti,
#   key_gene_sets$Death_pro
# )
# # keep only genes present in et_focus (key_genes already intersected)
# gene_order <- gene_order[gene_order %in% unique(df_long_g$gene)]
# # de-duplicate but keep first occurrence
# gene_order <- gene_order[!duplicated(gene_order)]
#
# # Optional: add a "module" column for faceting / annotation (useful in figure legend)
# gene2module <- data.frame(
#   gene = unlist(key_gene_sets),
#   module = rep(names(key_gene_sets), times = vapply(key_gene_sets, length, integer(1))),
#   stringsAsFactors = FALSE
# )
# gene2module <- gene2module[gene2module$gene %in% gene_order, , drop = FALSE]
#
# # ---- 2) Dot summary: mean expr + pct>0 ----
# df_long_g$expr_pos <- df_long_g$expr > 0
#
# agg_mean <- aggregate(expr ~ group + gene, data = df_long_g, FUN = mean)
# agg_pct  <- aggregate(expr_pos ~ group + gene, data = df_long_g, FUN = mean)
# colnames(agg_pct)[3] <- "pct_pos"
#
# df_dot <- merge(agg_mean, agg_pct, by = c("group","gene"), all = TRUE)
# df_dot$group <- factor(df_dot$group, levels = wanted_groups)
#
# # attach module (optional; not required for ordering)
# df_dot <- merge(df_dot, gene2module, by = "gene", all.x = TRUE)
#
# # enforce gene ordering (functional)
# df_dot$gene <- factor(df_dot$gene, levels = rev(gene_order))  # rev: put first group on top
#
# # ---- 3) Export tables (unchanged) ----
# wide <- reshape(df_dot[, c("group","gene","expr")], idvar = "gene", timevar = "group", direction = "wide")
#
# if (all(c("control1","s1_hIL","s40kD") %in% wanted_groups)) {
#   if ("expr.control1" %in% colnames(wide) && "expr.s1_hIL" %in% colnames(wide)) {
#     wide$logFC_s1_hIL_vs_control1 <- wide$expr.s1_hIL - wide$expr.control1
#   }
#   if ("expr.control1" %in% colnames(wide) && "expr.s40kD" %in% colnames(wide)) {
#     wide$logFC_s40kD_vs_control1 <- wide$expr.s40kD - wide$expr.control1
#   }
#   if ("expr.s1_hIL" %in% colnames(wide) && "expr.s40kD" %in% colnames(wide)) {
#     wide$logFC_s1_hIL_vs_s40kD <- wide$expr.s1_hIL - wide$expr.s40kD
#   }
# }
#
# utils::write.csv(df_dot, file.path(out_key, "KeyGenes_dot_summary_mean_pct.csv"), row.names = FALSE)
# utils::write.csv(wide,   file.path(out_key, "KeyGenes_groupMean_logFC_descriptive.csv"), row.names = FALSE)
#
# # ---- 4) Dot plot with 3-stop continuous gradient (low/mid/high) ----
# # Choose mid-point for the gradient: median of mean expression (robust & automatic)
# mid_point <- stats::median(df_dot$expr, na.rm = TRUE)
#
# p_dot <- ggplot(df_dot, aes(x = group, y = gene)) +
#   geom_point(aes(size = pct_pos, color = expr), alpha = 0.95) +
#   scale_size_continuous(
#     range = c(0.5, 5),
#     breaks = c(0.1, 0.3, 0.5, 0.7),
#     labels = scales::percent_format(accuracy = 1),
#     name = "% cells\n(expr > 0)"
#   ) +
#   scale_color_gradient2(
#     low  = "#F2F2F2",
#     mid  = "#A5C6DC",
#     high = "#F4C09B",
#     midpoint = mid_point,
#     name = "Mean expression\n(continuous gradient)"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.y = element_text(size = 9),
#     axis.text.x = element_text(size = 10),
#     legend.position = "right"
#   ) +
#   labs(
#     title = "eTreg-like | Key regulators (ordered by biological function)",
#     x = "Group",
#     y = NULL
#   )
#
# ggsave(file.path(out_key, "KeyGenes_dotplot_like.pdf"), p_dot, width = 10, height = 8)
#
# cat("\n[EXT++/KEY] Key gene dotplot written (functional order + 3-stop gradient): ",
#     file.path(out_key, "KeyGenes_dotplot_like.pdf"), "\n", sep = "")
#   cat("\n[EXT++/KEY] Key gene outputs: ", out_key, "\n", sep="")
# }
#
# # ============================================================
# # 3.4 IL-2 ligand source in ALL cells (optional but strong)
# #   - load combined_* if exists; summarize Il2 expression by celltype/cluster per group
# # ============================================================
#
# # Try load combined object
# if (!exists("combined")) {
#   cand_rds <- c(
#     "./temp/combined_after_marker_6.rds",
#     "./temp/combined_umap_5.rds",
#     "./temp/combined_normalization_4.rds",
#     "./temp/combined_raw_data_1.rds"
#   )
#   cand_rds <- cand_rds[file.exists(cand_rds)]
#   if (length(cand_rds) > 0) {
#     combined <- readRDS(cand_rds[1])
#     DefaultAssay(combined) <- "RNA"
#     cat("\n[EXT++/IL2] Loaded combined from: ", cand_rds[1], "\n", sep="")
#   }
# }
#
# detect_group_col2 <- function(md) {
#   cand <- c("group","orig.ident","sample","Sample","condition","treatment")
#   hit <- cand[cand %in% colnames(md)]
#   if (length(hit) == 0) return(NA_character_)
#   hit[1]
# }
# detect_celltype_col <- function(md) {
#   cand <- c("celltype","cell_type","CellType","annotation","Annotation","SingleR_label","SingleR.labels","predicted.id","seurat_clusters")
#   hit <- cand[cand %in% colnames(md)]
#   if (length(hit) == 0) return(NA_character_)
#   hit[1]
# }
#
# if (exists("combined")) {
#   gcol2 <- detect_group_col2(combined@meta.data)
#   ccol2 <- detect_celltype_col(combined@meta.data)
#
#   if (is.na(gcol2)) {
#     cat("\n[EXT++/IL2] combined has no group-like column; skip ligand source.\n")
#   } else if (!("Il2" %in% rownames(combined))) {
#     cat("\n[EXT++/IL2] combined has no Il2 gene (symbol mismatch). skip ligand source.\n")
#   } else {
#     if (is.na(ccol2)) {
#       # fallback: use seurat_clusters if exists else one bucket
#       if ("seurat_clusters" %in% colnames(combined@meta.data)) ccol2 <- "seurat_clusters"
#       else ccol2 <- NULL
#     }
#
#     vars <- c(gcol2, "Il2")
#     if (!is.null(ccol2)) vars <- c(vars, ccol2)
#
#     df_il2 <- Seurat::FetchData(combined, vars = vars)
#     df_il2$group <- as.character(df_il2[[gcol2]])
#     df_il2 <- df_il2[df_il2$group %in% wanted_groups, , drop = FALSE]
#     df_il2$Il2_pos <- df_il2$Il2 > 0
#
#     if (!is.null(ccol2)) {
#       df_il2$celltype <- as.character(df_il2[[ccol2]])
#     } else {
#       df_il2$celltype <- "ALL"
#     }
#
#     # summarize by (group, celltype)
#     agg_n <- aggregate(Il2 ~ group + celltype, data = df_il2, FUN = length)
#     colnames(agg_n)[3] <- "n_cells"
#     agg_pct <- aggregate(Il2_pos ~ group + celltype, data = df_il2, FUN = mean)
#     colnames(agg_pct)[3] <- "pct_Il2_pos"
#     agg_mean <- aggregate(Il2 ~ group + celltype, data = df_il2, FUN = mean)
#     colnames(agg_mean)[3] <- "mean_Il2"
#
#     df_sum <- merge(merge(agg_n, agg_pct, by=c("group","celltype")), agg_mean, by=c("group","celltype"))
#     utils::write.csv(df_sum, file.path(out_il2src, "Il2_ligand_source_summary_by_group_celltype.csv"), row.names = FALSE)
#
#     # plot: top celltypes by pct_Il2_pos per group (take top 15 each group)
#     top_list <- do.call(rbind, lapply(wanted_groups, function(g) {
#       subd <- df_sum[df_sum$group == g, , drop = FALSE]
#       subd <- subd[order(subd$pct_Il2_pos, decreasing = TRUE), , drop = FALSE]
#       head(subd, 15)
#     }))
#     top_list$group <- factor(top_list$group, levels = wanted_groups)
#
#     p_il2 <- ggplot(top_list, aes(x = reorder(celltype, pct_Il2_pos), y = pct_Il2_pos)) +
#       geom_col() +
#       coord_flip() +
#       facet_wrap(~group, scales = "free_y") +
#       theme_bw() +
#       labs(title = "IL-2 ligand source: % Il2+ by celltype (top15 per group)", x = "Celltype/cluster", y = "% Il2+")
#     ggsave(file.path(out_il2src, "Il2_ligand_source_top15_pctpos.pdf"), p_il2, width = 14, height = 7)
#
#     cat("\n[EXT++/IL2] IL-2 ligand-source outputs: ", out_il2src, "\n", sep="")
#   }
# } else {
#   cat("\n[EXT++/IL2] No combined_* found -> skip IL-2 ligand source (3.4).\n")
# }
#
# cat("\n[EXT++] Added 3.1/3.2/3.4 blocks DONE.\n")
#
#        # # ---------- 4. 输出发表用：UMAP + Treg marker 面板 + FeaturePlot + module scores ----------
# # # 4.1 UMAP（按 clean cluster）
# # p_umap <- DimPlot(treg_clean, reduction = "umap", group.by = clean_cluster_col, label = TRUE) +
# #   ggtitle("Strict Treg (clean) UMAP")
# # ggsave("./temp/Treg_clean_UMAP.pdf", p_umap, width = 11, height = 8)
# #
# # # 4.2 发表级 marker DotPlot（Treg core / activated / checkpoint / resting / cycling）
# # markers_pub <- c(
# #   "Foxp3","Ikzf2","Il2ra","Ctla4",
# #   "Icos","Tnfrsf18","Tnfrsf4","Tnfrsf9",
# #   "Pdcd1","Lag3","Tigit",
# #   "Tcf7","Lef1","Il7r",
# #   "Mki67","Top2a"
# # )
# # markers_pub <- base::intersect(markers_pub, rownames(treg_clean))
# #
# # p_dot_pub <- DotPlot(treg_clean, features = markers_pub, group.by = clean_cluster_col) +
# #   RotatedAxis() + ggtitle("Strict Treg (clean) markers (publication panel)")
# # ggsave("./temp/Treg_clean_markers_DotPlot.pdf", p_dot_pub, width = 18, height = 7)
# #
# # # 4.3 关键 FeaturePlot（空间分布核查）
# # fp_genes <- base::intersect(c("Foxp3","Il2ra","Ctla4","Pdcd1","Lag3","Tcf7","Il7r","Mki67"), rownames(treg_clean))
# # p_fp <- FeaturePlot(treg_clean, features = fp_genes, ncol = 4) +
# #   ggtitle("Strict Treg (clean) key markers (FeaturePlot)")
# # ggsave("./temp/Treg_clean_keymarkers_FeaturePlot.pdf", p_fp, width = 16, height = 10)
# #
# # # 4.4 module score（rTreg/eTreg/checkpoint/cycling）
# # treg_sets <- list(
# #   rTreg = c("Tcf7","Lef1","Il7r"),
# #   eTreg = c("Ctla4","Icos","Tnfrsf18","Tnfrsf4","Tnfrsf9"),
# #   Checkpoint = c("Pdcd1","Lag3","Tigit"),
# #   Cycling = c("Mki67","Top2a")
# # )
# # # 过滤掉不在数据里的基因
# # treg_sets <- lapply(treg_sets, function(g) base::intersect(g, rownames(treg_clean)))
# # treg_sets <- treg_sets[sapply(treg_sets, length) >= 1]
# #
# # for (nm in names(treg_sets)) {
# #   treg_clean <- AddModuleScore(treg_clean, features = list(treg_sets[[nm]]), name = paste0(nm, "_Score"))
# # }
# #
# # score_cols <- paste0(names(treg_sets), "_Score1")
# # p_scores <- VlnPlot(treg_clean, features = score_cols, group.by = clean_cluster_col, pt.size = 0, ncol = 2) +
# #   ggtitle("Strict Treg (clean) module scores by cluster")
# # ggsave("./temp/Treg_clean_module_scores.pdf", p_scores, width = 12, height = 8)
# #
# # cat("DONE. Outputs are saved under ./temp/\n")
#
#
#
#
# # p_umap <- DimPlot(treg_strict, label = TRUE) + ggtitle("Strict Treg UMAP (reclustered)")
# # ggsave("./temp/Treg_strict_UMAP.pdf", p_umap, width = 10, height = 8)
#
# # #一组“审稿人不太好挑刺”的 marker 面板（DotPlot + FeaturePlot）
# # markers_treg_pub <- base::intersect(
# #   c(
# #     # Treg core
# #     "Foxp3","Ikzf2","Il2ra","Ctla4",
# #     # activated / effector Treg
# #     "Tnfrsf18","Tnfrsf4","Tnfrsf9","Icos",
# #     # checkpoint-high Treg / dysfunction-like
# #     "Pdcd1","Lag3","Tigit",
# #     # resting/naive-like
# #     "Tcf7","Lef1","Il7r",
# #     # cycling
# #     "Mki67","Top2a"
# #   ),
# #   rownames(treg_strict)
# # )
# #
# # p_dot <- DotPlot(treg_strict, features = markers_treg_pub) + RotatedAxis() +
# #   ggtitle("Strict Treg markers (publication panel)")
# # ggsave("./temp/Treg_strict_markers_DotPlot.pdf", p_dot, width = 18, height = 7)
# # #彻底杜绝后续同类报错（推荐放脚本最前面一次）
# # conflicted::conflicts_prefer(base::intersect)
# # conflicted::conflicts_prefer(base::setdiff)
# # conflicted::conflicts_prefer(base::union)
# # conflicted::conflicts_prefer(base::unname)
# #
# # p_fp <- FeaturePlot(treg_strict, features = intersect(c("Foxp3","Il2ra","Ctla4","Pdcd1","Lag3","Mki67"), rownames(treg_strict)),
# #                     ncol = 3) + ggtitle("Strict Treg key markers (FeaturePlot)")
# # ggsave("./temp/Treg_strict_keymarkers_FeaturePlot.pdf", p_fp, width = 12, height = 8)
# #
# # # 用 module score 给亚群“定量命名”（可写入论文方法）
# # treg_sets2 <- list(
# #   rTreg = intersect(c("Tcf7","Lef1","Il7r"), rownames(treg_strict)),                  # resting/naive-like
# #   eTreg = intersect(c("Ctla4","Icos","Tnfrsf18","Tnfrsf4","Tnfrsf9"), rownames(treg_strict)), # activated/effector
# #   Checkpoint = intersect(c("Pdcd1","Lag3","Tigit"), rownames(treg_strict)),           # checkpoint-high
# #   Cycling = intersect(c("Mki67","Top2a"), rownames(treg_strict))                      # proliferating
# # )
# # treg_sets2 <- treg_sets2[sapply(treg_sets2, length) >= 1]
# #
# # for (nm in names(treg_sets2)) {
# #   treg_strict <- AddModuleScore(treg_strict, features = list(treg_sets2[[nm]]), name = paste0(nm, "_Score"))
# # }
# #
# # score_cols <- paste0(names(treg_sets2), "_Score1")
# # p_scores <- VlnPlot(treg_strict, features = score_cols, pt.size = 0, ncol = 2) +
# #   ggtitle("Strict Treg module scores by subclusters")
# # ggsave("./temp/Treg_strict_module_scores.pdf", p_scores, width = 12, height = 8)
#
# cat("（7.2） 完成【Treg分析】\n")

#
# # --------------------- 7.3. 开始分析CD8细胞 ---------------------
cat("7.3. 开始【结束分析CD8细胞】\n")

#
# # --------------------- 7. 使用SingleR进行细胞类型注释 ---------------------
# cat("7. 开始【使用SingleR进行细胞类型注释】\n")
# # 1）确保Seurat当前身份是cluster（否则WhichCells/levels会错）
# Idents(combined) <- "seurat_clusters"
# clusters <- levels(Idents(combined))
# # 加载参考数据集 (使用小鼠免疫细胞数据集作为示例，根据实际需要调整)
# ref_data <- celldex::ImmGenData()
# # 提取表达矩阵用于SingleR注释
# #用 counts 层，基因名最可靠；再做 log1p 作为 SingleR 输入
# expr_counts <- LayerData(combined, assay = "RNA", layer = "counts")
# # 基本检查：counts 必须有基因名和细胞名
# stopifnot(!is.null(rownames(expr_counts)))
# stopifnot(!is.null(colnames(expr_counts)))
# # log1p(counts) 得到“近似 log-normalized”，SingleR 可用
# expr_matrix <- log1p(expr_counts)
# # 节省内存：确保是稀疏矩阵
# if (!inherits(expr_matrix, "dgCMatrix")) expr_matrix <- as(expr_matrix, "dgCMatrix")
# # 4) 基础一致性检查：行名（基因名）必须存在且与参考库有交集
# if (is.null(rownames(expr_matrix)) || anyNA(rownames(expr_matrix))) {
#   stop("expr_matrix 没有有效的基因行名（rownames）。SingleR需要基因名来匹配参考库。")
# }
# common_genes <- base::intersect(rownames(expr_matrix), rownames(ref_data))
# cat("SingleR 基因交集数量:", length(common_genes), "\n")
# if (length(common_genes) < 500) {
#   warning("与参考库交集基因过少（<500）。可能是基因命名体系不一致（例如Ensembl vs Symbol），注释可靠性会显著下降。")
# }
# # 直接对单细胞进行注释
# singler_results <- SingleR(test = expr_matrix, ref = ref_data,
#                            labels = ref_data$label.main, de.method = "wilcox")
# # 6) 关键：使用pruned.labels（低置信度会变成NA），比直接labels更“发表级”
# #    同时保留原始labels便于排查
# combined$singler_label_raw <- singler_results$labels
# combined$singler_label_pruned <- singler_results$pruned.labels
# # 7) 结果数量检查：必须与细胞数匹配
# cat("SingleR结果数量(raw):", length(combined$singler_label_raw), "\n")
# cat("SingleR结果数量(pruned):", length(combined$singler_label_pruned), "\n")
# cat("细胞数量:", ncol(combined), "\n")
# if (length(combined$singler_label_raw) != ncol(combined)) stop("SingleR输出长度与细胞数不一致，请检查expr_matrix维度是否为“基因 x 细胞”。")
# # 8) cluster级别注释：对每个cluster做“多数投票”，并计算纯度
# cluster_labels <- data.frame(cluster = clusters, cell_type = NA_character_, purity = NA_real_, stringsAsFactors = FALSE)
#
# for (cl in clusters) {
#   cells_in_cluster <- WhichCells(combined, idents = cl)
#   labs <- combined$singler_label_pruned[cells_in_cluster]
#   labs <- labs[!is.na(labs)]  # 去掉低置信度的NA
#
#   if (length(labs) == 0) {
#     # 该cluster大多低置信度：标为Unassigned
#     cluster_labels[cluster_labels$cluster == cl, c("cell_type","purity")] <- list("Unassigned", 0)
#     next
#   }
#
#   tab <- sort(table(labs), decreasing = TRUE)
#   main_label <- names(tab)[1]
#   purity <- as.numeric(tab[1]) / sum(tab)
#
#   # 纯度阈值：<0.6 说明cluster内部标签混杂，建议不直接定性
#   if (purity < 0.60) {
#     cluster_labels[cluster_labels$cluster == cl, c("cell_type","purity")] <- list("Unassigned", purity)
#   } else {
#     cluster_labels[cluster_labels$cluster == cl, c("cell_type","purity")] <- list(main_label, purity)
#   }
# }
# # 9) 映射回Seurat对象：不依赖plyr，直接用命名向量索引（更稳）
# cluster_map <- setNames(cluster_labels$cell_type, cluster_labels$cluster)
# combined$cell_type_singler_cluster <- base::unname(cluster_map[as.character(combined$seurat_clusters)])
#
# # 10) 输出对照表（用于发表/补充材料：cluster vs label）
# comparison_table <- table(combined$seurat_clusters, combined$singler_label_pruned)
#
# # 11) 画图：cluster级别（推荐用于主图） + 单细胞pruned（用于补图排查）
# p_cluster <- DimPlot(combined, reduction = "umap", group.by = "cell_type_singler_cluster", label = TRUE) +
#   ggtitle("SingleR 注释（Cluster-level, pruned + majority vote）")
# ggsave("./temp/UMAP_celltype_singler_cluster_level.pdf", plot = p_cluster, width = 12, height = 8)
#
# p_cell <- DimPlot(combined, reduction = "umap", group.by = "singler_label_pruned", label = FALSE) +
#   ggtitle("SingleR 注释（Cell-level, pruned.labels）") +
#   theme(legend.position = "right")
# ggsave("./temp/UMAP_celltype_singler_cell_level_pruned.pdf", plot = p_cell, width = 12, height = 8)
# cat("7. 完成【使用SingleR进行细胞类型注释】\n")
# # 12) 保存结果（对象 + cluster注释表 + 对照矩阵）
# saveRDS(cluster_labels, file = "./temp/singler_cluster_labels_res7.rds")
# saveRDS(combined, file = "./temp/singleR_cell_annotation_7.rds")
#

# # --------------------- 8. 提取T cells和NKT cell并进行再分析 ---------------------
# cat("8. 开始【提取T细胞和NK细胞进行深入分析】\n") # 如果标签不同，请调整以下代码中的名称
# target_cell_types <- c("T cells", "NKT")
#
# # 确认目标细胞类型在数据中存在
# cat("数据中的细胞类型标签:", unique(combined$cell_type_singler), "\n")
# # 提取T细胞和NK细胞
# t_nk_cells <- subset(combined, subset = cell_type_singler %in% target_cell_types)
# cat("8. 完成【提取T细胞和NK细胞进行深入分析】。提取的细胞数量:", ncol(t_nk_cells), "\n")
# saveRDS(t_nk_cells, file = "./temp/t_nkt_cells_subset_8.rds")
# t_nk_cells <- readRDS("./temp/t_nkt_cells_subset_8.rds")
#
#
# # --------------------- 9. 对T细胞和NK细胞子集重新进行标准化和降维 ---------------------
# cat("9. 开始【对T细胞和NK细胞子集重新进行标准化和降维】\n")
# # 重置Idents以便后续分析
# Idents(t_nk_cells) <- "cell_type_singler"
# # 对子集重新进行标准化
# t_nk_cells <- NormalizeData(t_nk_cells)
# t_nk_cells <- FindVariableFeatures(t_nk_cells, selection.method = "vst", nfeatures = 2000)
# t_nk_cells <- ScaleData(t_nk_cells)
# # PCA降维
# t_nk_cells <- RunPCA(t_nk_cells, npcs = 30)
# # 确定应使用的PC数量
# # 绘制弯道图(Elbow plot)
# elbow_plot <- ElbowPlot(t_nk_cells, ndims = 30)
# ggsave("./temp/t_nk_elbow_plot.pdf", plot = elbow_plot, width = 8, height = 6)
# # 选择前20个PC进行后续分析（可以根据弯道图调整）
# pc_dims <- 1:20
# # 细胞聚类
# t_nk_cells <- FindNeighbors(t_nk_cells, dims = pc_dims)
# t_nk_cells <- FindClusters(t_nk_cells, resolution = 0.6) # 可以调整resolution以获得不同粒度的聚类
# # UMAP降维可视化
# t_nk_cells <- RunUMAP(t_nk_cells, dims = pc_dims)
# # 可视化聚类结果
# p1 <- DimPlot(t_nk_cells, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
#   ggtitle("T和NK细胞子集聚类")
# ggsave("./temp/t_nk_clustering.pdf", plot = p1, width = 10, height = 8)
# # 按原始细胞类型标签展示
# p2 <- DimPlot(t_nk_cells, reduction = "umap", label = TRUE, group.by = "cell_type_singler") +
#   ggtitle("T和NK细胞子集原始标签")
# ggsave("./temp/t_nk_original_labels.pdf", plot = p2, width = 10, height = 8)
# cat("9. 完成【对T细胞和NK细胞子集重新进行标准化和降维】\n")
#
#
# # --------------------- 10. 对T和NK细胞子集进行SingleR细分类型注释 ---------------------
# cat("10. 开始【对T和NK细胞子集进行SingleR细分类型注释】\n")
# # 使用更精细的参考数据集，例如ImmGen数据集的细分类型
# ref_data <- celldex::ImmGenData()
# # 提取表达矩阵用于SingleR注释
# t_nk_expr_matrix <- LayerData(t_nk_cells, assay = "RNA", layer = "data")
# # 运行优化版SingleR
# t_nk_singler_results <- SingleR(
#   test = t_nk_expr_matrix, ref = ref_data,
#   labels = ref_data$label.fine, de.method = "wilcox")
# # 添加细分类型标签到对象
# t_nk_cells$singler_fine_labels <- t_nk_singler_results$labels
# # 为每个新聚类确定主要细胞类型
# new_clusters <- unique(t_nk_cells$seurat_clusters)
# new_cluster_labels <- data.frame()
# for (cl in new_clusters) {
#   cells_in_cluster <- WhichCells(t_nk_cells, idents = cl)
#   # 获取该cluster中所有细胞的标签
#   labels_in_cluster <- t_nk_cells$singler_fine_labels[colnames(t_nk_cells) %in% cells_in_cluster]
#   # 找出最常见的标签
#   label_counts <- table(labels_in_cluster)
#   main_label <- names(label_counts)[which.max(label_counts)]
#
#   new_cluster_labels <- rbind(new_cluster_labels, data.frame(
#     cluster = cl,
#     cell_type = main_label,
#     proportion = max(label_counts) / sum(label_counts)
#   ))
# }
# # 映射细分类型标签到聚类
# new_cluster_map <- setNames(new_cluster_labels$cell_type, new_cluster_labels$cluster)
# t_nk_cells$fine_cell_type <- plyr::mapvalues(
#   x = t_nk_cells$seurat_clusters,
#   from = names(new_cluster_map),
#   to = as.character(new_cluster_map)
# )
# # 可视化细分类型
# p3 <- DimPlot(t_nk_cells, reduction = "umap", group.by = "fine_cell_type", label = TRUE) +
#   ggtitle("T和NK细胞细分类型") +
#   theme(legend.position = "right")
# ggsave("./temp/t_nk_fine_cell_types.pdf", plot = p3, width = 12, height = 8)
# cat("10. 完成【对T和NK细胞子集进行SingleR细分类型注释】\n")
# saveRDS(t_nk_cells, file = "./temp/t_nkt_cells_singleR_10.rds")
# t_nk_cells <- readRDS("./temp/t_nkt_cells_singleR_10.rds")
#
#
# # --------------------- 11. 可视化指定基因在T和NK细胞子集中的表达 ---------------------
# cat("11. 可视化指定基因在T和NK细胞子集中的表达\n")
# # 设置目标基因列表
# target_genes <- c('CCR7', 'SELL', 'PDCD1', 'CTLA4', 'ENTPD1', 'HAVCR2', 'LAG3', 'LAYN', 'SPRY1', 'IL7R',
#                   'LYAR', 'KLRG1', 'GZMK', 'EOMES', 'KLF2', 'TCF7', 'ZNF683', 'CXCR6', 'ITGAE', 'ITGA1',
#                   'FGFBP2', 'CX3CR1', 'TBX21', 'IFNG', 'GZMA', 'GZMB', 'GNLY', 'PRF1', 'IFIT1', 'IFIT3',
#                   'IFI44L', 'IFI6', 'STMN1', 'PCLAF', 'PCNA', 'MKI67', 'CCL4L2', 'XCL1', 'XCL2', 'CCL4',
#                   'PRDM1', 'JUNB', 'FOSB')
# # 大小写不敏感地查找基因
# gene_results <- find_genes_case_insensitive(target_genes, rownames(t_nk_cells))
# genes_exist <- gene_results$matched
# genes_missing <- gene_results$missing
# # 打印匹配情况
# cat("找到的基因:", length(genes_exist), "个\n")
# cat("匹配的基因对应表:\n")
# for (i in 1:length(gene_results$name_map)) {
#   original <- names(gene_results$name_map)[i]
#   matched <- gene_results$name_map[[i]]
#   if (original != matched) {
#     cat("目标基因:", original, "-> 数据集中的基因:", matched, "\n")
#   }
# }
# if (length(genes_missing) > 0) {
#   cat("以下基因不在数据集中:", paste(genes_missing, collapse = ", "), "\n")
# }
# # 为存在的基因生成表达图
# if (length(genes_exist) > 0) {
#   cat("为以下基因生成表达图:", paste(genes_exist, collapse = ", "), "\n")
#   # 在UMAP上显示基因表达
#   p4 <- FeaturePlot(t_nk_cells, features = genes_exist, ncol = 4, reduction = "umap")
#   ggsave("./temp/t_nk_genes_umap.pdf", plot = p4, width = 16, height = 4 * ceiling(length(genes_exist) / 4))
#   # 设置当前识别为细分类型，为热图绘制准备
#   Idents(t_nk_cells) <- "fine_cell_type"
#   # 绘制热图 - 不同细分类型中目标基因的平均表达
#   p5 <- DoHeatmap(t_nk_cells, features = genes_exist, group.by = "fine_cell_type") +
#     ggtitle("目标基因在T和NK细胞亚型中的表达热图")
#   ggsave("./temp/t_nk_genes_heatmap_by_fine_celltype.pdf", plot = p5, width = 12, height = max(10, length(genes_exist) / 3))
#   # 按聚类展示目标基因的表达热图
#   Idents(t_nk_cells) <- "seurat_clusters"
#   p6 <- DoHeatmap(t_nk_cells, features = genes_exist, group.by = "seurat_clusters") +
#     ggtitle("目标基因在不同聚类中的表达热图")
#   ggsave("./temp/t_nk_genes_heatmap_by_cluster.pdf", plot = p6, width = 12, height = max(10, length(genes_exist) / 3))
#   # 小提琴图 - 按细分类型展示基因表达
#   Idents(t_nk_cells) <- "fine_cell_type"
#   # 如果基因数量较多，分组绘制小提琴图
#   max_genes_per_plot <- 4
#   for (i in 1:ceiling(length(genes_exist) / max_genes_per_plot)) {
#     start_idx <- (i - 1) * max_genes_per_plot + 1
#     end_idx <- min(i * max_genes_per_plot, length(genes_exist))
#     gene_subset <- genes_exist[start_idx:end_idx]
#
#     p7 <- VlnPlot(t_nk_cells, features = gene_subset,
#                   ncol = length(gene_subset), pt.size = 0) # 不显示点
#     ggsave(paste0("./temp/t_nk_genes_violin_", i, ".pdf"),
#            plot = p7, width = 4 * length(gene_subset), height = 6) }
#   # DotPlot展示基因表达与细胞类型的关系
#   p8 <- DotPlot(t_nk_cells, features = genes_exist, dot.scale = 8) +
#     ggtitle("目标基因在不同细胞类型中的表达") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   ggsave("./temp/t_nk_genes_dotplot.pdf", plot = p8,
#          width = max(12, length(genes_exist) * 0.4), height = 8)
#   # 添加聚类标记基因分析，帮助理解不同簇的特征
#   Idents(t_nk_cells) <- "seurat_clusters"
#   markers <- FindAllMarkers(t_nk_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#   # 保存所有标记基因
#   write.csv(markers, "./temp/t_nk_cluster_markers.csv", row.names = FALSE)
#   # 提取目标基因在不同簇中的表达统计
#   target_gene_stats <- markers[markers$gene %in% genes_exist,]
#   if (nrow(target_gene_stats) > 0) {
#     write.csv(target_gene_stats, "./temp/t_nk_target_gene_markers.csv", row.names = FALSE)
#     # 绘制目标基因中显著标记基因的热图
#     if (nrow(target_gene_stats) >= 3) {
#       top_markers <- target_gene_stats %>%
#         group_by(cluster) %>%
#         top_n(3, avg_log2FC)
#       p9 <- DoHeatmap(t_nk_cells, features = unique(top_markers$gene)) +
#         ggtitle("目标基因中的显著标记基因热图")
#       ggsave("./temp/t_nk_target_markers_heatmap.pdf", plot = p9, width = 12, height = 10)
#     }
#   }
#   if (require(tidyverse)) {
#     # 计算每个基因在每个聚类中的表达比例
#     gene_prop_df <- data.frame()
#     for (cl in levels(Idents(t_nk_cells))) {
#       cl_cells <- WhichCells(t_nk_cells, idents = cl)
#       # 使用Seurat V5兼容的方法获取表达矩阵
#       # 首先获取原始数据
#       cl_data <- FetchData(t_nk_cells, vars = genes_exist, cells = cl_cells)
#       # 计算每个基因的阳性比例
#       gene_props <- sapply(cl_data, function(x) mean(x > 0))
#       temp_df <- data.frame(
#         gene = names(gene_props),
#         proportion = as.numeric(gene_props),
#         cluster = cl
#       )
#       gene_prop_df <- rbind(gene_prop_df, temp_df)
#     }
#     # 绘制比例热图
#     p10 <- ggplot(gene_prop_df, aes(x = gene, y = cluster, fill = proportion)) +
#       geom_tile() +
#       scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(0, 1), name = "表达比例") +
#       theme_minimal() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       ggtitle("目标基因在各聚类中的表达比例") +
#       xlab("基因") +
#       ylab("聚类")
#     ggsave("./temp/t_nk_gene_proportion_heatmap.pdf", plot = p10,
#            width = max(10, length(genes_exist) * 0.4), height = 8)
#   }
# }
# # 添加类似示例图的Z-score标准化表达热图，包含层次聚类树状图
# if (length(genes_exist) > 0) {
#   cat("生成分组展示的Z-score标准化表达热图，包含层次聚类树状图\n")
#   # 获取基因表达数据
#   gene_data <- FetchData(t_nk_cells, vars = c(genes_exist, "seurat_clusters", "orig.ident"))
#   # 确保orig.ident包含你的三个实验组
#   if (!"orig.ident" %in% colnames(gene_data)) {
#     cat("警告: 数据中没有实验组信息 (orig.ident)，请确认实验组的列名\n")
#   } else {
#     # 更新为您的三个组别: "control1", "s1_hIL", "s40kD"
#     group_order <- c("control1", "s1_hIL", "s40kD")
#     # 按实验组和聚类计算每个基因的平均表达
#     avg_exp <- gene_data %>%
#       group_by(orig.ident, seurat_clusters) %>%
#       summarise(across(all_of(genes_exist), mean), .groups = "drop")
#     # 转换为长格式以便于计算Z-score
#     avg_exp_long <- avg_exp %>%
#       pivot_longer(cols = all_of(genes_exist),
#                    names_to = "gene",
#                    values_to = "expression")
#     # 按基因计算Z-score
#     avg_exp_long <- avg_exp_long %>%
#       group_by(gene) %>%
#       mutate(zscore = scale(expression)[, 1]) %>%
#       ungroup()
#     # 确保使用正确的组顺序
#     if (all(group_order %in% unique(avg_exp_long$orig.ident))) {
#       avg_exp_long$orig.ident <- factor(avg_exp_long$orig.ident, levels = group_order)
#     } else {
#       cat("警告: 指定的组名与数据不匹配，使用数据中的原始组名\n")
#     }
#     # 对基因进行聚类以获得合理的排序
#     # 首先创建宽格式的表达矩阵用于聚类
#     gene_matrix <- avg_exp_long %>%
#       pivot_wider(id_cols = "gene",
#                   names_from = c("orig.ident", "seurat_clusters"),
#                   values_from = "zscore") %>%
#       as.data.frame()
#     row_names <- gene_matrix$gene
#     gene_matrix <- gene_matrix[, -1]
#     rownames(gene_matrix) <- row_names
#     # 进行层次聚类
#     gene_hclust <- hclust(dist(gene_matrix))
#     gene_order <- rownames(gene_matrix)[gene_hclust$order]
#
#     # 应用聚类顺序到数据
#     avg_exp_long$gene <- factor(avg_exp_long$gene, levels = gene_order)
#     # 同样，对聚类进行聚类（如果需要的话）
#     # 创建一个按聚类整理的矩阵
#     cluster_matrix <- avg_exp_long %>%
#       group_by(orig.ident, seurat_clusters, gene) %>%
#       summarise(zscore = mean(zscore), .groups = "drop") %>%
#       pivot_wider(id_cols = c("orig.ident", "seurat_clusters"),
#                   names_from = "gene", values_from = "zscore") %>%
#       as.data.frame()
#     # 为每个实验组创建聚类树状图
#     plots <- list()
#     dendro_plots <- list()
#     # 获取所有组别
#     all_groups <- unique(avg_exp_long$orig.ident)
#     # 获取所有聚类
#     all_clusters <- unique(avg_exp_long$seurat_clusters)
#     # 创建用于存储基因聚类树的数据结构
#
#     # 创建每个实验组的热图和树状图
#     for (group in all_groups) {
#       # 筛选当前组的数据
#       conflicts_prefer(dplyr::filter)
#       group_data <- avg_exp_long %>% filter(orig.ident == group)
#       # 为当前组的数据重新计算层次聚类（可选，如果想为每个组单独聚类）
#       group_gene_matrix <- group_data %>%
#         pivot_wider(id_cols = "gene", names_from = "seurat_clusters", values_from = "zscore") %>%
#         as.data.frame()
#       gene_names <- group_gene_matrix$gene
#       group_gene_matrix <- group_gene_matrix[, -1]
#       rownames(group_gene_matrix) <- gene_names
#       # 进行层次聚类
#       group_gene_hclust <- hclust(dist(group_gene_matrix))
#       # 转换聚类为可绘制的数据
#       gene_ddata <- dendro_data(group_gene_hclust, type = "rectangle")
#       # 调整树状图标签以匹配热图的基因顺序
#       # 创建一个映射，将聚类顺序索引转换为因子顺序
#       order_map <- setNames(1:length(gene_names), gene_names[group_gene_hclust$order])
#       # 创建树状图
#       dendro_plot <- ggplot() +
#         geom_segment(data = segment(gene_ddata),
#                      aes(x = x, y = y, xend = xend, yend = yend)) +
#         theme_minimal() +
#         theme(
#           axis.text = element_blank(),
#           axis.title = element_blank(),
#           panel.grid = element_blank(),
#           plot.margin = margin(0, 0, 0, 0)
#         ) +
#         scale_y_reverse() +  # 反转y轴以使树状图正确显示
#         coord_flip()  # 翻转坐标以与热图匹配
#       dendro_plots[[group]] <- dendro_plot
#       # 热图标题，更新为您的组名
#       title_text <- paste0("NK/T cells\n", group)
#       # 创建热图，保持与现有代码一致
#       p <- ggplot(group_data, aes(x = seurat_clusters, y = gene, fill = zscore)) +
#         geom_tile() +
#         scale_fill_gradientn(
#           colors = colorRampPalette(c("#00008B", "#4169E1", "white", "#FF6347", "#8B0000"))(100),
#           limits = c(-2, 2), name = "Z-score"
#         ) +
#         theme_minimal() +
#         labs(title = title_text, x = "", y = "") +
#         theme(
#           axis.text.y = element_text(size = 8),
#           axis.text.x = element_text(size = 8),
#           plot.title = element_text(size = 10, hjust = 0.5),
#           legend.position = "none",  # 最终图只在右侧显示一个图例
#           panel.grid = element_blank(),
#           panel.border = element_rect(fill = NA, color = "black", size = 0.5)
#         )
#       plots[[group]] <- p
#     }
#     # 使用patchwork组合多个热图和树状图
#     if (require(patchwork)) {
#       # 对于每个组，创建树状图+热图的组合
#       combined_group_plots <- list()
#       for (group in all_groups) {
#         # 设置树状图的宽度比例（通常树状图宽度约为热图的1/4或1/5）
#         combined_group_plots[[group]] <- dendro_plots[[group]] +
#           plots[[group]] +
#           plot_layout(widths = c(1, 5))
#       }
#       # 然后将所有组的组合图水平排列
#       combined_plot <- combined_group_plots[[1]] +
#         combined_group_plots[[2]] +
#         combined_group_plots[[3]] +
#         plot_layout(guides = "collect") &
#         theme(legend.position = "right")
#       # 保存组合热图
#       ggsave("./temp/combined_zscore_heatmap_with_dendro.pdf",
#              plot = combined_plot, width = 15,  # 增加宽度以适应树状图
#              height = max(6, length(genes_exist) * 0.3))
#       # 也可以尝试使用ComplexHeatmap包创建更专业的热图（推荐）
#       if (require(ComplexHeatmap)) {
#         cat("使用ComplexHeatmap创建更专业的热图，包含完整的树状图...\n")
#         # 为每个组创建一个ComplexHeatmap对象
#         heatmap_list <- list()
#         for (i in seq_along(all_groups)) {
#           group <- all_groups[i]
#           # 提取当前组的数据并转换为矩阵
#           group_data <- avg_exp_long %>%
#             filter(orig.ident == group) %>%
#             pivot_wider(id_cols = "gene", names_from = "seurat_clusters", values_from = "zscore") %>%
#             as.data.frame()
#           row_names <- group_data$gene
#           group_data <- group_data[, -1]
#           rownames(group_data) <- row_names
#           # 转换为矩阵
#           mat <- as.matrix(group_data)
#           # 创建聚类对象
#           row_dend <- hclust(dist(mat))
#           col_dend <- hclust(dist(t(mat)))
#           # 设置标题，更新为您的组名
#           title_text <- paste0("NK/T cells - ", group)
#           # 设置聚类组颜色
#           cluster_colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(colnames(mat)))
#           names(cluster_colors) <- colnames(mat)
#           column_anno <- HeatmapAnnotation(Cluster = colnames(mat), col = list(Cluster = cluster_colors), show_legend = TRUE)
#           # 创建热图
#           ht <- Heatmap(
#             mat,
#             name = paste0("Z-score_", group),
#             cluster_rows = row_dend,
#             cluster_columns = col_dend,
#             show_row_dend = TRUE,
#             show_column_dend = TRUE,
#             column_title = title_text,
#             row_names_gp = gpar(fontsize = 8),
#             column_names_gp = gpar(fontsize = 8),
#             heatmap_legend_param = list(title = "Z-score"),
#             col = colorRamp2(c(-2, -1, 0, 1, 2), c("#00008B", "#4169E1", "white", "#FF6347", "#8B0000")),
#             top_annotation = column_anno
#           )
#           heatmap_list[[i]] <- ht
#         }
#         # 组合所有热图
#         combined_heatmap <- NULL
#         for (i in seq_along(heatmap_list)) {
#           if (is.null(combined_heatmap)) {
#             combined_heatmap <- heatmap_list[[i]]
#           } else {
#             combined_heatmap <- combined_heatmap + heatmap_list[[i]]
#           }
#         }
#         # 保存ComplexHeatmap版本
#         pdf("./temp/complex_heatmap_with_dendro.pdf", width = 15, height = max(8, length(genes_exist) * 0.3))
#         draw(combined_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
#         dev.off()
#       }
#     }
#   }
# }
# cat("T和NK细胞深入分析完成，结果已保存\n")
# saveRDS(t_nk_cells, file = "./temp/t_nk_cells_visulization_11.rds")
# t_nk_cells <- readRDS("./temp/t_nk_cells_visulization_11.rds")
#
#
# # --------------------- 12. GO/KEGG 富集分析 & 柱状图 ---------------------
# cat("12. 开始 GO 和 KEGG 富集分析\n")
# # 12.1 准备基因列表
# all_markers <- readRDS("./temp/marker_analysis_6.rds")
# markers_dt <- as.data.table(all_markers, keep.rownames = 'gene')
# # 将 markers_dt 从宽表 melt 成长表，并提取所有基因
# long_dt <- melt(markers_dt, measure = patterns(gene = "\\.gene$"), value.name = "gene")
# all_genes <- unique(long_dt$gene)
#
# # 12.2 SYMBOL 转 ENTREZID
# entrez_map <- mapIds(org.Mm.eg.db, keys = all_genes,
#                      column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
# entrez_map <- na.omit(entrez_map)   # 去掉无法匹配的
#
# # 12.3 运行 GO 富集 (BP)
# df_bp <- as.data.frame(
#   enrichGO(gene = entrez_map,
#                    OrgDb = org.Mm.eg.db,
#                    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
#                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
# )
# df_bp$Category <- "BP"
#
# df_cc <- as.data.frame( # 用 ont="CC" 的 enrichGO 对象
#   enrichGO(gene = entrez_map,
#            OrgDb = org.Mm.eg.db,
#            keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH",
#            pvalueCutoff = 0.05, qvalueCutoff = 0.2)
# )
# df_cc$Category <- "CC"
#
# df_mf <- as.data.frame( # 用 ont="MF" 的 enrichGO 对象
#   enrichGO(gene = entrez_map,
#            OrgDb = org.Mm.eg.db,
#            keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH",
#            pvalueCutoff = 0.05, qvalueCutoff = 0.2)
# )
# df_mf$Category <- "MF"
#
# # 2. 合并三者，并选 top N 条（这里示例取每类前 5 项）
#
# topn <- 10
# df_top <- bind_rows(df_bp, df_cc, df_mf) %>%
#   group_by(Category) %>%
#   slice_min(order_by = p.adjust, n = topn) %>%
#   ungroup() %>%
#   # 计算 -log10(p.adjust) 做为横轴
#   mutate(logP = -log10(p.adjust)) %>%
#   # 为了让条按高度排好序，我们要为 Term 建一个因子，按 logP 降序
#   arrange(logP) %>%
#   mutate(Term = factor(Description, levels = unique(Description)))
#
# # 3. 绘图
# p <- ggplot(df_top, aes(x = logP, y = Term, fill = Category)) +
#   geom_col(width = 0.7) +
#   scale_fill_manual(
#     values = c(BP = "#2E8B57", CC = "#FFA500", MF = "#1E90FF")
#   ) +
#   labs(
#     x = expression(-log[10]("adjusted P-value")),
#     y = NULL,
#     title = "Marker 基因 GO 富集（BP/CC/MF）"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "right",
#         axis.text.y = element_text(size = 10),
#         panel.grid.major.y = element_blank())
#
# print(p)
# ggsave("./temp/go_combined_BP_CC_MF_barplot.pdf", p, width = 8, height = 6)
#
# # 12.4 运行 KEGG 富集
# ekegg <- enrichKEGG(gene = as.character(entrez_map),
#                     organism = "mmu", pAdjustMethod = "BH",
#                     pvalueCutoff = 0.05, qvalueCutoff = 0.2)
# # 12.6 可视化 KEGG
# p3 <- barplot(ekegg, showCategory = 10) +
#   ggtitle("Marker 基因 KEGG 通路富集") +
#   theme(axis.text.y = element_text(size = 10))
# p4 <- dotplot(ekegg, showCategory = 10) +
#   ggtitle("Marker 基因 KEGG 气泡图") +
#   theme(axis.text.y = element_text(size = 10))
# # 12.7 打印并保存所有图
# print(p3); ggsave("./temp/kegg_bar_marker.pdf", p3, width = 8, height = 6)
# print(p4); ggsave("./temp/kegg_dot_marker.pdf", p4, width = 8, height = 6)
#
# cat("12. GO 和 KEGG 富集分析完成，图表已保存到 ./temp/ 目录\n")
#
