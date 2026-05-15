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
# --------------------- 1. 读取多个组数据并构建 Seurat 对象 ---------------------
cat("1. 开始【读取多个组别数据并构建Seurat对象】\n")
group_paths <- list(control1 = "./data/control1/",
                    s1_hIL = "./data/s1-hIL/",
                    s40kD = "./data/s40kD/")
seurat_list <- list()
for (grp in names(group_paths)) {
  path <- group_paths[[grp]]
  cat("正在加载组别：", grp, "，数据路径：", path, "\n")

  # 使用 Seurat 内置函数 ReadMtx 读取 matrix.mtx、features.tsv 和 barcodes.tsv
  counts <- ReadMtx(
    mtx = paste0(path, "matrix.mtx"),
    features = paste0(path, "features.tsv"),
    cells = paste0(path, "barcodes.tsv")
  )
  # 构建 Seurat 对象
  seu_obj <- CreateSeuratObject(counts = counts, project = grp, min.cells = 3, min.features = 0)
  seu_obj$group <- grp
  seurat_list[[grp]] <- seu_obj
}
# 合并各组对象
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "CombinedData")
cat("1. 完成【读取多个组别数据并构建Seurat对象】\n")
saveRDS(combined, file = "./temp/combined_raw_data_1.rds")
combined <- readRDS("./temp/combined_raw_data_1.rds")

# --------------------- 2. 数据过滤 ---------------------
cat("2. 开始【数据过滤】\n")
# 基因过滤：保留表达大于0的细胞数大于3的基因
counts_matrix <- LayerData(combined, assay = "RNA", layer = "counts")
genes_cell_counts <- rowSums(counts_matrix > 0)
genes_to_keep <- names(genes_cell_counts[genes_cell_counts > 3])
cat("基因过滤前：", nrow(combined), "个基因\n")
combined <- subset(combined, features = genes_to_keep)
cat("基因过滤后：", nrow(combined), "个基因\n")
# 计算线粒体基因比例
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
# 血红蛋白基因
combined[["percent.hb"]] <- PercentageFeatureSet(combined, pattern = "^Hba-|^Hbb-")
# 核糖体基因
combined[["percent.rb"]] <- PercentageFeatureSet(combined, pattern = "^Rps|^Rp")
# # 风琴图，后续 QC 阈值提供参考
# p1= VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),layer = "counts",pt.size = 0.1,raster = FALSE, ncol = 5)
# p1_group= VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),layer = "counts",group.by = "group",pt.size = 0.1,raster = FALSE, ncol = 5)
# p_qc <- p1 / p1_group
# ggsave(filename = "temp/QC_violin_overall_and_group.pdf",plot = p_qc, width = 14, height = 10, dpi = 300)
# 后续 QC 阈值提供参考，决定阈值的关键证据，nCount_RNA vs percent.mt（坏细胞通常：低 nCount + 高 mt）；nCount_RNA vs nFeature_RNA（doublet 通常：高 nCount + 高 nFeature）
df_qc <- combined@meta.data
p2 <- ggplot(df_qc, aes(x = nCount_RNA, y = percent.mt)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + labs(x = "nCount_RNA (UMI)", y = "percent.mt") + theme_bw()
p3 <- ggplot(df_qc, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + labs(x = "nCount_RNA (UMI)", y = "nFeature_RNA (genes)") + theme_bw()
p_qc2 <- p2 / p3
ggsave(filename = "temp/QC_scatter_preQC.pdf",plot = p_qc2, width = 14, height = 10, dpi = 300)
# ---- 2) Print pre-QC cell counts per group ----让你知道过滤前每组规模，便于判断过滤是否过度。
cat("==== Pre-QC cell counts ====\n")
print(table(combined$group))
cat("Total cells:", ncol(combined), "\n\n")
# 计算“percent.mt > 60% 的细胞占比”，快速感知你的数据里极端坏细胞多不多
extreme_mt <- mean(combined$percent.mt > 60) * 100
cat(sprintf("Diagnostic: percent.mt > 60%% cells = %.2f%%\n\n", extreme_mt))
# ---- 3) Define medium thresholds ----设置“中等阈值”
min_features <- 500
max_features <- 6500
min_counts   <- 1500   # 如果你发现过滤后细胞数太少，最常调的是 min_counts（1000→800）
max_counts   <- 30000
max_mt       <- 15     # 如果你发现过滤后细胞数太少，最常调的是 max_mt（30→35/40）。
# ---- 4) HB/RB extreme cutoff (per-group 99th percentile) ----hb/rb 不要设死值，而是“每组只砍掉最极端的 1%”
# This avoids one group dominating the global quantile.
md <- combined@meta.data %>%
  mutate(cell = rownames(.))
hb_cut_by_group <- md %>%  #每组情况可能不同，用“每组各自 top1%”更公平，不会出现某组被误砍更多。
  group_by(group) %>%
  summarise(hb_q99 = quantile(percent.hb, probs = 0.99, na.rm = TRUE), rb_q99 = quantile(percent.rb, probs = 0.99, na.rm = TRUE), .groups = "drop")
cat("==== Per-group 99% quantiles for HB/RB ====\n") #把上面每组阈值打印出来，让你知道 hb/rb “砍掉极端值”到底砍到哪面
print(hb_cut_by_group)
cat("\n")
# Join thresholds back to metadata 这样每个细胞都知道“自己所在组的阈值是多少”，方便下一步筛选
md2 <- md %>%
  left_join(hb_cut_by_group, by = "group")
cells_keep <- md2 %>%  #真正执行 QC 过滤：选出要保留的细胞3v
  dplyr::filter(nFeature_RNA >= min_features, nFeature_RNA <= max_features, nCount_RNA >= min_counts, nCount_RNA   <= max_counts, percent.mt   <= max_mt, percent.hb   <= hb_q99, percent.rb   <= rb_q99) %>%
  dplyr::pull(cell)  #把保留细胞的名字抽出来做成一个向量
cat("Cells kept after QC:", length(cells_keep), " / ", ncol(combined), "\n\n") #看每组保留率是否合理、是否某组被砍得特别严重
combined_qc <- subset(combined, cells = cells_keep)
# ---- 5) Print post-QC cell counts per group ----
cat("==== Post-QC cell counts ====\n")
print(table(combined_qc$group))
cat("Total cells:", ncol(combined_qc), "\n\n")
# ---- 6) Post-QC plots (violin + scatter) ----
# (A) Violin plots: set raster=FALSE to avoid rasterization-related internal errors
p_vln_all <- VlnPlot(combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "percent.rb"), pt.size = 0.1, raster = FALSE, ncol = 5)
p_vln_grp <- VlnPlot(combined_qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "percent.rb"), group.by = "group", pt.size = 0.1, raster = FALSE, ncol = 5)
p_vln <- p_vln_all / p_vln_grp
ggsave(filename = "temp/QC_violin_postQC_overall_and_group.pdf", plot = p_vln, width = 14, height = 10)
# (B) Scatter plots using ggplot facet (stable for large data)
df_qc <- combined_qc@meta.data
p_sc1 <- ggplot(df_qc, aes(x = nCount_RNA, y = percent.mt)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + scale_x_log10() + geom_hline(yintercept = max_mt, linetype = 2) + labs(x = "nCount_RNA (UMI, log10)", y = "percent.mt") + theme_bw()
p_sc2 <- ggplot(df_qc, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point(size = 0.2) + facet_wrap(~ group, ncol = 3, scales = "free") + scale_x_log10() + geom_hline(yintercept = c(min_features, max_features), linetype = 2) + geom_vline(xintercept = c(min_counts, max_counts), linetype = 2) + labs(x = "nCount_RNA (UMI, log10)", y = "nFeature_RNA (genes)") + theme_bw()
p_sc <- p_sc1/p_sc2
ggsave("temp/QC_scatter_postQC.pdf", plot = p_sc, width = 12, height = 6)
cat("Post-QC plots saved to:\n", file.path(getwd(), "temp/QC_violin_postQC_overall_and_group.pdf"), "\n", file.path(getwd(), "temp/QC_scatter_postQC_nCount_vs_percentMT.pdf"), "\n", file.path(getwd(), "temp/QC_scatter_postQC_nCount_vs_nFeature.pdf"), "\n")
# 细胞过滤：保留 nFeature_RNA 在 [200,7500] 内、线粒体比例小于 25%
# cat("细胞过滤前：", ncol(combined), "个细胞\n")
# combined <- subset(combined, subset = nFeature_RNA >= 200 &
#   nFeature_RNA <= 7500 &
#   percent.mt < 25)
# cat("细胞过滤后：", ncol(combined), "个细胞\n")
cat("2. 完成【数据过滤】\n")

## --------------------- 2. 数据过滤 ---------------------
#cat("2. 开始【数据过滤】\n")
## 基因过滤：保留表达大于0的细胞数大于3的基因
#counts_matrix <- LayerData(combined, assay = "RNA", layer = "counts")
#genes_cell_counts <- rowSums(counts_matrix > 0)
#genes_to_keep <- names(genes_cell_counts[genes_cell_counts > 3])
#cat("基因过滤前：", nrow(combined), "个基因\n")
#combined <- subset(combined, features = genes_to_keep)
#cat("基因过滤后：", nrow(combined), "个基因\n")
## 计算线粒体基因比例
#combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
## 细胞过滤：保留 nFeature_RNA 在 [200,7500] 内、线粒体比例小于 25%
#cat("细胞过滤前：", ncol(combined), "个细胞\n")
#combined <- subset(combined, subset = nFeature_RNA >= 200 &
#  nFeature_RNA <= 7500 &
#  percent.mt < 25)
#cat("细胞过滤后：", ncol(combined), "个细胞\n")
#cat("2. 完成【数据过滤】\n")


# --------------------- 3. 双细胞检测 ---------------------
cat("3. 开始【双细胞检测及过滤】\n")
# 预处理数据以便进行双细胞检测
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)
# 仅缩放变量基因以节省内存
combined <- ScaleData(combined, features = VariableFeatures(combined))
combined <- RunPCA(combined, npcs = 30)
# 使用 scDblFinder 进行双细胞检测
counts <- LayerData(combined, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays = list(counts = counts))
sce <- scDblFinder(sce, dbr = 0.06)  # 设置预期双细胞比例为6%
# 将双细胞信息添加回 Seurat 对象
combined$scDblFinder.class <- sce$scDblFinder.class
# 过滤掉被标记为双细胞的细胞
combined <- subset(combined, subset = scDblFinder.class == "singlet")
cat("3. 完成【双细胞检测及过滤】。剩余细胞数:", ncol(combined), "\n")


# --------------------- 4. 数据预处理 (归一化) ---------------------
cat("4. 开始【数据预处理】\n")
# 标准归一化和对数转换
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 1e4)
# 找出高变异基因
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
# 优化内存使用的缩放数据
# 1. 仅缩放变量基因而非所有基因
# 2. 分批处理缩放数据
var_genes <- VariableFeatures(combined)
cat("正在缩放数据（仅变量基因），节省内存使用...\n")
combined <- ScaleData(combined, features = var_genes, vars.to.regress = "percent.mt", clip.max = 10)
cat("4. 完成【数据预处理】\n")
saveRDS(combined, file = "./temp/combined_normalization_4.rds")
combined <- readRDS("./temp/combined_normalization_4.rds")


# --------------------- 5. 降维 & 聚类 ---------------------
cat("5. 开始【降维与聚类分析】\n")
# 使用PCA结果构建KNN图
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:15, k.param = 20)
# Louvain社区聚类，resolution控制聚类粒度
combined <- FindClusters(combined, resolution = 0.5)
# 计算UMAP（非线性降维可视化）
combined <- RunUMAP(combined, reduction = "pca", dims = 1:15)
# 绘制UMAP聚类结果
p <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP 聚类结果")
ggsave("./temp/UMAP_louvain.pdf", plot = p, width = 10, height = 8)
cat("5. 完成【降维与聚类分析】\n")
saveRDS(combined, file = "./temp/combined_umap_5.rds")
combined <- readRDS("./temp/combined_umap_5.rds")


