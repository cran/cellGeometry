## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)

## ----eval=FALSE---------------------------------------------------------------
# library(zellkonverter)
# library(SingleCellExperiment)
# library(cellGeometry)

## ----eval=FALSE---------------------------------------------------------------
# # 33 Gb
# # 2,480,956 cells
# brain <- readH5AD("../c2f66cd5-4ff4-4578-876c-55783a57cf8f.h5ad",
#                   use_hdf5 = TRUE, reader = "R")
# 
# mat <- brain@assays@data$X
# rownames(mat) <- rownames(brain)  # need to add rownames (genes)
# meta <- brain@colData@listData

## ----eval=FALSE---------------------------------------------------------------
# # shows possible cell cluster subclasses and groups
# sort(table(meta$roi))
# sort(table(meta$supercluster_term))
# 
# mk <- cellMarkers(mat, subclass = meta$roi,
#                   cellgroup = meta$supercluster_term,
#                   dual_mean = TRUE, cores = 8)
# # 27 mins (intel 8 cores)

## ----eval=FALSE---------------------------------------------------------------
# library(AnnotationHub)
# ah <- AnnotationHub()
# ensDb_v110 <- ah[["AH113665"]]
# 
# mk <- gene2symbol(mk, ensDb_v110)

## ----eval=FALSE---------------------------------------------------------------
# signature_heatmap(mk, text = FALSE, show_row_names = FALSE,
#                   row_title_rot = 0, column_title_rot = 45)

## ----eval=FALSE---------------------------------------------------------------
# # non-neuronal cells
# # 888,263 cells
# # 4 Gb
# brainNN <- readH5AD("../99f27be8-9fac-451e-9723-9e4c7191589e.h5ad",
#                     use_hdf5 = TRUE, reader = "R")
# 
# mat2 <- brainNN@assays@data$X
# rownames(mat2) <- rownames(brainNN)  # add rownames (genes)
# meta2 <- brainNN@colData@listData
# 
# sort(table(meta2$supercluster_term))
# sort(table(meta2$cell_type))
# 
# mkNN <- cellMarkers(mat2, subclass = meta2$cell_type,
#                     cellgroup = meta2$supercluster_term,
#                     dual_mean = TRUE, cores = 8)
# # 9 mins (intel 8 cores)
# 
# mkNN <- gene2symbol(mkNN, ensDb_v110)
# 
# signature_heatmap(mkNN)

## ----eval=FALSE---------------------------------------------------------------
# mkm <- mergeMarkers(mk, mkNN, transform = "none")
# mkm <- updateMarkers(mkm, expfilter = 0.2)

## ----eval=FALSE---------------------------------------------------------------
# signature_heatmap(mkm, top = 5,
#                   text = FALSE, show_row_names = FALSE,
#                   row_title_rot = 0, column_title_rot = 75,
#                   scale = "sphere")

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("brain_sig.png")

## ----eval=FALSE---------------------------------------------------------------
# # view signature for just the Amygdata excitatory group
# signature_heatmap(mkm, subset = "Amygdala excitatory")

## ----out.width='50%', fig.align = 'center', echo=FALSE------------------------
knitr::include_graphics("brain_sig_amygdala.png")

## ----eval=FALSE---------------------------------------------------------------
# # generate neuronal cell count samples
# set.seed(3)
# sim_counts <- generate_samples(mk, 30)
# sim_percent <- sim_counts / rowSums(sim_counts) * 100
# 
# # simulate neuronal bulk
# # sample from neuronal count matrix
# # 35 mins (Intel)
# sim_sampled <- simulate_bulk(mat, sim_counts, meta$roi,
#                              times = 1)
# 
# rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)

## ----eval=FALSE---------------------------------------------------------------
# fit <- deconvolute(mk, sim_sampled,
#                    arith_mean = TRUE,
#                    use_filter = FALSE)
# mset <- metric_set(sim_percent, fit$subclass$percent)
# summary(mset)
# 
# pdf("../brain_sim_neuron.pdf",
#     width = 12, height = 12.5)
# plot_set(sim_counts, fit$subclass$output, show_zero = T,
#          mfrow = c(8, 8))
# plot_set(sim_percent, fit$subclass$percent, show_zero = T,
#          mfrow = c(8, 8))
# dev.off()

## ----eval=FALSE---------------------------------------------------------------
# # generate non-neuronal cell count samples
# set.seed(3)
# sim_countsNN <- generate_samples(mkNN, 30)
# sim_percentNN <- sim_countsNN / rowSums(sim_countsNN) * 100
# 
# # simulate non-neuronal bulk
# # sample from non-neuronal count matrix
# # 6.32 mins (Intel)
# sim_sampledNN <- simulate_bulk(mat2, sim_countsNN, meta2$cell_type,
#                                times = 1)
# 
# rownames(sim_sampledNN) <- gene2symbol(rownames(sim_sampledNN), ensDb_v110)

## ----eval=FALSE---------------------------------------------------------------
# # check genenames are the same in both datasets
# identical(rownames(sim_sampled), rownames(sim_sampledNN))
# ## TRUE
# 
# # merge pseudobulk counts
# sim_sampled_merge <- sim_sampled + sim_sampledNN
# 
# # merge sample cell counts (ground truth)
# sim_counts_merge <- cbind(sim_counts, sim_countsNN)
# sim_percent_merge <- sim_counts_merge / rowSums(sim_counts_merge) * 100

## ----eval=FALSE---------------------------------------------------------------
# fitm <- deconvolute(mkm, sim_sampled_merge,
#                     arith_mean = TRUE, use_filter = FALSE, cores = 8)
# # 19.2 secs (ARM, 1 core)
# # 3.22 secs (ARM, 8 cores)

## ----eval=FALSE---------------------------------------------------------------
# mset <- metric_set(sim_percent_merge, fitm$subclass$percent)
# summary(mset)
# ##  pearson.rsq          Rsq               RMSE
# ## Min.   :0.6634   Min.   :-0.5747   Min.   :0.02479
# ## 1st Qu.:0.9298   1st Qu.: 0.9146   1st Qu.:0.05244
# ## Median :0.9789   Median : 0.9718   Median :0.07411
# ## Mean   :0.9452   Mean   : 0.8969   Mean   :0.10484
# ## 3rd Qu.:0.9917   3rd Qu.: 0.9855   3rd Qu.:0.11979
# ## Max.   :0.9999   Max.   : 0.9958   Max.   :0.51257

## ----eval=FALSE---------------------------------------------------------------
# # scatter plots
# pdf("../brain_sim_merge.pdf",
#     width = 12, height = 12.5)
# plot_set(sim_counts_merge, fitm$subclass$output, show_zero = TRUE,
#          mfrow = c(8, 8))
# plot_set(sim_percent_merge, fitm$subclass$percent, show_zero = TRUE,
#          mfrow = c(8, 8))
# dev.off()

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("brain_plotset.png")

