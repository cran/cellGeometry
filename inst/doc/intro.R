## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE
)

## ----eval = FALSE-------------------------------------------------------------
# # Bioconductor must be installed +/- updated first
# BiocManager::install(version = "3.xx")  # set to latest version
# 
# # minimum necessary Bioconductor packages to install cellGeometry package
# BiocManager::install(c("ensembldb", "DelayedArray"))
# 
# # packages needed to read h5ad files
# BiocManager::install(c("zellkonverter", "rhdf5", "HDF5Array"))
# 
# # optional, if you are using Seurat
# install.packages("Seurat")
# 
# # package needed to convert ensembl gene ids to symbols
# BiocManager::install("AnnotationHub")

## ----eval = FALSE-------------------------------------------------------------
# devtools::install_github("myles-lewis/cellGeometry")

## ----out.width='100%', echo=FALSE---------------------------------------------
knitr::include_graphics("workflow.png")

## ----eval = FALSE-------------------------------------------------------------
# library(zellkonverter)
# library(SingleCellExperiment)
# library(cellGeometry)
# 
# typist_h5 <- readH5AD("2ac906a5-9725-4258-8e36-21a9f6c0302a.h5ad",
#                       use_hdf5 = TRUE, reader = "R")

## ----eval = FALSE-------------------------------------------------------------
# mat <- typist_h5@assays@data$X
# rownames(mat) <- rownames(typist_h5)
# meta <- typist_h5@colData@listData

## ----eval=FALSE---------------------------------------------------------------
# library(Seurat)
# typist <- readRDS("2ac906a5-9725-4258-8e36-21a9f6c0302a.rds")  # 15.5 GB in memory
# 
# mat <- typist@assays$RNA$counts
# meta <- typist@meta.data

## ----eval = FALSE-------------------------------------------------------------
# table(meta$Majority_voting_CellTypist)
# 
# subcl <- meta$Majority_voting_CellTypist
# cellgrp <- meta$Majority_voting_CellTypist_high
# 
# # reduce dataset to only blood (optional)
# subcl[meta$tissue != "blood"] <- NA
# cellgrp[meta$tissue != "blood"] <- NA

## ----eval = FALSE-------------------------------------------------------------
# mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
#                   dual_mean = TRUE, cores = 2)

## ----eval = FALSE-------------------------------------------------------------
# # example code using future for parallelisation on windows
# library(future)
# plan(multisession, workers = 4)
# 
# mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp,
#                   use_future = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# library(AnnotationHub)
# ah <- AnnotationHub()
# ensDb_v110 <- ah[["AH113665"]]
# mk <- gene2symbol(mk, ensDb_v110)

## ----eval = FALSE-------------------------------------------------------------
# signature_heatmap(mk)  # visualise whole signature
# signature_heatmap(mk, top = 5)  # show top 5 genes for each subclass

## ----out.width='90%', echo=FALSE----------------------------------------------
knitr::include_graphics("typist_sig.png")

## ----eval = FALSE-------------------------------------------------------------
# spillover_heatmap(mk)

## ----out.width='80%', echo=FALSE----------------------------------------------
knitr::include_graphics("typist_spillover.png")

## ----eval = FALSE-------------------------------------------------------------
# cs <- cos_similarity(mk)
# ComplexHeatmap::Heatmap(cs)

## ----out.width='80%', echo=FALSE----------------------------------------------
knitr::include_graphics("typist_cos_sim.png")

## ----eval = FALSE-------------------------------------------------------------
# rank_angle(cs)

## ----eval = FALSE-------------------------------------------------------------
# diagnose(mk)

## ----eval = FALSE-------------------------------------------------------------
# mk <- updateMarkers(mk,
#                     remove_subclass = c("Helper T cells", "Cytotoxic T cells"))

## ----eval = FALSE-------------------------------------------------------------
# # simulated bulk
# set.seed(3)
# sim_counts <- generate_samples(mk, 25)
# sim_percent <- sim_counts / rowSums(sim_counts) * 100
# sim_pseudo <- simulate_bulk(mk, sim_counts)

## ----eval = FALSE-------------------------------------------------------------
# # mode 1: (perfect deconvolution)
# fit <- deconvolute(mk, sim_pseudo,
#                    use_filter = FALSE)
# plot_set(sim_counts, fit$subclass$output)
# plot_set(sim_percent, fit$subclass$percent)
# 
# metric_set(sim_percent, fit$subclass$percent)  # table of results

## ----eval = FALSE-------------------------------------------------------------
# # mode 2: sample from original sc count matrix
# # 1.43 mins (Intel); 45 secs (ARM)
# set.seed(99)
# times <- 1  # can be increased
# sim_sampled <- simulate_bulk(mat, sim_counts, subcl, times = times)
# 
# # fix rownames
# rownames(sim_sampled) <- gene2symbol(rownames(sim_sampled), ensDb_v110)
# 
# # near optimal deconvolution of counts sampled from the original scRNA-Seq
# fit2 <- deconvolute(mk, sim_sampled,
#                     use_filter = FALSE, arith_mean = TRUE)
# 
# # plot results
# plot_set(sim_counts, fit2$subclass$output / times)  # adjust for oversampling using `times`
# plot_set(sim_percent, fit2$subclass$percent)
# 
# metric_set(sim_percent, fit2$subclass$percent)

## ----out.width='90%', echo=FALSE----------------------------------------------
knitr::include_graphics("typist_sim.png")

## ----eval = FALSE-------------------------------------------------------------
# mk <- updateMarkers(mk, bulkdata = my_bulk_matrix)

## ----eval = FALSE-------------------------------------------------------------
# res <- tune_deconv(mk, sim_sampled, sim_counts * times,
#                    grid = list(nsubclass = c(5, 10, 15, 25, 50, 100, 200, 500, 1000),
#                                expfilter = c(0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5),
#                                weight_method = c("none", "equal")),
#                    arith_mean = TRUE,
#                    use_filter = FALSE,
#                    cores = 8)
# 
# ## Tuning parameters: nsubclass, expfilter, weight_method
# ##   |==========================================================| 100%, Elapsed 00:07
# ## Best tune:
# ##   nsubclass  expfilter  weight_method  mean.RMSE
# ##        1000       0.25          equal      35.21

## ----eval = FALSE-------------------------------------------------------------
# plot_tune(res, xvar = "nsubclass", group = "weight_method")
# plot_tune(res, xvar = "nsubclass", group = "expfilter")
# plot_tune(res, xvar = "expfilter", group = "weight_method")
# plot_tune(res, xvar = "expfilter", group = "nsubclass")

## ----out.width='90%', echo=FALSE----------------------------------------------
knitr::include_graphics("typist_tune.png")

## ----eval = FALSE-------------------------------------------------------------
# # simple gaussian noise applied to counts
# sim_noise <- add_noise(sim_sampled)
# 
# # noise applied to log2 counts
# sim_noise <- log_noise(sim_sampled)
# 
# # noise applied to sqrt transformed counts
# sim_noise <- sqrt_noise(sim_sampled)
# 
# # whole genes are scaled up/down by a random amount
# # this simulates differences in chemistry
# sim_noise <- shift_noise(sim_sampled)

## ----eval = FALSE-------------------------------------------------------------
# res2 <- tune_deconv(mk, sim_noise, sim_counts,
#                    grid = list(nsubclass = c(5, 10, 15, 25, 50, 100, 200, 500, 1000),
#                                expfilter = c(0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5),
#                                weight_method = c("none", "equal")),
#                    arith_mean = TRUE,
#                    use_filter = FALSE,
#                    cores = 8)

