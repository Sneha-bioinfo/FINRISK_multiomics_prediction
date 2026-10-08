# ==============================================================================
# LOAD PACKAGES
# ==============================================================================
library("SummarizedExperiment")
library("MultiAssayExperiment")
library("TreeSummarizedExperiment")
library("S4Vectors")
library("data.table")
library("survival")
library("mia")
library('ggplot2')
library('cowplot')
library("splitTools")
library("IntegratedLearner")

# ==============================================================================
# TAXA (FINRISK taxonomic assay)
# ==============================================================================
tse_taxa <- readRDS("data/tse_metaphlan4_rel_ab_SGB.rds")

tse_species <- altExp(tse_taxa, "Species")
assay_taxa <- assay(tse_species, "rel_ab")
taxa_rownames <- sub("^s__", "", rownames(assay_taxa))
assay_taxa <- apply(assay_taxa, 2, as.numeric)
rownames(assay_taxa) <- taxa_rownames

taxa_coldata <- as.data.frame(colData(tse_taxa))
stopifnot(identical(rownames(taxa_coldata), as.character(taxa_coldata$Barcode)))
stopifnot(identical(rownames(taxa_coldata), colnames(tse_taxa)))
rownames(taxa_coldata) <- colnames(tse_taxa)
colnames(assay_taxa) <- as.character(taxa_coldata[colnames(assay_taxa), "FID"]) #setting colnames of assay to FID

# ==============================================================================
# METABOLOMICS (FINRISK metabolomics assay)
# ==============================================================================
metab <- fread("data/FR02_Microbiome_NMR_2022-06-14.txt")
metab$FID <- as.character(metab$FID)

meta_cols <- c("FID", "SAMPLE_COLLECTION", "SPECTROMETER")
metab_mat <- metab[, setdiff(colnames(metab), meta_cols), with = FALSE]
metab_mat <- as.matrix(metab_mat)
metab_mat <- t(metab_mat)
colnames(metab_mat) <- metab$FID
stopifnot(all(colnames(metab_mat) == metab$FID))
rn <- rownames(metab_mat)
remove_idx <- grepl("^TOTAL", rn) | grepl("TAG$", rn) | grepl("UNSATURATION$", rn) |
  grepl("SIZE$", rn) | grepl("PCT$", rn) | grepl("BY", rn)
metab_mat <- metab_mat[!remove_idx, ]
metab_mat <- metab_mat[rowSums(is.na(metab_mat)) < ncol(metab_mat), ]

metab_feature_names <- rownames(metab_mat)
fid_names <- colnames(metab_mat)
metab_mat <- matrix(as.numeric(metab_mat), nrow = nrow(metab_mat), ncol = ncol(metab_mat),
                    dimnames = list(metab_feature_names, fid_names))


# ==============================================================================
# REMOVING SAMPLES THAT HAVE NA FOR ENTIRE FEATURE SET
# ==============================================================================
cat("Taxa NA samples for entire feature set:", sum(colSums(!is.na(assay_taxa)) == 0), "\n") 
cat("Metab NA samples for entire feature set:", sum(colSums(!is.na(metab_mat)) == 0), "\n") 

# Removing samples that have NA for entire feature set
assay_taxa <- assay_taxa[, colSums(!is.na(assay_taxa))>0]
metab_mat <- metab_mat[, colSums(!is.na(metab_mat))>0]
cat("Metab samples after removing samples with all NA", dim(metab_mat), "\n") 

#Sanity checks
# check to see if any sample with NA for entire feature set is still present
cat("Metab NA samples for entire feature set:", sum(colSums(!is.na(metab_mat)) == 0), "\n")

# Checking if any sample in taxa has NA for some of the feature but not all
cat("Taxa - samples with SOME (not all) NA features remaining:", sum(colSums(is.na(assay_taxa)) > 0 & colSums(is.na(assay_taxa)) < nrow(assay_taxa)), "\n") 

# Checking if any sample in Metab has NA for some of the feature but not all
cat("Metab - samples with SOME (not all) NA features remaining:", sum(colSums(is.na(metab_mat)) > 0 & colSums(is.na(metab_mat)) < nrow(metab_mat)), "\n") 

# Converting NA to 0 for those samples which may have NA for few features but not ALL in metab
metab_mat[is.na(metab_mat)] <- 0


# ==============================================================================
# FINRISK AMR (phenotype/label source only); Cause-of-death labels (K_TPKS) are extracted from this assay. Samples are then matched to the taxonomy and metabolomics assays by FID, restricting the analysis to samples common to all layers.
# ==============================================================================

tse_amr       <- readRDS("data/TSE_FINRISK_AMR.rds")
amr_coldata   <- as.data.frame(colData(tse_amr))
rownames(amr_coldata) <- as.character(amr_coldata$Row.names)
stopifnot(!any(duplicated(rownames(amr_coldata))))

amr_coldata$FID <- taxa_coldata[rownames(amr_coldata), "FID"] # map AMR Barcode -> FID via taxa_coldata
cat("AMR - unmatched to taxa:", sum(is.na(amr_coldata$FID)), "\n")

AMR_valid_samples <- !is.na(amr_coldata$FID)
amr_coldata <- amr_coldata[AMR_valid_samples, ]
rownames(amr_coldata) <- as.character(amr_coldata$FID) # re-key by FID

# ==============================================================================
# TAXA-METAB COMMON SAMPLES ( FID key)
# ==============================================================================
final_fids <- intersect(colnames(assay_taxa), colnames(metab_mat))
cat("Common samples across taxa + metabolomics:", length(final_fids), "\n")

# ==============================================================================
# ENDPOINT - Multiclass cause of death 
# Logic: K_TPKS checked; no match -> NA -> dropped 
# ==============================================================================
meta <- amr_coldata
meta <- meta[intersect(final_fids, rownames(meta)), ] # restrict to taxa/metab-covered samples that also have AMR-derived labels
cat("meta (multiclass, from AMR) samples before cause filtering:", nrow(meta), "\n")

stopifnot("K_TPKS" %in% colnames(meta))

cat("Total entries in K_TPKS:", length(meta$K_TPKS), "\n")
cat("Total NA in K_TPKS:", sum(is.na(meta$K_TPKS)), "\n")
cat("Non-NA entries:", sum(!is.na(meta$K_TPKS)), "\n")

meta$K_TPKS_letter <- gsub("[0-9]+", "", meta$K_TPKS)

meta$cause_multiclass <- dplyr::case_when(
  meta$K_TPKS_letter == "C" ~ "Cancer",
  meta$K_TPKS_letter == "I" ~ "Cardiovascular",
  !is.na(meta$K_TPKS_letter) ~ "Other",
  TRUE ~ NA_character_
)

cat("Distribution before dropping NA (unmatched/alive/censored):\n")
print(table(meta$cause_multiclass, useNA = "ifany"))

meta <- meta[!is.na(meta$cause_multiclass), ]
meta$cause_multiclass <- factor(meta$cause_multiclass)

cat("Cause-of-death class distribution:\n")
print(table(meta$cause_multiclass))

final_fids <- rownames(meta)
event_vec  <- meta$cause_multiclass
cat("Samples retained after cause filtering:", length(final_fids), "\n")

# ==============================================================================
# SUBSET TAXA + METAB TO FINAL CAUSE-MATCHED SAMPLES (from AMR)
# ==============================================================================
assay_taxa_s <- assay_taxa[, final_fids]
metab_mat_s  <- metab_mat[, final_fids]

stopifnot(identical(colnames(assay_taxa_s), final_fids))
stopifnot(identical(colnames(metab_mat_s), final_fids))

# ==============================================================================
# WRAP IN TSES
# ==============================================================================
taxa_tse_il  <- TreeSummarizedExperiment(assays = list(taxa_relative_abundance = assay_taxa_s))
colData(taxa_tse_il)$subjectID <- final_fids; colData(taxa_tse_il)$Y <- event_vec

metab_tse_il <- TreeSummarizedExperiment(assays = list(metabolite_abundance = metab_mat_s))
colData(metab_tse_il)$subjectID <- final_fids; colData(metab_tse_il)$Y <- event_vec

# ==============================================================================
# PREVALENCE FILTERING
# ==============================================================================
taxa_tse_il  <- subsetByPrevalent(taxa_tse_il, prevalence = 0.10, detection = 0.001, assay.type = "taxa_relative_abundance")
metab_tse_il <- subsetByPrevalent(metab_tse_il, prevalence = 0.10, detection = 0.001, assay.type = "metabolite_abundance")

cat("After prevalence filtering - Taxa:", nrow(taxa_tse_il), "Metab:", nrow(metab_tse_il), "\n")

# ==============================================================================
# BUILD MAE
# ==============================================================================
cd <- S4Vectors::DataFrame(Y = event_vec, subjectID = final_fids, row.names = final_fids)

experiment_names <- c("taxonomy", "metabolomics")
smap <- S4Vectors::DataFrame(
  assay   = rep(experiment_names, each = length(final_fids)),
  primary = rep(final_fids, times = length(experiment_names)),
  colname = rep(final_fids, times = length(experiment_names))
)

mae <- MultiAssayExperiment(
  experiments = ExperimentList(taxonomy = taxa_tse_il, metabolomics = metab_tse_il),
  colData = cd, sampleMap = smap
)
all_ids <- rownames(colData(mae))
stopifnot(identical(all_ids, final_fids))

# ==============================================================================
# OUTPUT DIRECTORY
# ==============================================================================
set.seed(1) # change to shuffle train-test split
output_dir <- "outputs/multiclass/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# STRATIFICATION - sex + cause_multiclass ; MEN
# ==============================================================================
combined <- multi_strata(
  data.frame(sex   = meta[all_ids, "MEN"],
             event = meta[all_ids, "cause_multiclass"]),
  k = 4
)
inds      <- partition(combined,
                       p = c(train = 0.7, test = 0.3),
                       split_into_list = FALSE,
                       type = "stratified",
                       seed = 1)
train_ids <- all_ids[inds == "train"]
valid_ids <- all_ids[inds == "test"]

# Sanity Checks to check gender balance in train-test data and within classes

cat("=== Stratified by sex + cause_multiclass ===\n")
cat("Train class distribution:\n")
print(table(meta[train_ids, "cause_multiclass"]))
cat("Test class distribution:\n")
print(table(meta[valid_ids, "cause_multiclass"]))

cat("Train set Gender ratio: \n")
print(prop.table(table(meta[train_ids, "MEN"])))
cat("Test set Gender ratio: \n")
print(prop.table(table(meta[valid_ids, "MEN"])))
cat("Full dataset gender ratio: \n")
print(prop.table(table(meta[all_ids,"MEN"])))

# --- Gender balance inside the classes-
cat("Gender composition within each cause - Full:\n")
print(round(prop.table(table(meta[all_ids, "MEN"], meta[all_ids, "cause_multiclass"]), margin = 2), 3))

cat("Gender composition within each cause - Train:\n")
print(round(prop.table(table(meta[train_ids, "MEN"], meta[train_ids, "cause_multiclass"]), margin = 2), 3))

cat("Gender composition within each cause - Test:\n")
print(round(prop.table(table(meta[valid_ids, "MEN"], meta[valid_ids, "cause_multiclass"]), margin = 2), 3))


# ==============================================================================
# RUN MODEL - Multiclass classification
# ==============================================================================
mae_train <- mae[, train_ids]
mae_valid <- mae[, valid_ids]

fit_all <- IntegratedLearner(
  MAE_train = mae_train, MAE_valid = mae_valid,
  experiment = c("taxonomy","metabolomics"),
  assay.type = c("taxa_relative_abundance","metabolite_abundance"),
  folds = 5, base_learner = "xgboost", meta_learner = "xgboost",
  filter_method = NULL, filter_pct = NULL, run_screening = FALSE,
  drop_poor_performing_layers = FALSE, family = stats::binomial(),
  verbose = TRUE
)

save(fit_all, file = file.path(output_dir, "taxa_metab_xgboost.RData"))
capture.output(fit_all, file = file.path(output_dir, "taxa_metab_xgboost_output.txt"))

#-- plot---#
plot_obj <- IntegratedLearner:::plot.learner(fit_all)

ggsave(
  filename = file.path(output_dir, "model_performance_plot.png"),
  plot = plot_obj$plot,
  width = 8,
  height = 6,
  dpi = 300
)

sessionInfo()
