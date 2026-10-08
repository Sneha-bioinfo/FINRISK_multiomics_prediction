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
# COMMON SAMPLES ACROSS BOTH LAYERS (FID is the shared key)
# ==============================================================================
cat("Taxa original no of samples", dim(assay_taxa), "\n") 
cat("Metab original no of sample", dim(metab_mat), "\n") 
final_fids <- intersect(colnames(assay_taxa), colnames(metab_mat))
cat("Common samples across taxa + metabolomics:", length(final_fids), "\n") 
assay_taxa_s <- assay_taxa[, final_fids] 
metab_mat_s  <- metab_mat[, final_fids] 

stopifnot(identical(colnames(assay_taxa_s), final_fids))
stopifnot(identical(colnames(metab_mat_s), final_fids))

# ==============================================================================
# ENDPOINT (the event/outcome of interest)
# ==============================================================================
meta <- taxa_coldata
rownames(meta) <- meta$FID
meta <- meta[final_fids, ]
cat("meta after subset to common_ids:", dim(meta), "\n")

time_vec  <- meta$DEATH_AGEDIFF #FINRISK survival time
event_vec <- meta$DEATH #FINRISK all cause mortality- coded as DEATH- 1; CENSOR- 0

cat("Samples:", length(final_fids), "\n")
cat("Events: ", sum(event_vec == 1, na.rm = TRUE), "\n")
cat("Censored:", sum(event_vec == 0, na.rm = TRUE), "\n")

# ==============================================================================
# WRAP IN TSEs (building two separate tse objects for taxonomy and metabolomics using only the required columns from the original assays)
# ==============================================================================
taxa_tse_il <- TreeSummarizedExperiment(assays = list(taxa_relative_abundance = assay_taxa_s))
colData(taxa_tse_il)$subjectID <- final_fids
colData(taxa_tse_il)$Y         <- Surv(time_vec, event_vec)
colData(taxa_tse_il)$time      <- time_vec
colData(taxa_tse_il)$event     <- event_vec

metab_tse_il <- TreeSummarizedExperiment(assays = list(metabolite_abundance = metab_mat_s))
colData(metab_tse_il)$subjectID <- final_fids
colData(metab_tse_il)$Y         <- Surv(time_vec, event_vec)
colData(metab_tse_il)$time      <- time_vec
colData(metab_tse_il)$event     <- event_vec

# ==============================================================================
# PREVALENCE FILTERING- 
# ==============================================================================
taxa_tse_il  <- subsetByPrevalent(taxa_tse_il,  prevalence = 0.10, detection = 0.001, assay.type = "taxa_relative_abundance")
metab_tse_il <- subsetByPrevalent(metab_tse_il, prevalence = 0.10, detection = 0.001, assay.type = "metabolite_abundance")

cat("After prevalence filtering:\n")
cat("Taxa:", nrow(taxa_tse_il), "Metab:", nrow(metab_tse_il), "\n")

# ==============================================================================
# TRANSFORMATION
# ==============================================================================
taxa_tse_il  <- transformAssay(taxa_tse_il,  method = "pa",       assay.type = "taxa_relative_abundance", pseudocount = FALSE,  name = "taxa_pa") # To use taxa_pa assay if needed to evaluate tool on taxonomy presence absence asssay later

# ==============================================================================
# BUILD MAE
# ==============================================================================
cd <- S4Vectors::DataFrame(time      = time_vec,
                           event     = event_vec,
                           subjectID = final_fids,
                           row.names = final_fids)
cd$Y <- Surv(time_vec, event_vec)

experiment_names <- c("taxonomy", "metabolomics")
smap <- S4Vectors::DataFrame(
  assay   = as.character(rep(experiment_names, each  = length(final_fids))),
  primary = as.character(rep(final_fids,       times = length(experiment_names))),
  colname = as.character(rep(final_fids,       times = length(experiment_names)))
)

mae <- MultiAssayExperiment(
  experiments = ExperimentList(taxonomy = taxa_tse_il, metabolomics = metab_tse_il),
  colData   = cd,
  sampleMap = smap
)
all_ids <- rownames(colData(mae))
stopifnot(identical(all_ids, final_fids))

# ==============================================================================
# OUTPUT DIRECTORY - change as needed
# ==============================================================================
set.seed(1) # change to shuffle train-test samples per run
output_dir <- "outputs/survival/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# STRATIFICATION - sex + event (to maintain FINRISK cohort's 45-55% men to women ratio, as well as event ratio in the train and test cohorts too); cohort split into 30-70% train and test data respectively; MEN- 1- indicates MEN, 0 indicates WOMEN, DEATH- 1 event (death); 0- censored
# ==============================================================================
combined <- multi_strata(
  data.frame(sex   = meta[all_ids, "MEN"],
             event = meta[all_ids, "DEATH"]),
  k = 4
)
inds      <- partition(combined,
                       p    = c(train = 0.7, test = 0.3),
                       split_into_list = FALSE,
                       type = "stratified",
                       seed = 1)
train_ids <- all_ids[inds == "train"]
valid_ids <- all_ids[inds == "test"]

cat("=== Stratified by sex + event ===\n")
cat("Train - Women:", sum(meta[train_ids, "MEN"] == 0),
    "Men:",           sum(meta[train_ids, "MEN"] == 1),
    "Events:",        sum(meta[train_ids, "DEATH"] == 1),
    "Censored:",      sum(meta[train_ids, "DEATH"] ==0), "\n")
cat("Test  - Women:", sum(meta[valid_ids, "MEN"] == 0),
    "Men:",           sum(meta[valid_ids, "MEN"] == 1),
    "Events:",        sum(meta[valid_ids, "DEATH"] == 1),
    "Censored:",      sum(meta[valid_ids, "DEATH"] ==0), "\n")

# ==============================================================================
# RUN MODEL (IntegratedLearner framework)
# ==============================================================================
mae_train <- mae[, train_ids]
mae_valid <- mae[, valid_ids]

fit_all <- IntegratedLearner(
  MAE_train=mae_train, MAE_valid=mae_valid,
  experiment=c("taxonomy","metabolomics"),
  assay.type=c("taxa_relative_abundance","metabolite_abundance"),
  folds=5, base_learner = "surv.xgboost.aft", filter_method=NULL, 
  filter_pct=NULL, run_screening=FALSE, drop_poor_performing_layers = FALSE, verbose=TRUE)
save(fit_all, file=file.path(output_dir, "taxa_metab_xgboost.aft.RData"))



# -- plot---#
plot_obj <- IntegratedLearner:::plot.learner(fit_all)

ggsave(
  filename = file.path(output_dir, "model_performance_plot.png"),
  plot = plot_obj$plot,
  width = 8,
  height = 6,
  dpi = 300
)

sessionInfo()
