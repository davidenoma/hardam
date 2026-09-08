### ──────────────────────────────────────────────────────────────────────────────
### 📌 Main Pipeline#!/usr/bin/env Rscript
### ──────────────────────────────────────────────────────────────────────────────
### 📌 Load Required Libraries
### ──────────────────────────────────────────────────────────────────────────────
library(bigsnpr)    # Only used for recoding PLINK if needed
library(glmnet)
library(grpreg)
library(dplyr)
library(data.table)
library(stringr)

"%&%" <- function(a, b) paste(a, b, sep = "")

FIXED_SEED <- 42L
## Make RNG reproducible globally
set.seed(FIXED_SEED)
### ──────────────────────────────────────────────────────────────────────────────
### 📌 Define Command-Line Arguments
### ──────────────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript pipeline.R <plink_prefix> <har_bed_file> <gene_annot_file> <gene_expr_file> <chrom> [mode]")
}

plink_prefix       <- args[1]      # GTEx genotype prefix (e.g., "gtex_geno")
har_bed_file       <- args[2]      # HAR locations in BED format
gene_annot_file    <- args[3]      # Gene annotation file (tab-delimited)
gene_expr_file     <- args[4]      # GTEx gene expression file (with metadata in first 4 columns)
chrom              <- as.character(args[5])  # Chromosome number (e.g., "1")
mode               <- ifelse(length(args) >= 6, args[6], "har_only")  # Analysis mode

chrom <- gsub("chr", "", chrom)

# Validate mode argument
valid_modes <- c("har_only", "har_cis", "cis_only")
if (!mode %in% valid_modes) {
  stop("Invalid mode. Use one of: ", paste(valid_modes, collapse = ", "))
}

# Set flags based on mode
include_cis <- mode == "har_cis"
cis_only <- mode == "cis_only"

message("🔧 Configuration:")
message("  - Analysis mode: ", toupper(mode))
message("  - Include cis-SNPs: ", include_cis)
message("  - Cis-only mode: ", cis_only)

### ──────────────────────────────────────────────────────────────────────────────
### 📌 Define Output Filenames (Chromosome-specific)
### ──────────────────────────────────────────────────────────────────────────────
tissue_name  <- tools::file_path_sans_ext(basename(gene_expr_file))
tissue_name <- gsub("\\.v8\\.normalized_expression", "", tissue_name)  # Clean GTEx naming
# Unique prefix per chromosome
prefix_name  <- paste0("chr", chrom, "_", tissue_name)
results_file <- paste0("model_results_", prefix_name, ".txt")
weights_file <- paste0("snp_weights_", prefix_name, ".txt")
pred_file    <- paste0("predicted_expression_", prefix_name, ".txt")


### 📌 Functions to Load PLINK Data

load_bim_file <- function(bim_file) {
  bim_dt <- fread(bim_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bim_dt) <- c("chrom", "snp_id", "cm_pos", "physical_pos", "allele1", "allele2")
  return(bim_dt)
}

# Recode PLINK genotype files to raw numeric format.
plink_recode_raw <- function(geno_prefix,chrom) {
  message("🔄 Converting PLINK files to numeric (.raw) format...")

  new_prefix <- paste0(geno_prefix, "_", tissue_name,"_",chrom)
  cmd <- paste("plink --bfile", geno_prefix, "--chr",chrom,"--recode A --out", new_prefix)

  system(cmd)
  rawfile <- paste0(geno_prefix, "_raw.raw")
  message("✅ .raw file created: ", rawfile)
  return(rawfile)
}

load_plink_genotype <- function(geno_prefix,chrom) {
  message(paste("🟢 Loading PLINK genotype:", geno_prefix))

  bim <- load_bim_file(paste0(geno_prefix, ".bim"))
  fam <- fread(paste0(geno_prefix, ".fam"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
  colnames(fam) <- c("FID", "IID", "father", "mother", "sex", "phenotype")

  raw_file <- plink_recode_raw(geno_prefix,chrom)
  geno_mat_dt <- fread(raw_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
  snp_cols <- colnames(geno_mat_dt)[7:ncol(geno_mat_dt)]
  geno_mat <- as.matrix(geno_mat_dt[, ..snp_cols])
  clean_snp_names <- sub("_[ATGCatgc]+$", "", snp_cols)
  colnames(geno_mat) <- clean_snp_names
  rownames(geno_mat) <- as.character(geno_mat_dt$IID)
  message(paste("✅ Loaded genotype matrix:", nrow(geno_mat), "individuals x", ncol(geno_mat), "SNPs"))

  return(list(geno = geno_mat, bim = bim, fam = fam))
}

### ──────────────────────────────────────────────────────────────────────────────
### 📌 ilter GTEx BIM Based on Target SNPs
### ──────────────────────────────────────────────────────────────────────────────
filter_gtex_by_target_snps <- function(gtex_prefix, target_bim) {
  message("🔍 Filtering GTEx genotype using target SNPs (PLINK)...")
  tmp_snp_list <- paste0(gtex_prefix, "_target_snps.txt")
  fwrite(target_bim[, .(snp_id)], file = tmp_snp_list, col.names = FALSE, row.names = FALSE, quote = FALSE)
  out_prefix <- paste0(gtex_prefix, "_filtered")
  cmd <- paste("plink --bfile", gtex_prefix,
               "--extract", tmp_snp_list,
               "--make-bed --out", out_prefix)
  system(cmd)
  message("✅ PLINK filtering complete: ", out_prefix)
  return(out_prefix)
}

### ──────────────────────────────────────────────────────────────────────────────
### 📌 Gene Annotation & Expression Functions
### ──────────────────────────────────────────────────────────────────────────────
get_gene_annotation <- function(gene_annot_file, chrom) {
  message("📄 Loading Gene Annotations...")

  gene_df <- read.table(gene_annot_file, header = FALSE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE)
  gene_df1 <- dplyr::filter(gene_df, V3 %in% "gene")
  geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
  genename <- gsub("gene_name (\\S+);", "\\1", str_extract(gene_df1[,9], "gene_name (\\S+);"), perl = TRUE)
  gene_used <- as.data.frame(cbind(geneid, genename, gene_df1[, c(1,4,5,3)]))
  colnames(gene_used) <- c("geneid", "genename", "chr", "start", "end", "anno")
  gtf_used <- dplyr::filter(gene_used, gene_used[,3] %in% ('chr' %&% chrom))
  return(gtf_used)
}

get_gene_expression <- function(gene_expr_file, gene_annot) {
  message("🔍 Loading GTEx gene expression file...")
  expr_df <- fread(gene_expr_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(expr_df)[1:4] <- c("chr", "start", "end", "gene_id")
  expr_df <- dplyr::filter(expr_df, gene_id %in% gene_annot$geneid)

  sample_ids <- colnames(expr_df)[5:ncol(expr_df)]
  mat_expr <- t(as.matrix(expr_df[, 5:ncol(expr_df), with = FALSE]))
  colnames(mat_expr) <- expr_df$gene_id
  rownames(mat_expr) <- sample_ids
  message("✅ Loaded expression matrix: ", nrow(mat_expr), " samples x ", ncol(mat_expr), " genes.")
  return(mat_expr)
}

match_sample_ids <- function(geno_prefix, expr_mat) {
  message("🔍 Matching sample IDs between genotype and expression data...")

  # Create unique suffix for this job using SLURM_JOB_ID (if available) and process ID
  unique_suffix <- paste0(
    Sys.getenv("SLURM_JOB_ID", ""),
    "_",
    Sys.getpid()
  )

  # Setup cleanup function for intermediate files
  cleanup_files <- function(prefix) {
    files_to_remove <- paste0(prefix, c(".bed", ".bim", ".fam", ".log", ".nosex"))
    unlink(files_to_remove)
    message("✨ Cleaned up intermediate files")
  }

  tryCatch({
    fam <- fread(paste0(geno_prefix, ".fam"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
    colnames(fam) <- c("FID", "IID", "father", "mother", "sex", "phenotype")
    geno_samples <- as.character(fam$IID)
    expr_samples <- rownames(expr_mat)
    common_samples <- intersect(geno_samples, expr_samples)
    if (length(common_samples) == 0) {
      stop("🚨 No matching sample IDs found!")
    }
    message("✅ Matched samples: ", length(common_samples))

    # Use unique suffix for intermediate files
    matched_prefix <- paste0(geno_prefix, "_", tissue_name, "_", unique_suffix)
    matched_fam_file <- paste0(matched_prefix, ".fam")

    # Write temporary .fam file with matched samples
    matched_fam <- fam[fam$IID %in% common_samples, ]
    fwrite(matched_fam, file = matched_fam_file, sep = " ", col.names = FALSE)

    # Create PLINK files with matched samples
    cmd <- paste("plink --bfile", geno_prefix, "--keep", matched_fam_file,
                 "--make-bed --out", matched_prefix)
    system_status <- system(cmd)
    if (system_status != 0) {
      stop("PLINK command failed")
    }

    message("✅ New PLINK genotype file with matched samples: ", matched_prefix)

    # Register cleanup on exit (will run even if error occurs later)
    reg.finalizer(environment(), function(e) {
      cleanup_files(matched_prefix)
    }, onexit = TRUE)

    return(matched_prefix)
  }, error = function(e) {
    # Clean up any files that might have been created before the error
    cleanup_files(paste0(geno_prefix, "_", tissue_name, "_", unique_suffix))
    stop(e)  # Re-throw the error after cleanup
  })
}

extend_har_regions <- function(bim_dt, har_bed_file, chr) {
  message("🔄 Extending HAR regions to ensure at least 10 SNPs...")
  har_regions <- fread(har_bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(har_regions) <- c("chrom", "start", "end", "har_id")
  har_regions <- har_regions %>%
    mutate(chrom = gsub("^chr", "", chrom)) %>%
    filter(chrom == chr)
  total_hars <- nrow(har_regions)
  for (i in seq_len(total_hars)) {
    har_id <- har_regions$har_id[i]
    start <- har_regions$start[i]
    end <- har_regions$end[i]
    region_snps <- bim_dt %>% filter(chrom == chr & physical_pos >= start & physical_pos <= end)
    while (nrow(region_snps) < 10) {
      start <- max(0, start - 500)
      end <- end + 500
      if ((end - start) > 5e4) {
        # message("⚠️ Skipping HAR: ", har_id, " due to excessive expansion (>50kb).")
        break
      }
      region_snps <- bim_dt %>% filter(chrom == chr & physical_pos >= start & physical_pos <= end)
    }
    har_regions$start[i] <- start
    har_regions$end[i] <- end
  }
  message("✅ HAR regions extended successfully! Total: ", nrow(har_regions))
  return(har_regions)
}

filter_har_snps <- function(bim_dt, har_regions) {
  message("🔍 Filtering SNPs within the extended HAR regions...")
  if (is.null(har_regions) || nrow(har_regions) == 0) {
    message("⚠️ No HAR regions available for filtering.")
    return(NULL)
  }
  bim_dt <- bim_dt %>% mutate(chrom = as.character(chrom))
  har_regions <- har_regions %>% mutate(chrom = as.character(chrom))

  har_snp_list <- list()
  for (i in seq_len(nrow(har_regions))) {
    region <- har_regions[i, ]
    snps <- bim_dt %>% filter(chrom == region$chrom &
                                physical_pos >= region$start &
                                physical_pos <= region$end)
    if (nrow(snps) == 0) next
    snps$har_id <- region$har_id
    har_snp_list[[i]] <- snps
  }
  if (length(har_snp_list) == 0) {
    message("⚠️ No HAR SNPs found.")
    return(NULL)
  }
  har_snp_dt <- rbindlist(har_snp_list, use.names = TRUE, fill = TRUE)
  har_snp_dt <- har_snp_dt %>%
    group_by(har_id) %>%
    mutate(har_group_id = paste0(har_id, "_", row_number())) %>%
    ungroup()
  message("✅ Total HAR SNPs Found: ", nrow(har_snp_dt))
  return(as.data.frame(har_snp_dt))
}

apply_group_lasso <- function(geno_har, expr_vec, har_snp_info) {
  message("🔍 Running Group Lasso for HAR grouping...")
  matching_idx <- match(colnames(geno_har), har_snp_info$snp_id)
  valid_idx <- which(!is.na(matching_idx))
  if (length(valid_idx) == 0) {
    stop("🚨 No genotype SNP IDs could be matched to HAR SNP info!")
  }
  geno_har <- geno_har[, valid_idx, drop = FALSE]
  matching_idx <- matching_idx[valid_idx]

  group_factor <- har_snp_info$har_id[matching_idx]
  group_factor <- as.factor(group_factor)
  message("Group factor levels: ", paste(levels(group_factor), collapse = ", "))

  if (length(unique(group_factor)) < 2) {
    message("⚠️ Not enough HAR groups (<2) for Group Lasso. Skipping...")
    return(character(0))
  }

  cvfit <- tryCatch({
    set.seed(FIXED_SEED)
    cv.grpreg(X = as.matrix(geno_har), y = expr_vec, group = group_factor, nfolds = 10)
  }, error = function(e) {
    message("Error in Group Lasso: ", e$message)
    return(NULL)
  })

  if (is.null(cvfit)) return(character(0))

  message("Group Lasso Coefficients at lambda.min:")
  print(coef(cvfit))

  selected_groups <- predict(cvfit, type = "groups", lambda = cvfit$lambda.min)
  selected_idx <- which(group_factor %in% selected_groups)
  selected_snps <- colnames(geno_har)[selected_idx]
  message("Selected SNPs from nonzero groups: ", paste(selected_snps, collapse = ", "))
  return(selected_snps)
}


### ──────────────────────────────────────────────────────────────────────────────
### 📌 Run Elastic Net on Selected HAR SNPs
do_elastic_net <- function(geno, expr, alpha = 0.5) {
  message("🚀 Running Elastic Net...")
  set.seed(FIXED_SEED)
  fit <- cv.glmnet(as.matrix(geno), expr, alpha = alpha, nfolds = 10, type.measure = "mse")
  best_lambda <- fit$lambda.min
  return(list(cv_fit = fit, best_lambda = best_lambda))
}
### 📌 Function to Extract Cis-SNPs for Genes
### ──────────────────────────────────────────────────────────────────────────────
extract_cis_snps <- function(bim_dt, gene_info, cis_window = 1e6) {
  message("🔍 Extracting cis-SNPs for gene: ", gene_info$genename)

  gene_chrom <- gsub("^chr", "", gene_info$chr)
  gene_start <- gene_info$start
  gene_end <- gene_info$end

  # Define cis region (gene boundaries +/- window)
  cis_start <- max(0, gene_start - cis_window)
  cis_end <- gene_end + cis_window

  # Filter SNPs within cis region
  cis_snps <- bim_dt %>%
    filter(chrom == gene_chrom &
           physical_pos >= cis_start &
           physical_pos <= cis_end)

  if (nrow(cis_snps) > 0) {
    cis_snps$cis_type <- "cis"
    message("✅ Found ", nrow(cis_snps), " cis-SNPs for gene: ", gene_info$genename)
  } else {
    message("⚠️ No cis-SNPs found for gene: ", gene_info$genename)
  }

  return(cis_snps)
}
### ──────────────────────────────────────────────────────────────────────────────
### 📌 Function to Combine HAR and Cis-SNPs
### ──────────────────────────────────────────────────────────────────────────────

combine_har_cis_snps <- function(har_selected_snps, cis_snps, har_snp_info) {
  message("🔄 Combining selected HAR SNPs and cis-SNPs...")

  # Subset HAR info for selected SNPs
  har_selected_info <- har_snp_info[har_snp_info$snp_id %in% har_selected_snps, , drop = FALSE]
  if (nrow(har_selected_info) > 0) {
    har_selected_info$snp_source <- "har"
  }

  # Annotate cis SNPs
  if (!is.null(cis_snps) && nrow(cis_snps) > 0) {
    cis_snps$snp_source <- "cis"
    if (!"har_group_id" %in% names(cis_snps)) cis_snps$har_group_id <- paste0("cis_", seq_len(nrow(cis_snps)))
  }

  # Collect non-empty parts
  parts <- list()
  if (nrow(har_selected_info) > 0) parts[[length(parts) + 1]] <- data.table::as.data.table(har_selected_info)
  if (!is.null(cis_snps) && nrow(cis_snps) > 0) parts[[length(parts) + 1]] <- data.table::as.data.table(cis_snps)

  # If both empty, return empty HAR structure
  if (length(parts) == 0) {
    combined_snp_info <- har_snp_info[0, , drop = FALSE]
    message("✅ Combined SNPs: 0 total (0 HAR + 0 cis)")
    return(combined_snp_info)
  }

  # Robust row bind with fill for missing columns
  combined <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)

  # Prefer HAR over cis on duplicates
  if ("snp_id" %in% names(combined)) {
    o <- order(factor(combined$snp_source, levels = c("har", "cis")))
    combined <- combined[o, , drop = FALSE]
    combined <- combined[!duplicated(combined$snp_id), , drop = FALSE]
  }

  combined_snp_info <- as.data.frame(combined)

  message("✅ Combined SNPs: ", nrow(combined_snp_info), " total (",
          sum(combined_snp_info$snp_source == "har", na.rm = TRUE), " HAR + ",
          sum(combined_snp_info$snp_source == "cis", na.rm = TRUE), " cis)")
  return(combined_snp_info)
}

### ──────────────────────────────────────────────────────────────────────────────
### 📌 Evaluate Model Performance
evaluate_performance <- function(geno, expr, fit, best_lambda) {
  pred_expr <- as.vector(predict(fit$cv_fit$glmnet.fit, as.matrix(geno), s = best_lambda))

  # Preserve sample names
  if (is.null(names(pred_expr)) || length(names(pred_expr)) == 0) {
    names(pred_expr) <- rownames(geno)
  }

  # Align lengths
  pred_expr <- pred_expr[rownames(geno)]
  expr <- expr[rownames(geno)]

  R2 <- 1 - sum((expr - pred_expr)^2) / sum((expr - mean(expr))^2)

  nonzero_coef <- coef(fit$cv_fit$glmnet.fit, s = best_lambda)
  selected_snps <- rownames(nonzero_coef)[nonzero_coef[, 1] != 0]

  return(list(R2 = R2, best_lambda = best_lambda, selected_snps = selected_snps, pred_expr = pred_expr))
}


### ──────────────────────────────────────────────────────────────────────────────
### 📌 Save Model Results, SNP Weights, & Predicted Expression
### ──────────────────────────────────────────────────────────────────────────────

### ──────────────────────────────────────────────────────────────────────────────
### 📌 Save Model Results, SNP Weights, & Predicted Expression
### ──────────────────────────────────────────────────────────────────────────────

save_results <- function(gene, gene_name, R2, best_lambda, num_snps) {
  # Create output directory if it doesn't exist
  output_dir <- "output_results_seed"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add suffix based on analysis mode
  suffix <- switch(mode,
                  "har_only" = "",
                  "har_cis" = "_with_cis",
                  "cis_only" = "_cis_only")

  shared_file <- file.path(output_dir, paste0("model_results_", tissue_name, suffix, ".txt"))

  header <- "Gene\tGeneName\tR2\tBestLambda\tNumSNPs"
  if (!file.exists(shared_file)) {
    write(header, file = shared_file, append = FALSE)
  }
  out_line <- paste(gene, gene_name, R2, best_lambda, num_snps, sep = "\t")
  write(out_line, file = shared_file, append = TRUE)
}

save_weights <- function(gene, fit, best_lambda, snp_info) {
  output_dir <- "output_results_seed"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  suffix <- switch(mode,
                  "har_only" = "",
                  "har_cis" = "_with_cis",
                  "cis_only" = "_cis_only")
  shared_file <- file.path(output_dir, paste0("snp_weights_", tissue_name, suffix, ".txt"))

  co <- coef(fit$cv_fit$glmnet.fit, s = best_lambda)
  w_df <- data.frame(SNP = rownames(co), Weight = as.numeric(co), stringsAsFactors = FALSE)
  w_df <- w_df[w_df$Weight != 0, ]

  if (nrow(w_df) > 0) {
    merged <- merge(w_df, snp_info, by.x = "SNP", by.y = "snp_id", all.x = TRUE)
    merged$Gene <- gene

    if (mode == "cis_only") {
      header <- "SNP\tWeight\tchrom\tcm_pos\tphysical_pos\tallele1\tallele2\tsnp_source\tGene"
      cols <- c("SNP", "Weight", "chrom", "cm_pos", "physical_pos", "allele1", "allele2", "snp_source", "Gene")
    } else if (mode == "har_cis") {
      header <- "SNP\tWeight\tchrom\tcm_pos\tphysical_pos\tallele1\tallele2\thar_id\thar_group_id\tsnp_source\tGene"
      cols <- c("SNP", "Weight", "chrom", "cm_pos", "physical_pos", "allele1", "allele2", "har_id", "har_group_id", "snp_source", "Gene")
    } else {
      # fallback for har_only or other
      header <- "SNP\tWeight\tchrom\tcm_pos\tphysical_pos\tallele1\tallele2\thar_id\thar_group_id\tGene"
      cols <- c("SNP", "Weight", "chrom", "cm_pos", "physical_pos", "allele1", "allele2", "har_id", "har_group_id", "Gene")
    }

    if (!file.exists(shared_file)) {
      write(header, file = shared_file, append = FALSE)
    }
    write.table(merged[, cols, drop = FALSE],
                file      = shared_file,
                append    = TRUE,
                sep       = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote     = FALSE)
  }
}

save_predictions <- function(gene, gene_name, chrom, pred_expr, actual_expr) {
  # Create output directory if it doesn't exist
    output_dir <- "output_results_seed"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add suffix based on analysis mode
  suffix <- switch(mode,
                  "har_only" = "",
                  "har_cis" = "_with_cis",
                  "cis_only" = "_cis_only")

  shared_file <- file.path(output_dir, paste0("predicted_expression_", tissue_name, suffix, ".txt"))

  df <- data.frame(
    Gene               = gene,
    GeneName           = gene_name,
    Chrom              = chrom,
    Sample             = names(pred_expr),
    PredictedExpression = as.numeric(pred_expr),
    ActualExpression   = as.numeric(actual_expr[names(pred_expr)]),
    stringsAsFactors   = FALSE
  )
  header <- "Gene\tGeneName\tChrom\tSample\tPredictedExpression\tActualExpression"
  if (!file.exists(shared_file)) {
    write(header, file = shared_file, append = FALSE)
  }
  write.table(df,
              file      = shared_file,
              append    = TRUE,
              sep       = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}


### ──────────────────────────────────────────────────────────────────────────────
### 📌 Main Pipeline
run_pipeline <- function() {
  message("🔄 Starting the Pipeline...")

  # Load GTEx BIM file
  gtex_bim <- load_bim_file(paste0(plink_prefix, ".bim"))

  # Extend HAR regions using the GTEx BIM
  har_regions <- extend_har_regions(gtex_bim, har_bed_file, chrom)

  # Filter HAR SNPs based on extended HAR regions & label with HAR group IDs
  har_snps <- filter_har_snps(gtex_bim, har_regions)
  if (is.null(har_snps) || nrow(har_snps) == 0) {
    stop("🚨 No HAR SNPs found. Exiting.")
  }


  # Load gene annotations and gene expression data
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df <- get_gene_expression(gene_expr_file, gene_annot)

  # Construct the expected matched file prefix from plink_prefix and tissue_name
  expected_matched_prefix <- paste0(plink_prefix, "_", tissue_name, "_matched")
  matched_bed <- paste0(expected_matched_prefix, ".bed")

  if (!file.exists(matched_bed)) {
    message("Matched PLINK file not found. Running match_sample_ids...")
    matched_prefix <- match_sample_ids(plink_prefix, expr_df)
  } else {
    message("Matched PLINK file found. Using: ", expected_matched_prefix)
    matched_prefix <- expected_matched_prefix
  }

  # Construct the raw file path produced by PLINK recoding
  raw_file <- paste0(matched_prefix, "_raw.raw")
  if (!file.exists(raw_file)) {
    message("Raw file not found. Recoding PLINK file...")
    system(paste("plink --bfile", matched_prefix, "--recode A --out", paste0(matched_prefix, "_raw")))
  } else {
    message("Raw file already exists: ", raw_file)
  }
  geno_mat_dt <- fread(raw_file, header = TRUE, sep = " ", stringsAsFactors = FALSE)
  snp_cols <- colnames(geno_mat_dt)[7:ncol(geno_mat_dt)]
  geno_mat <- as.matrix(geno_mat_dt[, ..snp_cols])
  clean_snp_names <- sub("_[ATGCatgc]+$", "", snp_cols)
  colnames(geno_mat) <- clean_snp_names
  rownames(geno_mat) <- as.character(geno_mat_dt$IID)
  message("Genotype matrix dimensions: ", paste(dim(geno_mat), collapse = " x "))

    # ─────────────── Clean up intermediate genotype files ───────────────
  message("🧹 Cleaning up intermediate PLINK and .raw files...")
  cleanup_patterns <- c(".bed", ".bim", ".fam", ".log", ".nosex", ".raw")
  for (ext in cleanup_patterns) {
    files_to_remove <- paste0(matched_prefix, ext)
    if (file.exists(files_to_remove)) unlink(files_to_remove)
  }
  message("✅ Temporary genotype files removed.")

  # Filter genotype matrix for HAR SNPs using SNP IDs from BIM
  har_snp_ids <- har_snps$snp_id
  message("HAR SNP IDs to filter:")

  common_snp_ids <- intersect(colnames(geno_mat), har_snp_ids)
  if (length(common_snp_ids) == 0) {
    stop("🚨 No HAR SNPs found in genotype matrix!")
  }
  geno_har <- geno_mat[, common_snp_ids, drop = FALSE]
  message("Filtered HAR genotype matrix dimensions: ", paste(dim(geno_har), collapse = " x "))

  # Do not rename the genotype matrix columns; use har_snp_info to create a grouping factor.
  har_snp_info <- har_snps[har_snps$snp_id %in% common_snp_ids, ]
  message("HAR SNP info (first 5 rows):")


  for (gene in colnames(expr_df)) {
    gene_info <- gene_annot %>% filter(geneid == gene)
    gene_name <- ifelse(nrow(gene_info) > 0, gene_info$genename[1], gene)
    expr_vec <- as.numeric(expr_df[, gene])
    names(expr_vec) <- rownames(expr_df)

    message("🧬 Processing gene: ", gene_name, " (", toupper(mode), ")")

    # Skip if gene annotation is missing
    if (cis_only && nrow(gene_info) == 0) {
      message("⚠️ No annotation for ", gene_name, "; skipping.")
      next
    }

    # ─────────────── Mode 1: CIS-ONLY ───────────────
    if (cis_only) {
      cis_snps <- extract_cis_snps(gtex_bim, gene_info[1, ], cis_window = 5e5)
      if (nrow(cis_snps) == 0) next

      common_snps <- intersect(colnames(geno_mat), cis_snps$snp_id)
      if (length(common_snps) == 0) next

      samples <- intersect(rownames(geno_mat), names(expr_vec))
      final_geno <- geno_mat[samples, common_snps, drop = FALSE]
      expr_vec <- expr_vec[samples]

      final_snp_info <- cis_snps[cis_snps$snp_id %in% common_snps, ]
      final_snp_info$snp_source <- "cis"
      final_snp_info$har_id <- "cis_region"
      final_snp_info$har_group_id <- paste0("cis_", seq_len(nrow(final_snp_info)))

# ─────────────── Mode 2–3: HAR-ONLY / HAR+CIS ───────────────
    } else {
      #Performing GROUP LASSO on HAR SNPs
      samples <- intersect(rownames(geno_har), names(expr_vec))
      expr_vec <- expr_vec[samples]
      geno_subset <- geno_har[samples, , drop = FALSE]

      selected_snps <- apply_group_lasso(geno_subset, expr_vec, har_snp_info)
      if (length(selected_snps) == 0) {
        message("⚠️ No HAR groups selected for ", gene_name)
        next
      }

      # ─────────────── HAR+CIS ───────────────
      if (include_cis) {
        cis_snps <- extract_cis_snps(gtex_bim, gene_info[1, ], cis_window = 5e5)
        combined_info <- combine_har_cis_snps(selected_snps, cis_snps, har_snp_info)

        # Ensure ≥1 HAR SNP remains
        if (!any(combined_info$snp_source == "har", na.rm = TRUE)) {
          message("⚠️ No HAR SNPs left after combining; skipping ", gene_name)
          next
        }

        final_snps <- intersect(colnames(geno_mat), combined_info$snp_id)
        final_geno <- geno_mat[samples, final_snps, drop = FALSE]
        final_snp_info <- combined_info[combined_info$snp_id %in% final_snps, ]

      # ─────────────── HAR-ONLY ───────────────
      } else {
        final_geno <- geno_subset[, selected_snps, drop = FALSE]
        final_snp_info <- har_snp_info[har_snp_info$snp_id %in% selected_snps, ]
      }
    }

    # ─────────────── Run Elastic Net & Save ───────────────
    if (ncol(final_geno) < 2) {
      message("⚠️ Fewer than 2 SNPs selected for ", gene_name, "; skipping elastic net.")
      next
    }
    message("🚀 Elastic Net: ", nrow(final_geno), " samples × ", ncol(final_geno), " SNPs")
    model <- do_elastic_net(final_geno, expr_vec)
    perf <- evaluate_performance(final_geno, expr_vec, model, model$best_lambda)

    save_results(gene, gene_name, perf$R2, model$best_lambda, ncol(final_geno))
    save_weights(gene, model, model$best_lambda, final_snp_info)
    save_predictions(gene, gene_name, chrom, perf$pred_expr, expr_vec)
  }
}

run_pipeline()
message("🎯 Pipeline complete for chromosome ", chrom, " (", toupper(mode), ")")
