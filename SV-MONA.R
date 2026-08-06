# ============================================================================
# SV-MONA: Single-Variable Multi-Omics Association and Network Analysis

## Integrative multi-omics association and network analysis
# using average CPITN scores. The workflow can also be applied to other continuous or
# ordinal single variables by replacing the CPITN input.
# Author: Xiaojia He, Ph.D.
# Institution: UL Research Institutes' Chemical Insights
#              
# Version: 1.0.0
# Date: 2026-08-06
#
# License: MIT
#
# GitHub:
# https://github.com/OrangeBluee/SV-MONA

# Description:
#   Performs Spearman association analysis between average CPITN scores and
#   omics_1, omics_2, and omics_3 features; selects associated features;
#   constructs a CPITN-centered multi-omics network; calculates network
#   topology; and exports files for Cytoscape and downstream visualization.
#
# Input structure:
#   - Rows: features
#   - Columns: samples
#   - First column: feature identifier
#   - CPITN input: exactly one feature row (for example, single_variable)
#
# Required packages:
#   - dplyr
#   - igraph
#
# Repository layout:
#   project_root/
#   ├── input/
#   │   ├── single_variable.csv
#   │   ├── omics_1.csv
#   │   ├── omics_2.csv
#   │   └── omics_3.csv
#   └── output/
#
# Notes:
#   - Omics values are analyzed without transformation by this script.
#   - Spearman correlation is used for primary associations.
#   - No maximum feature limit is imposed.
# ============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "igraph"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
})

# -----------------------------------------------------------------------------
# 1. Project paths and input files
# -----------------------------------------------------------------------------

# Run the script from the repository root, or set project_dir explicitly.
project_dir <- normalizePath(".", mustWork = TRUE)
input_dir <- file.path(project_dir, "input")
network_dir <- file.path(project_dir, "output")

file_cpitn <- file.path(input_dir, "single_variable.csv")
file_omics_1 <- file.path(input_dir, "omics_1.csv")
file_omics_2 <- file.path(input_dir, "omics_2.csv")
file_omics_3 <- file.path(input_dir, "omics_3.csv")

if (!dir.exists(network_dir)) {
  dir.create(
    network_dir,
    recursive = TRUE
  )
}

#----------------------------------------------------------------------
# 2. Confirm that all files exist
#----------------------------------------------------------------------

input_files <- c(
  CPITN = file_cpitn,
  omics_1 = file_omics_1,
  omics_2 = file_omics_2,
  omics_3 = file_omics_3
)

file_check <- data.frame(
  dataset = names(input_files),
  file = unname(input_files),
  exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

print(file_check)

if (any(!file_check$exists)) {
  
  missing <- file_check$file[
    !file_check$exists
  ]
  
  stop(
    paste0(
      "The following files were not found:\n",
      paste(missing, collapse = "\n")
    )
  )
}

#----------------------------------------------------------------------
# 3. Function to read a feature-by-sample matrix
#----------------------------------------------------------------------

read_feature_matrix <- function(file) {
  
  x <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(x) < 2) {
    stop(
      "The file must contain a feature-name column ",
      "and at least one sample column: ",
      file
    )
  }
  
  feature_names <- trimws(
    as.character(x[[1]])
  )
  
  x <- x[, -1, drop = FALSE]
  
  # Remove leading and trailing spaces from sample names
  colnames(x) <- trimws(
    colnames(x)
  )
  
  # Convert sample values to numeric
  x[] <- lapply(
    x,
    function(z) {
      suppressWarnings(
        as.numeric(z)
      )
    }
  )
  
  # Ensure unique feature identifiers
  rownames(x) <- make.unique(
    feature_names
  )
  
  output <- as.matrix(x)
  
  storage.mode(output) <- "numeric"
  
  output
}

#----------------------------------------------------------------------
# 4. Read all datasets
#----------------------------------------------------------------------

cpitn_mat <- read_feature_matrix(
  file_cpitn
)

omics_1_mat <- read_feature_matrix(
  file_omics_1
)

omics_2_mat <- read_feature_matrix(
  file_omics_2
)

omics_3_mat <- read_feature_matrix(
  file_omics_3
)

cat(
  "\nDimensions before sample alignment:\n"
)

cat(
  "CPITN:",
  nrow(cpitn_mat),
  "features x",
  ncol(cpitn_mat),
  "samples\n"
)

cat(
  "omics_1:",
  nrow(omics_1_mat),
  "features x",
  ncol(omics_1_mat),
  "samples\n"
)

cat(
  "omics_2:",
  nrow(omics_2_mat),
  "features x",
  ncol(omics_2_mat),
  "samples\n"
)

cat(
  "omics_3:",
  nrow(omics_3_mat),
  "features x",
  ncol(omics_3_mat),
  "samples\n"
)

#----------------------------------------------------------------------
# 5. Extract CPITN average
#----------------------------------------------------------------------

if (nrow(cpitn_mat) != 1) {
  
  stop(
    "The CPITN file must contain exactly one feature row. ",
    "It currently contains ",
    nrow(cpitn_mat),
    " rows."
  )
}

cpitn <- as.numeric(
  cpitn_mat[1, ]
)

names(cpitn) <- colnames(
  cpitn_mat
)

cpitn_name <- rownames(
  cpitn_mat
)[1]

if (
  is.na(cpitn_name) ||
  cpitn_name == ""
) {
  cpitn_name <- "single_variable"
}

cat(
  "\nCPITN feature:",
  cpitn_name,
  "\n"
)

#----------------------------------------------------------------------
# 6. Check shared samples
#----------------------------------------------------------------------

shared_samples_all <- Reduce(
  intersect,
  list(
    names(cpitn),
    colnames(omics_1_mat),
    colnames(omics_2_mat),
    colnames(omics_3_mat)
  )
)

cat(
  "Samples shared across all four datasets:",
  length(shared_samples_all),
  "\n"
)

print(shared_samples_all)

if (length(shared_samples_all) < 5) {
  stop(
    "Fewer than five samples are shared across all datasets. ",
    "Check the sample column names."
  )
}

#----------------------------------------------------------------------
# 7. Analysis parameters
#----------------------------------------------------------------------

# Minimum complete samples required for an association
min_complete <- 10

# CPITN-feature association threshold
cpitn_rho_threshold <- 0.4
cpitn_p_threshold <- 0.05

# Omics-feature pairwise association threshold
omics_rho_threshold <- 0.4
omics_p_threshold <- 0.05

# Pairwise block size
#
# Larger values may run faster but require more memory.
# Values between 250 and 1000 are usually reasonable.
pairwise_block_size <- 500

# No maximum number of features is imposed.

#----------------------------------------------------------------------
# 8. Function: associate one omics block with CPITN
#----------------------------------------------------------------------

associate_block_with_cpitn <- function(
    feature_matrix,
    cpitn,
    omics_name,
    min_complete = 10
) {
  
  feature_matrix <- as.matrix(
    feature_matrix
  )
  
  common_samples <- intersect(
    names(cpitn),
    colnames(feature_matrix)
  )
  
  if (length(common_samples) < min_complete) {
    
    stop(
      omics_name,
      " has fewer than ",
      min_complete,
      " shared samples with CPITN."
    )
  }
  
  y <- cpitn[
    common_samples
  ]
  
  xmat <- feature_matrix[
    ,
    common_samples,
    drop = FALSE
  ]
  
  results <- vector(
    "list",
    nrow(xmat)
  )
  
  for (i in seq_len(nrow(xmat))) {
    
    x <- as.numeric(
      xmat[i, ]
    )
    
    keep <- complete.cases(
      x,
      y
    )
    
    n_complete <- sum(
      keep
    )
    
    feature_name <- rownames(
      xmat
    )[i]
    
    current_result <- data.frame(
      omics = omics_name,
      feature = feature_name,
      n = n_complete,
      rho = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
    
    if (n_complete < min_complete) {
      results[[i]] <- current_result
      next
    }
    
    x_complete <- x[keep]
    y_complete <- y[keep]
    
    if (
      length(unique(x_complete)) < 2 ||
      length(unique(y_complete)) < 2
    ) {
      results[[i]] <- current_result
      next
    }
    
    test <- suppressWarnings(
      cor.test(
        x_complete,
        y_complete,
        method = "spearman",
        exact = FALSE
      )
    )
    
    current_result$rho <- unname(
      test$estimate
    )
    
    current_result$p_value <- test$p.value
    
    results[[i]] <- current_result
  }
  
  results <- bind_rows(
    results
  )
  
  results$fdr_block <- p.adjust(
    results$p_value,
    method = "BH"
  )
  
  results$direction <- case_when(
    results$rho > 0 ~ "Positive",
    results$rho < 0 ~ "Negative",
    TRUE ~ NA_character_
  )
  
  results %>%
    arrange(
      p_value,
      desc(abs(rho))
    )
}

#----------------------------------------------------------------------
# 9. Run CPITN association analysis
#----------------------------------------------------------------------

cat(
  "\nRunning CPITN-omics_1 associations...\n"
)

omics_1_results <- associate_block_with_cpitn(
  feature_matrix = omics_1_mat,
  cpitn = cpitn,
  omics_name = "omics_1",
  min_complete = min_complete
)

cat(
  "Running CPITN-omics_2 associations...\n"
)

omics_2_results <- associate_block_with_cpitn(
  feature_matrix = omics_2_mat,
  cpitn = cpitn,
  omics_name = "omics_2",
  min_complete = min_complete
)

cat(
  "Running CPITN-omics_3 associations...\n"
)

omics_3_results <- associate_block_with_cpitn(
  feature_matrix = omics_3_mat,
  cpitn = cpitn,
  omics_name = "omics_3",
  min_complete = min_complete
)

#----------------------------------------------------------------------
# 10. Combine all CPITN association results
#----------------------------------------------------------------------

all_results <- bind_rows(
  omics_1_results,
  omics_2_results,
  omics_3_results
)

all_results$fdr_global <- p.adjust(
  all_results$p_value,
  method = "BH"
)

all_results <- all_results %>%
  mutate(
    abs_rho = abs(rho),
    
    passes_cpitn_threshold =
      !is.na(rho) &
      !is.na(p_value) &
      abs(rho) >= cpitn_rho_threshold &
      p_value < cpitn_p_threshold,
    
    passes_block_fdr =
      !is.na(rho) &
      !is.na(fdr_block) &
      abs(rho) >= cpitn_rho_threshold &
      fdr_block < 0.05,
    
    passes_global_fdr =
      !is.na(rho) &
      !is.na(fdr_global) &
      abs(rho) >= cpitn_rho_threshold &
      fdr_global < 0.05
  ) %>%
  arrange(
    p_value,
    desc(abs_rho)
  )

#----------------------------------------------------------------------
# 11. Export all CPITN association results
#----------------------------------------------------------------------

write.csv(
  all_results,
  file.path(
    network_dir,
    "CPITN_all_omics_associations.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  omics_1_results,
  file.path(
    network_dir,
    "CPITN_omics_1_associations.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  omics_2_results,
  file.path(
    network_dir,
    "CPITN_omics_2_associations.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  omics_3_results,
  file.path(
    network_dir,
    "CPITN_omics_3_associations.csv"
  ),
  row.names = FALSE,
  na = ""
)

#----------------------------------------------------------------------
# 12. Select every feature passing CPITN thresholds
#
# No maximum feature limit is used.
#----------------------------------------------------------------------

candidate_features <- all_results %>%
  filter(
    passes_cpitn_threshold
  ) %>%
  arrange(
    omics,
    p_value,
    desc(abs_rho)
  )

candidate_summary <- candidate_features %>%
  count(
    omics,
    name = "selected_features"
  )

print(candidate_summary)

cat(
  "\nTotal selected features:",
  nrow(candidate_features),
  "\n"
)

write.csv(
  candidate_features,
  file.path(
    network_dir,
    "CPITN_selected_features.csv"
  ),
  row.names = FALSE,
  na = ""
)

if (nrow(candidate_features) == 0) {
  
  stop(
    "No features passed the CPITN thresholds. ",
    "Consider reducing cpitn_rho_threshold or ",
    "increasing cpitn_p_threshold."
  )
}

#----------------------------------------------------------------------
# 13. Extract selected features from each omics matrix
#----------------------------------------------------------------------

extract_selected_features <- function(
    feature_matrix,
    selected_features,
    omics_name
) {
  
  selected_features <- intersect(
    selected_features,
    rownames(feature_matrix)
  )
  
  if (length(selected_features) == 0) {
    return(NULL)
  }
  
  result <- feature_matrix[
    selected_features,
    ,
    drop = FALSE
  ]
  
  rownames(result) <- paste(
    omics_name,
    selected_features,
    sep = "::"
  )
  
  result
}

selected_omics_1 <- extract_selected_features(
  omics_1_mat,
  candidate_features$feature[
    candidate_features$omics == "omics_1"
  ],
  "omics_1"
)

selected_omics_2 <- extract_selected_features(
  omics_2_mat,
  candidate_features$feature[
    candidate_features$omics == "omics_2"
  ],
  "omics_2"
)

selected_omics_3 <- extract_selected_features(
  omics_3_mat,
  candidate_features$feature[
    candidate_features$omics == "omics_3"
  ],
  "omics_3"
)

selected_matrix_list <- list(
  omics_1 = selected_omics_1,
  omics_2 = selected_omics_2,
  omics_3 = selected_omics_3
)

selected_matrix_list <- selected_matrix_list[
  !vapply(
    selected_matrix_list,
    is.null,
    logical(1)
  )
]

#----------------------------------------------------------------------
# 14. Align selected feature matrices
#----------------------------------------------------------------------

common_selected_samples <- Reduce(
  intersect,
  c(
    list(names(cpitn)),
    lapply(
      selected_matrix_list,
      colnames
    )
  )
)

if (
  length(common_selected_samples) <
  min_complete
) {
  
  stop(
    "Fewer than ",
    min_complete,
    " common samples remain after alignment."
  )
}

selected_matrix_list <- lapply(
  selected_matrix_list,
  function(x) {
    
    x[
      ,
      common_selected_samples,
      drop = FALSE
    ]
  }
)

selected_feature_matrix <- do.call(
  rbind,
  selected_matrix_list
)

storage.mode(
  selected_feature_matrix
) <- "numeric"

cpitn_network <- cpitn[
  common_selected_samples
]

cat(
  "\nSelected network matrix:\n"
)

cat(
  nrow(selected_feature_matrix),
  "features x",
  ncol(selected_feature_matrix),
  "samples\n"
)

#----------------------------------------------------------------------
# 15. Create CPITN-feature edge table
#----------------------------------------------------------------------

clinical_node_id <- paste0(
  "Clinical::",
  cpitn_name
)

cpitn_edges <- candidate_features %>%
  transmute(
    source = clinical_node_id,
    target = paste(
      omics,
      feature,
      sep = "::"
    ),
    source_type = "Clinical",
    target_type = omics,
    edge_type = "CPITN-omics",
    association = rho,
    abs_association = abs(rho),
    sign = ifelse(
      rho > 0,
      "Positive",
      "Negative"
    ),
    weight = abs(rho),
    distance = 1 / pmax(
      abs(rho),
      0.0001
    ),
    n = n,
    p_value = p_value,
    fdr_block = fdr_block,
    fdr_global = fdr_global,
    pairwise_bonferroni = NA_real_,
    pairwise_fdr_retained = NA_real_
  )

#----------------------------------------------------------------------
# 16. Blockwise pairwise Spearman correlation
#
# This evaluates every selected feature pair.
# It does not impose a maximum number of features.
#
# Only edges satisfying the rho and nominal-p thresholds are retained
# in memory. A global Bonferroni value is calculated using the total
# number of possible feature-pair tests.
#----------------------------------------------------------------------

blockwise_spearman_edges <- function(
    feature_matrix,
    rho_threshold = 0.70,
    p_threshold = 0.05,
    min_complete = 10,
    block_size = 500
) {
  
  feature_matrix <- as.matrix(
    feature_matrix
  )
  
  storage.mode(
    feature_matrix
  ) <- "numeric"
  
  feature_names <- rownames(
    feature_matrix
  )
  
  number_features <- nrow(
    feature_matrix
  )
  
  if (number_features < 2) {
    return(data.frame())
  }
  
  total_pairwise_tests <- choose(
    number_features,
    2
  )
  
  block_starts <- seq(
    1,
    number_features,
    by = block_size
  )
  
  retained_blocks <- list()
  
  output_counter <- 1L
  
  for (
    block_i in seq_along(block_starts)
  ) {
    
    i_start <- block_starts[
      block_i
    ]
    
    i_end <- min(
      i_start + block_size - 1,
      number_features
    )
    
    index_i <- i_start:i_end
    
    matrix_i <- feature_matrix[
      index_i,
      ,
      drop = FALSE
    ]
    
    available_i <- !is.na(
      matrix_i
    )
    
    for (
      block_j in block_i:length(block_starts)
    ) {
      
      j_start <- block_starts[
        block_j
      ]
      
      j_end <- min(
        j_start + block_size - 1,
        number_features
      )
      
      index_j <- j_start:j_end
      
      matrix_j <- feature_matrix[
        index_j,
        ,
        drop = FALSE
      ]
      
      available_j <- !is.na(
        matrix_j
      )
      
      # Spearman correlation between feature blocks
      correlation_block <- suppressWarnings(
        cor(
          t(matrix_i),
          t(matrix_j),
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
      
      # Number of complete paired observations
      n_block <- available_i %*%
        t(available_j)
      
      keep <- !is.na(
        correlation_block
      ) &
        abs(correlation_block) >=
        rho_threshold &
        n_block >= min_complete
      
      # For within-block comparisons, keep only
      # the upper triangle
      if (block_i == block_j) {
        
        keep <- keep &
          upper.tri(
            correlation_block,
            diag = FALSE
          )
      }
      
      locations <- which(
        keep,
        arr.ind = TRUE
      )
      
      if (nrow(locations) == 0) {
        next
      }
      
      rho_values <- correlation_block[
        locations
      ]
      
      n_values <- n_block[
        locations
      ]
      
      # Approximate two-sided p-values
      t_statistics <- rho_values *
        sqrt(
          (n_values - 2) /
            pmax(
              1 - rho_values^2,
              .Machine$double.eps
            )
        )
      
      p_values <- 2 * pt(
        -abs(t_statistics),
        df = n_values - 2
      )
      
      pass_p <- !is.na(p_values) &
        p_values < p_threshold
      
      if (!any(pass_p)) {
        next
      }
      
      locations <- locations[
        pass_p,
        ,
        drop = FALSE
      ]
      
      rho_values <- rho_values[
        pass_p
      ]
      
      n_values <- n_values[
        pass_p
      ]
      
      p_values <- p_values[
        pass_p
      ]
      
      current_edges <- data.frame(
        source = feature_names[
          index_i[
            locations[, 1]
          ]
        ],
        target = feature_names[
          index_j[
            locations[, 2]
          ]
        ],
        n = as.integer(
          n_values
        ),
        association = as.numeric(
          rho_values
        ),
        p_value = as.numeric(
          p_values
        ),
        stringsAsFactors = FALSE
      )
      
      retained_blocks[[
        output_counter
      ]] <- current_edges
      
      output_counter <-
        output_counter + 1L
      
      cat(
        "Processed blocks",
        block_i,
        "and",
        block_j,
        "- retained",
        nrow(current_edges),
        "edges\n"
      )
    }
  }
  
  if (length(retained_blocks) == 0) {
    return(data.frame())
  }
  
  result <- bind_rows(
    retained_blocks
  )
  
  # Conservative adjustment across all possible pairs
  result$pairwise_bonferroni <- pmin(
    result$p_value *
      total_pairwise_tests,
    1
  )
  
  # BH adjustment among retained strong candidate edges
  #
  # This is not equivalent to BH correction across every
  # possible pair because only strong candidate edges are retained.
  result$pairwise_fdr_retained <- p.adjust(
    result$p_value,
    method = "BH"
  )
  
  result$total_pairwise_tests <-
    total_pairwise_tests
  
  result
}

#----------------------------------------------------------------------
# 17. Run all selected feature-pair correlations
#----------------------------------------------------------------------

cat(
  "\nRunning blockwise omics-omics correlations...\n"
)

cat(
  "Total possible feature pairs:",
  choose(
    nrow(selected_feature_matrix),
    2
  ),
  "\n"
)

omics_pairwise_results <- blockwise_spearman_edges(
  feature_matrix =
    selected_feature_matrix,
  rho_threshold =
    omics_rho_threshold,
  p_threshold =
    omics_p_threshold,
  min_complete =
    min_complete,
  block_size =
    pairwise_block_size
)

#----------------------------------------------------------------------
# 18. Format omics-omics edges
#----------------------------------------------------------------------

if (nrow(omics_pairwise_results) > 0) {
  
  omics_edges <- omics_pairwise_results %>%
    mutate(
      source_type =
        sub(
          "::.*$",
          "",
          source
        ),
      
      target_type =
        sub(
          "::.*$",
          "",
          target
        ),
      
      edge_type = case_when(
        source_type ==
          target_type ~
          paste0(
            source_type,
            "-within"
          ),
        
        TRUE ~
          "Cross-omics"
      ),
      
      abs_association =
        abs(association),
      
      sign = ifelse(
        association > 0,
        "Positive",
        "Negative"
      ),
      
      weight =
        abs_association,
      
      distance =
        1 / pmax(
          abs_association,
          0.0001
        ),
      
      fdr_block =
        NA_real_,
      
      fdr_global =
        NA_real_
    ) %>%
    select(
      source,
      target,
      source_type,
      target_type,
      edge_type,
      association,
      abs_association,
      sign,
      weight,
      distance,
      n,
      p_value,
      fdr_block,
      fdr_global,
      pairwise_bonferroni,
      pairwise_fdr_retained,
      total_pairwise_tests
    )
  
} else {
  
  omics_edges <- data.frame()
}

#----------------------------------------------------------------------
# 19. Combine CPITN and omics-omics edges
#----------------------------------------------------------------------

cytoscape_edges <- bind_rows(
  cpitn_edges,
  omics_edges
) %>%
  distinct(
    source,
    target,
    edge_type,
    .keep_all = TRUE
  )

cat(
  "\nNetwork edge counts:\n"
)

cat(
  "CPITN-feature edges:",
  nrow(cpitn_edges),
  "\n"
)

cat(
  "Omics-omics edges:",
  nrow(omics_edges),
  "\n"
)

cat(
  "Total edges:",
  nrow(cytoscape_edges),
  "\n"
)

#----------------------------------------------------------------------
# 20. Create node table
#----------------------------------------------------------------------

omics_nodes <- data.frame(
  node_id = rownames(
    selected_feature_matrix
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    node_type =
      sub(
        "::.*$",
        "",
        node_id
      ),
    
    feature =
      sub(
        "^[^:]+::",
        "",
        node_id
      ),
    
    label = feature
  )

clinical_node <- data.frame(
  node_id = clinical_node_id,
  node_type = "Clinical",
  feature = cpitn_name,
  label = cpitn_name,
  stringsAsFactors = FALSE
)

cytoscape_nodes <- bind_rows(
  clinical_node,
  omics_nodes
)

#----------------------------------------------------------------------
# 21. Add CPITN association statistics to nodes
#----------------------------------------------------------------------

node_cpitn_statistics <- candidate_features %>%
  transmute(
    node_id = paste(
      omics,
      feature,
      sep = "::"
    ),
    cpitn_rho = rho,
    cpitn_abs_rho = abs(rho),
    cpitn_p_value = p_value,
    cpitn_fdr_block = fdr_block,
    cpitn_fdr_global = fdr_global,
    cpitn_direction = direction
  )

cytoscape_nodes <- cytoscape_nodes %>%
  left_join(
    node_cpitn_statistics,
    by = "node_id"
  )

cytoscape_nodes$cpitn_rho[
  cytoscape_nodes$node_type ==
    "Clinical"
] <- 1

cytoscape_nodes$cpitn_abs_rho[
  cytoscape_nodes$node_type ==
    "Clinical"
] <- 1

cytoscape_nodes$cpitn_direction[
  cytoscape_nodes$node_type ==
    "Clinical"
] <- "Clinical"

#----------------------------------------------------------------------
# 22. Build igraph object
#----------------------------------------------------------------------

network_graph <- graph_from_data_frame(
  d = cytoscape_edges,
  directed = FALSE,
  vertices = cytoscape_nodes
)

E(network_graph)$weight <-
  as.numeric(
    E(network_graph)$weight
  )

E(network_graph)$distance <-
  as.numeric(
    E(network_graph)$distance
  )

#----------------------------------------------------------------------
# 23. Calculate node-level topology
#----------------------------------------------------------------------

node_degree <- degree(
  network_graph,
  mode = "all"
)

node_strength <- strength(
  network_graph,
  mode = "all",
  weights =
    E(network_graph)$weight
)

node_betweenness <- betweenness(
  network_graph,
  directed = FALSE,
  weights =
    E(network_graph)$distance,
  normalized = TRUE
)

node_harmonic <- harmonic_centrality(
  network_graph,
  weights =
    E(network_graph)$distance,
  normalized = TRUE
)

node_eigenvector <- eigen_centrality(
  network_graph,
  directed = FALSE,
  weights =
    E(network_graph)$weight,
  scale = TRUE
)$vector

node_pagerank <- page_rank(
  network_graph,
  directed = FALSE,
  weights =
    E(network_graph)$weight
)$vector

node_clustering <- transitivity(
  network_graph,
  type = "local",
  isolates = "zero"
)

topology_table <- data.frame(
  node_id = V(network_graph)$name,
  degree = as.numeric(node_degree),
  strength = as.numeric(node_strength),
  betweenness =
    as.numeric(node_betweenness),
  harmonic_centrality =
    as.numeric(node_harmonic),
  eigenvector_centrality =
    as.numeric(node_eigenvector),
  pagerank =
    as.numeric(node_pagerank),
  clustering_coefficient =
    as.numeric(node_clustering),
  stringsAsFactors = FALSE
)

#----------------------------------------------------------------------
# 24. Community detection
#----------------------------------------------------------------------

if (
  vcount(network_graph) > 2 &&
  ecount(network_graph) > 1
) {
  
  community_result <- cluster_louvain(
    network_graph,
    weights =
      E(network_graph)$weight
  )
  
  topology_table$community <-
    as.integer(
      membership(
        community_result
      )[
        topology_table$node_id
      ]
    )
  
} else {
  
  topology_table$community <- 1L
}

#----------------------------------------------------------------------
# 25. Join topology with node metadata
#----------------------------------------------------------------------

cytoscape_nodes <- cytoscape_nodes %>%
  left_join(
    topology_table,
    by = "node_id"
  ) %>%
  arrange(
    desc(degree),
    desc(strength)
  )

#----------------------------------------------------------------------
# 26. Network-level topology summary
#----------------------------------------------------------------------

component_results <- components(
  network_graph
)

network_summary <- data.frame(
  metric = c(
    "Number_of_nodes",
    "Number_of_edges",
    "Network_density",
    "Number_of_components",
    "Largest_component_size",
    "Global_clustering_coefficient",
    "Average_degree",
    "Average_strength",
    "Average_path_length",
    "Network_diameter",
    "Selected_omics_1_features",
    "Selected_omics_2_features",
    "Selected_omics_3_features",
    "Total_selected_features",
    "Possible_omics_feature_pairs",
    "Retained_CPITN_edges",
    "Retained_omics_edges"
  ),
  
  value = c(
    vcount(network_graph),
    
    ecount(network_graph),
    
    edge_density(
      network_graph,
      loops = FALSE
    ),
    
    component_results$no,
    
    max(
      component_results$csize
    ),
    
    transitivity(
      network_graph,
      type = "global",
      isolates = "zero"
    ),
    
    mean(
      degree(network_graph)
    ),
    
    mean(
      strength(
        network_graph,
        weights =
          E(network_graph)$weight
      )
    ),
    
    if (is_connected(network_graph)) {
      mean_distance(
        network_graph,
        directed = FALSE,
        weights =
          E(network_graph)$distance
      )
    } else {
      NA_real_
    },
    
    if (is_connected(network_graph)) {
      diameter(
        network_graph,
        directed = FALSE,
        weights =
          E(network_graph)$distance
      )
    } else {
      NA_real_
    },
    
    sum(
      candidate_features$omics ==
        "omics_1"
    ),
    
    sum(
      candidate_features$omics ==
        "omics_2"
    ),
    
    sum(
      candidate_features$omics ==
        "omics_3"
    ),
    
    nrow(
      candidate_features
    ),
    
    choose(
      nrow(candidate_features),
      2
    ),
    
    nrow(
      cpitn_edges
    ),
    
    nrow(
      omics_edges
    )
  ),
  
  stringsAsFactors = FALSE
)

print(network_summary)

#----------------------------------------------------------------------
# 27. Create hub table
#----------------------------------------------------------------------

hub_table <- cytoscape_nodes %>%
  filter(
    node_type != "Clinical"
  ) %>%
  arrange(
    desc(degree),
    desc(betweenness),
    desc(strength)
  ) %>%
  select(
    node_id,
    label,
    node_type,
    community,
    degree,
    strength,
    betweenness,
    harmonic_centrality,
    eigenvector_centrality,
    pagerank,
    clustering_coefficient,
    cpitn_rho,
    cpitn_p_value,
    cpitn_fdr_block,
    cpitn_fdr_global
  )

#----------------------------------------------------------------------
# 28. Export Cytoscape and topology files
#----------------------------------------------------------------------

write.csv(
  cytoscape_nodes,
  file.path(
    network_dir,
    "CPITN_Cytoscape_nodes.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  cytoscape_edges,
  file.path(
    network_dir,
    "CPITN_Cytoscape_edges.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  topology_table,
  file.path(
    network_dir,
    "CPITN_node_topology.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  network_summary,
  file.path(
    network_dir,
    "CPITN_network_summary.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  hub_table,
  file.path(
    network_dir,
    "CPITN_network_hubs.csv"
  ),
  row.names = FALSE,
  na = ""
)

if (
  nrow(omics_pairwise_results) > 0
) {
  
  write.csv(
    omics_pairwise_results,
    file.path(
      network_dir,
      paste0(
        "CPITN_selected_feature_",
        "pairwise_associations.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
}

#----------------------------------------------------------------------
# 29. Export GraphML
#
# GraphML preserves network and attribute information for Cytoscape.
#----------------------------------------------------------------------

write_graph(
  network_graph,
  file.path(
    network_dir,
    "CPITN_multiomics_network.graphml"
  ),
  format = "graphml"
)

#----------------------------------------------------------------------
# 30. Export Cytoscape SIF
#----------------------------------------------------------------------

sif_table <- cytoscape_edges %>%
  transmute(
    source = source,
    interaction = paste(
      edge_type,
      sign,
      sep = "_"
    ),
    target = target
  )

write.table(
  sif_table,
  file.path(
    network_dir,
    "CPITN_multiomics_network.sif"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#----------------------------------------------------------------------
# 31. Network visualization attributes
#----------------------------------------------------------------------

node_metadata_ordered <- cytoscape_nodes[
  match(
    V(network_graph)$name,
    cytoscape_nodes$node_id
  ),
]

V(network_graph)$size <-
  5 +
  3 * sqrt(
    degree(network_graph) + 1
  )

# Show CPITN plus nodes with degree at or above the
# 90th percentile
degree_label_threshold <- as.numeric(
  quantile(
    degree(network_graph),
    probs = 0.90,
    na.rm = TRUE
  )
)

V(network_graph)$label <- ifelse(
  V(network_graph)$name ==
    clinical_node_id |
    degree(network_graph) >=
    degree_label_threshold,
  
  node_metadata_ordered$label,
  
  NA
)

V(network_graph)$shape <- ifelse(
  node_metadata_ordered$node_type ==
    "Clinical",
  "square",
  "circle"
)

# Node colors by omics block
V(network_graph)$color <- case_when(
  node_metadata_ordered$node_type ==
    "Clinical" ~ "#E41A1C",
  
  node_metadata_ordered$node_type ==
    "omics_1" ~ "#4DAF4A",
  
  node_metadata_ordered$node_type ==
    "omics_2" ~ "#377EB8",
  
  node_metadata_ordered$node_type ==
    "omics_3" ~ "#984EA3",
  
  TRUE ~ "#999999"
)

E(network_graph)$width <-
  0.5 +
  4 * E(network_graph)$weight

E(network_graph)$lty <- ifelse(
  E(network_graph)$sign ==
    "Positive",
  1,
  2
)

E(network_graph)$color <- ifelse(
  E(network_graph)$sign ==
    "Positive",
  "#555555",
  "#AAAAAA"
)

set.seed(2026)

network_layout <- layout_with_fr(
  network_graph,
  weights =
    E(network_graph)$weight,
  niter = 2000
)

#----------------------------------------------------------------------
# 32. Save PDF network
#----------------------------------------------------------------------

pdf(
  file.path(
    network_dir,
    "CPITN_multiomics_network.pdf"
  ),
  width = 16,
  height = 13
)

plot(
  network_graph,
  layout = network_layout,
  vertex.label.cex = 0.60,
  vertex.label.dist = 0.5,
  vertex.frame.color = "grey30",
  edge.curved = 0.05,
  main = paste0(
    "CPITN-Centered Multi-Omics Network\n",
    "CPITN-feature: |rho| >= ",
    cpitn_rho_threshold,
    ", p < ",
    cpitn_p_threshold,
    "; omics-omics: |rho| >= ",
    omics_rho_threshold,
    ", p < ",
    omics_p_threshold
  )
)

legend(
  "topleft",
  legend = c(
    "CPITN",
    "omics_1",
    "omics_2",
    "omics_3",
    "Positive association",
    "Negative association"
  ),
  pch = c(
    15,
    16,
    16,
    16,
    NA,
    NA
  ),
  lty = c(
    NA,
    NA,
    NA,
    NA,
    1,
    2
  ),
  col = c(
    "#E41A1C",
    "#4DAF4A",
    "#377EB8",
    "#984EA3",
    "#555555",
    "#AAAAAA"
  ),
  bty = "n"
)

dev.off()

#----------------------------------------------------------------------
# 33. Save high-resolution PNG network
#----------------------------------------------------------------------

png(
  file.path(
    network_dir,
    "CPITN_multiomics_network.png"
  ),
  width = 4800,
  height = 3900,
  res = 300
)

plot(
  network_graph,
  layout = network_layout,
  vertex.label.cex = 0.60,
  vertex.label.dist = 0.5,
  vertex.frame.color = "grey30",
  edge.curved = 0.05,
  main = paste0(
    "CPITN-Centered Multi-Omics Network\n",
    "CPITN-feature: |rho| >= ",
    cpitn_rho_threshold,
    ", p < ",
    cpitn_p_threshold,
    "; omics-omics: |rho| >= ",
    omics_rho_threshold,
    ", p < ",
    omics_p_threshold
  )
)

legend(
  "topleft",
  legend = c(
    "CPITN",
    "omics_1",
    "omics_2",
    "omics_3",
    "Positive association",
    "Negative association"
  ),
  pch = c(
    15,
    16,
    16,
    16,
    NA,
    NA
  ),
  lty = c(
    NA,
    NA,
    NA,
    NA,
    1,
    2
  ),
  col = c(
    "#E41A1C",
    "#4DAF4A",
    "#377EB8",
    "#984EA3",
    "#555555",
    "#AAAAAA"
  ),
  bty = "n"
)

dev.off()

#----------------------------------------------------------------------
# 34. Save R objects for later reuse
#----------------------------------------------------------------------

save(
  cpitn,
  omics_1_mat,
  omics_2_mat,
  omics_3_mat,
  all_results,
  candidate_features,
  selected_feature_matrix,
  cpitn_edges,
  omics_edges,
  cytoscape_edges,
  cytoscape_nodes,
  topology_table,
  network_summary,
  hub_table,
  network_graph,
  file = file.path(
    network_dir,
    "CPITN_multiomics_network_objects.RData"
  )
)

#----------------------------------------------------------------------
# 35. Final report
#----------------------------------------------------------------------

cat(
  "\n====================================================\n"
)

cat(
  "Analysis completed.\n"
)

cat(
  "Output directory:\n",
  network_dir,
  "\n"
)

cat(
  "\nMain Cytoscape files:\n"
)

cat(
  "1. CPITN_multiomics_network.graphml\n"
)

cat(
  "2. CPITN_Cytoscape_nodes.csv\n"
)

cat(
  "3. CPITN_Cytoscape_edges.csv\n"
)

cat(
  "4. CPITN_multiomics_network.sif\n"
)

cat(
  "\nTopology files:\n"
)

cat(
  "1. CPITN_node_topology.csv\n"
)

cat(
  "2. CPITN_network_hubs.csv\n"
)

cat(
  "3. CPITN_network_summary.csv\n"
)

cat(
  "\nNetwork figures:\n"
)

cat(
  "1. CPITN_multiomics_network.pdf\n"
)

cat(
  "2. CPITN_multiomics_network.png\n"
)

cat(
  "====================================================\n"
)