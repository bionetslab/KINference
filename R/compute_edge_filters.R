#' Compute CORR Filter
#'
#' This function computes a filter for edges based on correlation and adjusted p-values.
#' It removes edges where the correlation between the phosphorylation site of the kinase 
#' and the kinase is below a specified threshold (`delta`) and the adjusted p-value of 
#' the correlation is above a specified threshold (`epsilon`).
#'
#' @param KIN A data frame containing the kinase interaction network.
#' @param x A matrix or data frame containing the data for correlation computation.
#' @param delta A numeric value specifying the correlation threshold. Default is 0.8.
#' @param epsilon A numeric value specifying the p-value threshold. Default is 0.05.
#' @param multiple_testing_correction_method A character string specifying the method for 
#' multiple testing correction. Default is 'BH' (Benjamini-Hochberg).
#'
#' @return A data frame containing the filtered edges.
#' @noRd
compute_CORR_filter <- function(
    KIN,
    x,
    delta = 0.8,
    epsilon = 0.05,
    m = 10,
    multiple_testing_correction_method = 'BH',
    serine_threonine_kinase_data,
    tyrosine_kinase_data
) {

    # Computing mask of edges that should be removed
    # Only edges that can be removed are phosphorylated kinase -> phosphorylation site of a kinase edges
    # These edges are removed IFF the correlation between the phosphorylation site of the kinase and the kinase is below a threshold delta and the adjusted p-value of the correlation is above a threshold epsilon
    keep_edge_mask <- rep(TRUE, nrow(KIN))
    # don't consider kinase -> phosphorylation site of a kinase edges with less than m measurements
    not_enough_measurements_indices <- unname(which((ncol(x) - rowSums(is.na(x))) < m))

    # Iterate over all kinases
    for (kinase in unique(KIN$Source)) {

        # get all kinase phosphorylation sites: If there are no phosphorylation sites, skip the kinase
        source_phosphoSites <- unique(KIN[KIN$Target_Uniprot == kinase, ]$Target)
        if (length(source_phosphoSites) == 0) {
            next
        }

        # Get all kinase -> phosphorylation site of a kinase edges
        kinase_substrate_edge_indicies <- which(
            (KIN$Source == kinase) &
            (KIN$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings[["ACC#"]], tyrosine_kinase_data$kinase_name_mappings[["ACC#"]]))
        )

        # Default: Remove kinase -> phosphorylation site of a kinase edges
        # If a significant correlation is found, the edge is kept
        keep_edge_mask[kinase_substrate_edge_indicies] <- FALSE 
        
        for (kinase_substrate_edge_index in kinase_substrate_edge_indicies) {
            corr_target <- unlist(x[setdiff(which(rownames(x) == KIN$Target[kinase_substrate_edge_index]), not_enough_measurements_indices),])
            # If the corr target has not enough measurements -> Delete the edge
            if (length(corr_target) == 0) {
                next
            }
            
            p_values <- list()
            corr_values <- list()
            for (source_phosphoSite in source_phosphoSites) {
                corr_source <- unlist(x[setdiff(which(rownames(x) == source_phosphoSite), not_enough_measurements_indices),])
                # If the corr source has not enough measurements it is not considered
                if (length(corr_source) == 0) {
                    next
                }
                corr_test <- stats::cor.test(corr_source, corr_target)
                corr_values[[source_phosphoSite]] <- corr_test$estimate
                p_values[[source_phosphoSite]] <- corr_test$p.value
            }

            # Apply multiple testing correction
            adjusted_p_values <- stats::p.adjust(unname(p_values), method = multiple_testing_correction_method)
            names(adjusted_p_values) <- names(p_values)
            for (source_phosphoSite in names(adjusted_p_values)) {
                # Keep edge as soon as there is correlation found that surpasses the thresholds
                if ((adjusted_p_values[[source_phosphoSite]] <= epsilon) & (abs(corr_values[[source_phosphoSite]]) >= delta)) {
                    keep_edge_mask[kinase_substrate_edge_index] <- TRUE
                    break
                }
            }
        }
    }

    return(KIN[keep_edge_mask, ])
}

#' Compute PCST Filter
#'
#' This function computes a filter for edges based the Prize-Collecting Steiner Tree (PCST) algorithm.
#'
#' @param KIN A data frame containing the kinase interaction network.
#' @param gamma Gamma as defined in the paper. Default is 1.0.
#' @param log2FC_colname A string specifying the column name for log2 fold change values in the KIN data. Default is 'f'.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#'
#' @return A data frame containing the filtered edges.
#' @noRd
compute_PCST_filter <- function(
    KIN,
    gamma = 1.0,
    log2FC_colname = 'f',
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL
) {
    numpy <- reticulate::import("numpy", delay_load = TRUE, convert = FALSE)
    pcst_fast <- reticulate::import("pcst_fast", delay_load = TRUE)

    KIN <- Kinference::add_phosphosite_kinase_edges(KIN, serine_threonine_kinase_data, tyrosine_kinase_data)
    KIN[KIN$Edge_type == 'SK', log2FC_colname] <- gamma

    node_names <- unique(c(KIN$Source, KIN$Target))
    nodes <- 0:(length(node_names)-1)
    names(nodes) <- node_names

    KIN$Source_idx <- unname(sapply(KIN$Source, function(x) nodes[[x]]))
    KIN$Target_idx <- unname(sapply(KIN$Target, function(x) nodes[[x]]))

    edges <- numpy$array(KIN[, c('Source_idx', 'Target_idx')], dtype=numpy$int64)
    node_prizes <- numpy$array(abs(KIN[[log2FC_colname]]), dtype=numpy$float64)
    edge_costs <- numpy$full(shape=nrow(KIN), fill_value=gamma, dtype=numpy$float64)

    # defining pcst algorithm parameters
    root <- -1L
    num_clusters <- 1L
    pruning <- 'strong'
    verbosity_level <- 0L

    result <- pcst_fast$pcst_fast(
        edges,
        node_prizes,
        edge_costs,
        root,
        num_clusters,
        pruning,
        verbosity_level
    )

    return(KIN[result[[2]], ])
}


#' Compute Edge Filters for KIN Data
#'
#' This function applies edge filtering techniques to the provided KIN data.
#' It supports two types of filters: CORR filter and PCST filter.
#'
#' @param data A list containing data frames with the data to be filtered.
#' @param KIN A data frame representing the KIN data.
#' @param output.path A string specifying the directory where the results will be saved. Default is 'results'.
#' @param output.id A string specifying the identifier for the output files. Default is 'key'.
#' @param apply.CORR A logical value indicating whether to apply the CORR filter. Default is FALSE.
#' @param apply.PCST A logical value indicating whether to apply the PCST filter. Default is TRUE.
#' @param gamma Gamma as defined in the paper. Default is 1.0.
#' @param delta Delta as defined in the paper. Default is 0.8.
#' @param epsilon Epsilon as defined in the paper. Default is 0.05.
#' @param m Minimum measurement threshold. Default is 10.
#' @param multiple_testing_correction A string specifying the method for multiple testing correction in the CORR filter. Default is 'BH'.
#' @param log2FC_colname A string specifying the column name for log2 fold change values in the KIN data. Default is 'f'.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#' @param translate_uniprots A boolean value indicating if the Uniprot IDs of other species should be translated to human Uniprot IDs. Default is TRUE.
#'
#' @export
compute_edge_filters <- function(
    data,
    KIN,
    output.path = 'results',
    output.id = 'key',
    apply.CORR = FALSE,
    apply.PCST = TRUE,
    gamma = 1.0,
    delta = 0.8,
    epsilon = 0.05,
    m = 10,
    multiple_testing_correction = 'BH',
    log2FC_colname = 'f', 
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL,
    translate_uniprots = TRUE
) {

    dir.create(paste0(output.path, '/edge_filtered_KIN/'), recursive = TRUE)

    if (apply.CORR) {
        message("Computing CORR filter!")

        CORR_KINs <- list()
        CORR_KINs[["x0"]] <- compute_CORR_filter(
            KIN,
            data[["x0"]], 
            delta = delta, 
            epsilon = epsilon, 
            m = m,
            multiple_testing_correction_method = multiple_testing_correction,
            serine_threonine_kinase_data = serine_threonine_kinase_data,
            tyrosine_kinase_data = tyrosine_kinase_data
        )
        CORR_KINs[["x1"]] <- compute_CORR_filter(
            KIN,
            data[["x1"]], 
            delta = delta, 
            epsilon = epsilon, 
            m = m,
            multiple_testing_correction_method = multiple_testing_correction,
            serine_threonine_kinase_data = serine_threonine_kinase_data,
            tyrosine_kinase_data = tyrosine_kinase_data
        )
        if (translate_uniprots) {
            readr::write_tsv(
                Kinference::add_phosphosite_kinase_edges(CORR_KINs[['x0']], serine_threonine_kinase_data, tyrosine_kinase_data), 
                paste0(output.path, '/edge_filtered_KIN/', output.id, '_CORR_filtered_x0_KIN.tsv')
            )
            readr::write_tsv(
                Kinference::add_phosphosite_kinase_edges(CORR_KINs[['x1']], serine_threonine_kinase_data, tyrosine_kinase_data), 
                paste0(output.path, '/edge_filtered_KIN/', output.id, '_CORR_filtered_x1_KIN.tsv')
            )
        } else {
            readr::write_tsv(
                CORR_KINs[['x0']], 
                paste0(output.path, '/edge_filtered_KIN/', output.id, '_CORR_filtered_x0_KIN.tsv')
            )
            readr::write_tsv(
                CORR_KINs[['x1']], 
                paste0(output.path, '/edge_filtered_KIN/', output.id, '_CORR_filtered_x1_KIN.tsv')
            )
        }
    }    

    if (apply.PCST) {
        message("Computing PCST filter!")
        
        PCST_KIN <- compute_PCST_filter(
            KIN,
            gamma = gamma,
            log2FC_colname = log2FC_colname,
            serine_threonine_kinase_data = serine_threonine_kinase_data,
            tyrosine_kinase_data = tyrosine_kinase_data
        )

        dir.create(paste0(output.path, '/edge_filtered_KIN/'), recursive = TRUE)
        readr::write_tsv(
            PCST_KIN, 
            paste0(output.path, '/edge_filtered_KIN/', output.id, '_PCST_filtered_KIN.tsv')
        )
    }    
    return()
}