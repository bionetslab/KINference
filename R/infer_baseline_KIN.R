#' Infer Baseline Kinase Inference
#'
#' This function infers the baseline kinase activity for Serine/Threonine and Tyrosine kinases.
#'
#' @param md Metadata or input data required for the inference.
#' @param output.path The path where the results will be saved. Default is 'results'.
#' @param output.id An identifier for the output. Default is 'key'.
#' @param n Number of maximum edges to keep per kinases -> phosphorylation site interaction computation. Default is 15.
#' @param alpha Minimum percentile score to keep kinase -> phosphorylation site interaction. Default is 0.9.
#' @param left_phoSite_flank Number of amino acids to the left of the phosphorylation site. Default is 5. (The underlying kinase data has to be extended to increase this value)
#' @param right_phoSite_flank Number of amino acids to the right of the phosphorylation site. Default is 4. (The underlying kinase data has to be extended to increase this value)
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#' @param translate_uniprots A boolean value indicating if the Uniprot IDs of other species should be translated to human Uniprot IDs. Default is TRUE.
#'
#' @return A list containing the inferred baseline kinase activity for Serine/Threonine and Tyrosine kinases.
#' @export
infer_baseline_KIN <- function(
    md,
    output.path = 'results',
    output.id = 'key',
    n = 15,
    alpha = 0.9,
    left_phoSite_flank = 5,
    right_phoSite_flank = 4,
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL,
    translate_uniprots = TRUE
) {

    baseline_KIN <- list()

    baseline_KIN[['serine_threonine']] <- infer_kinase_specific_KIN(
        md = md,
        kinase_type = 'serine_threonine',
        n = n,
        alpha = alpha,
        left_phoSite_flank = left_phoSite_flank,
        right_phoSite_flank = right_phoSite_flank,
        serine_threonine_kinase_data = serine_threonine_kinase_data,
        tyrosine_kinase_data = tyrosine_kinase_data
    )

    baseline_KIN[['tyrosine']] <- infer_kinase_specific_KIN(
        md = md,
        kinase_type = 'tyrosine',
        n = n,
        alpha = alpha,
        left_phoSite_flank = left_phoSite_flank,
        right_phoSite_flank = right_phoSite_flank,
        serine_threonine_kinase_data = serine_threonine_kinase_data,
        tyrosine_kinase_data = tyrosine_kinase_data
    )

    baseline_KIN[['combined']] <- dplyr::bind_rows(baseline_KIN) 

    # Save results to TSV files
    dir.create(paste0(output.path, '/baseline_KIN/'), recursive = TRUE)
    
    if (translate_uniprots) {
        readr::write_tsv(
            Kinference::add_phosphosite_kinase_edges(baseline_KIN[['serine_threonine']], serine_threonine_kinase_data, tyrosine_kinase_data), 
            paste0(output.path, '/baseline_KIN/', output.id, '_serine_threonine_baseline_KIN.tsv')
        )
        readr::write_tsv(
            Kinference::add_phosphosite_kinase_edges(baseline_KIN[['tyrosine']], serine_threonine_kinase_data, tyrosine_kinase_data), 
            paste0(output.path, '/baseline_KIN/', output.id, '_tyrosine_baseline_KIN.tsv')
        )
        readr::write_tsv(
            Kinference::add_phosphosite_kinase_edges(baseline_KIN[['combined']], serine_threonine_kinase_data, tyrosine_kinase_data), 
            paste0(output.path, '/baseline_KIN/', output.id, '_KIN.tsv')
        )
    } else {
        readr::write_tsv(
            baseline_KIN[['serine_threonine']], 
            paste0(output.path, '/baseline_KIN/', output.id, '_serine_threonine_baseline_KIN.tsv')
        )
        readr::write_tsv(
            baseline_KIN[['tyrosine']], 
            paste0(output.path, '/baseline_KIN/', output.id, '_tyrosine_baseline_KIN.tsv')
        )
        readr::write_tsv(
            baseline_KIN[['combined']], 
            paste0(output.path, '/baseline_KIN/', output.id, '_KIN.tsv')
        )
    }
    return(baseline_KIN)
}


#' Query UniProt Web Service for Protein Sequences
#'
#' This function queries the UniProt web service for protein sequences based on 
#' a vector of UniProt accession numbers.
#'
#' @param uniprots A character vector of UniProt accession numbers to query.
#' @param fields A character vector of fields to retrieve from the UniProt web service. 
#'        Default is c("accession", "sequence").
#'
#' @return A data.frame containing the queried information from UniProt, with an 
#'         additional column 'query' that matches the input accession numbers.
#' @noRd
.query_uniprot_ws_for_sequence <- function(
    uniprots, 
    fields = c("accession","sequence")
) {
    queries <- paste0("accession:", uniprots)

    res <- UniProt.ws::queryUniProt(queries, fields = fields)
    res <- res %>%
        dplyr::mutate(query = Entry) %>%
        dplyr::right_join(data.frame(Entry = uniprots), by = "Entry")
    return(res)
}

#' Retrieve UniProt Sequences from Web
#'
#' This function retrieves protein sequences from the UniProt database for a given set of UniProt IDs.
#'
#' @param uniprotIDs A character vector of UniProt IDs for which sequences are to be retrieved.
#' @param chunkSize An integer specifying the number of UniProt IDs to query in each batch. Default is 25.
#' @param fields A character vector specifying the fields to retrieve from UniProt. Default is c("accession", "sequence").
#'
#' @return A data.frame containing the retrieved information for the specified UniProt IDs.
#' @noRd
uniprotSequencesFromWeb <- function(
    uniprotIDs, 
    chunkSize = 25, 
    fields = c("accession","sequence")
    )
{
    uniqueUniprots <- unique(uniprotIDs)
    chunks <- split(uniqueUniprots, (1:length(uniqueUniprots))%/%chunkSize)
    infoMapList <- pbapply::pblapply(chunks, .query_uniprot_ws_for_sequence, fields = fields)
    infoMap <- dplyr::bind_rows(infoMapList)
    return (infoMap)
}


#' Infer Kinase-Specific KIN
#'
#' This function infers kinase-specific kinase interaction networks (KIN) based on the provided metadata and kinase type.
#'
#' @param md A data frame containing metadata with columns 'Uniprot', 'Pos', 'Protein', and 'f'.
#' @param kinase_type A character string specifying the type of kinase. Must be either 'Serine_Threonine' or 'Tyrosine'.
#' @param n An integer specifying the maximum number of edges to retain in the network. Default is 15.
#' @param alpha A numeric value specifying the threshold for percentile scores. Default is 0.9.
#' @param left_phoSite_flank An integer specifying the number of amino acids to the left of the phosphorylation site. Default is 5. (The underlying kinase data has to be extended to increase this value)
#' @param right_phoSite_flank An integer specifying the number of amino acids to the right of the phosphorylation site. Default is 4. (The underlying kinase data has to be extended to increase this value)
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#'
#' @return A data frame containing the inferred kinase interaction network (KIN) with columns 'Source', 'Target', 'Target_Uniprot', 'ModifiedSequence', 'log2Score', 'percentileScore', 'percentileRank', 'f', and 'Type'.
#'
#' @details
#' The function performs the following steps:
#' 1. Loads kinase data based on the specified kinase type.
#' 2. Retrieves sequences from Uniprot based on the metadata.
#' 3. Extracts sequences surrounding phosphorylation sites.
#' 4. Computes log2 scores for the sequences using kinase motifs.
#' 5. Retrieves percentile scores for the sequences.
#' 6. Constructs the final data frame containing the inferred KIN.
#' 7. Saves the results to a TSV file.
#' @noRd
infer_kinase_specific_KIN <- function(
    md,
    kinase_type,
    n = 15,
    alpha = 0.9,
    left_phoSite_flank = 5,
    right_phoSite_flank = 4,
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL
) {
    
    if (kinase_type == 'serine_threonine') {
        kinase_data <- serine_threonine_kinase_data
        phosphorylated_aa <- c('S', 'T')
    } else if (kinase_type == 'tyrosine') {
        kinase_data <- tyrosine_kinase_data
        phosphorylated_aa <- c('Y')
    } else {
        stop("Invalid kinase type. Please choose either 'Serine_Threonine' or 'Tyrosine'.")
    }

    kinases <- kinase_data$kinase_name_mappings$`ACC#`
    # Loading sequences from Uniprot  
    sequences <- uniprotSequencesFromWeb(unique(md$Uniprot))

    get_modified_sequence <- function(seq, pos) {
        if (pos > nchar(seq) || !(stringr::str_sub(seq, pos, pos) %in% phosphorylated_aa)) {
            return(NULL)
        } else if ((pos > left_phoSite_flank) && ((pos + right_phoSite_flank) < nchar(seq))) {
            return(list(seq = stringr::str_sub(seq, pos - left_phoSite_flank, pos + right_phoSite_flank), zero_idx = left_phoSite_flank + 1))
        } else if (pos > left_phoSite_flank) {
            return(list(seq = stringr::str_sub(seq, pos - left_phoSite_flank), zero_idx = left_phoSite_flank + 1))
        } else if ((pos + right_phoSite_flank) < nchar(seq)) {
            n_seq <- stringr::str_sub(seq, 1, pos + right_phoSite_flank)
            return(list(seq = n_seq, zero_idx = nchar(n_seq) - right_phoSite_flank))
        } else {
            return(NULL)
        }
    }

    # Extracting sequences surrounding phosphorylation sites
    sequences_list <- lapply(1:nrow(md), function(i) {
        seq <- sequences$Sequence[which(sequences$query == md$Uniprot[i])]
        pos <- md$Pos[i]
        get_modified_sequence(seq, pos)
    })

    # Removing NULL values
    valid_indices <- which(!sapply(sequences_list, is.null))
    sequences_list <- sequences_list[valid_indices]
    md_filtered <- md[valid_indices,]

    # Implements the substrate scoring from Johnson et al. 2023 (https://doi.org/10.1038/s41586-022-05575-3)
    get_scores_log2 <- function(seq_info, kinase_data) {
        seq <- seq_info$seq
        zero_idx <- seq_info$zero_idx
        seq_split <- stringr::str_extract_all(seq, stringr::boundary('character'))[[1]]
        indices <- if (zero_idx == (left_phoSite_flank + 1)) {
            seq(-left_phoSite_flank, length(seq_split) - (left_phoSite_flank + 1), 1)
        } else {
            seq(-(length(seq_split) - left_phoSite_flank), right_phoSite_flank, 1)
        }
    
        aa <- paste0(indices, seq_split)[-zero_idx]
        scores_log2 <- log2(apply(kinase_data$kinase_motifs %>% dplyr::select(dplyr::all_of(aa)), 1, prod))
        names(scores_log2) <- kinases
        scores_log2
    }
    # Compute scores for all sequences
    scores_log2_list <- lapply(sequences_list, get_scores_log2, kinase_data = kinase_data)
  
    # Retrieves the percentile scores
    percentile_scores_list <- lapply(scores_log2_list, function(scores_log2) {
        sapply(names(kinase_data$kinase_aprior_distributions), function(x) { 
            kinase_data$kinase_aprior_distributions[[x]](scores_log2[[x]]) 
        })
    })

    # Returns the baseline KIN
    construct_final_dataframe <- function(i, scores_log2, percentile_scores) {
        percentile_scores <- percentile_scores[order(-percentile_scores)]
        edges_threshold <- min(n, length(which(percentile_scores >= alpha)))
        if (edges_threshold == 0) {
            return(NULL)
        }
        seq_info <- sequences_list[[i]]
        new_rows <- data.frame(
            Source = names(percentile_scores)[1:edges_threshold], 
            Target = rep(md_filtered$Protein[i], edges_threshold), 
            Target_Uniprot = rep(md_filtered$Uniprot[i], edges_threshold),
            ModifiedSequence = rep(seq_info$seq, edges_threshold), 
            log2Score = scores_log2[names(percentile_scores)[1:edges_threshold]], 
            percentileScore = percentile_scores[1:edges_threshold], 
            percentileRank = seq(1, edges_threshold, 1), 
            f = rep(md_filtered$f[i], edges_threshold)
        )
        return(new_rows)
    }

    # Generate baseline KIN
    result_list <- mapply(
        construct_final_dataframe, 
        i = seq_along(scores_log2_list), 
        scores_log2 = scores_log2_list, 
        percentile_scores = percentile_scores_list, 
        SIMPLIFY = FALSE
    )
    # Remove null rows
    result_list <- result_list[which(!sapply(result_list, is.null))]
    baseline_KIN <- dplyr::bind_rows(result_list)
    rownames(baseline_KIN) <- NULL
    return(baseline_KIN)
}
