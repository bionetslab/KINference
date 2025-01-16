#' Compute DIFF Filter
#'
#' This function computes the DIFF filter for serine/threonine and tyrosine kinases
#' based on the provided baseline kinase data. The results are saved as TSV files.
#'
#' @param KIN A kinase interaction network with a column 'log2fc_colname' containing the log2 fold change values.
#' @param output.path A character string specifying the directory where the results will be saved. Default is 'results'.
#' @param output.id A character string specifying the identifier for the output files. Default is 'key'.
#' @param gamma A numeric value specifying the threshold for filtering. Default is 1.0.
#' @param log2FC_colname A character string specifying the column name containing the log2 fold change values. Default is 'f'.
#'
#' @return A list containing the filtered kinase data for 'serine_threonine' and 'tyrosine'.
#' @noRd
compute_DIFF_filter <- function(
    KIN,
    output.path = 'results',
    output.id = 'key',
    gamma = 1.0,
    log2FC_colname = 'f'
) {
        
    DIFF_filtered_KIN <- KIN[which(abs(KIN[[log2FC_colname]]) >= gamma),]

    return(DIFF_filtered_KIN)
}


#' Compute Functional Score Filter
#'
#' This function computes the functional score filter for kinase inference.
#'
#' @param md A data frame containing the metadata with columns 'Uniprot', 'Pos', and 'AA'.
#' @param KIN A kinase interaction network.
#' @param output.path A character string specifying the output directory path. Default is 'results'.
#' @param output.id A character string specifying the output file identifier. Default is 'key'.
#' @param beta A numeric value specifying the threshold for the functional score filter. Default is 0.4.
#'
#' @return A list containing filtered kinase information based on the functional score.
#' @noRd
compute_FS_filter <- function(
    md,
    KIN,
    output.path = 'results',
    output.id = 'key',
    beta = 0.4
) {
    if (file.exists(paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/FS_scores.tsv'))) {
        all_scores <- readr::read_tsv(
            paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/FS_scores.tsv'),
            show_col_types = FALSE
        )
    } else {
        # Compute functional score using FunSCoR
        message("Computing FS filter! This may take a while...")
        data_targets <- data.frame(
            acc=md$Uniprot,
            position=md$Pos,
            residue=md$AA
        )

        merged_phosphoproteome <- unique(rbind(data_targets, Kinference::phosphoproteome))
        # removing duplicates (keeping reads from lab data)
        merged_phosphoproteome <- merged_phosphoproteome[!duplicated(merged_phosphoproteome[, c(1,2)]), ]
        annotated_phos <- funscoR::annotate_sites(merged_phosphoproteome)
        ST_features <- funscoR::preprocess_features(annotated_phos, "ST")
        Y_features <- funscoR::preprocess_features(annotated_phos, "Y")

        ## train new model
        ST_model <- funscoR::train_funscore(ST_features, "ST", funscoR::psp, ncores = 8)
        Y_model <- funscoR::train_funscore(Y_features, "Y", funscoR::psp, ncores = 8)
    
        ## predict funcscoR for all sites
        ST_scores <- funscoR::predict_funscore(ST_features, ST_model, ncores = 8)
        Y_scores <- funscoR::predict_funscore(Y_features, Y_model, ncores = 8)
    
        ## gather all predictions
        all_scores <- dplyr::bind_rows(ST_scores, Y_scores) %>% 
            dplyr::mutate(probabilities = funscoR::log_scaling(probabilities))
        
        readr::write_tsv(
            all_scores, 
            paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/FS_scores.tsv')
        )
    }

    # Add functional score to baseline KINs
    KIN[['FS']] <- sapply(KIN[['Target']], function(x) {
        phoSite <- paste0(strsplit(x, "_")[[1]][1], "_", stringr::str_sub(strsplit(x, "_")[[1]][2], 2))
        if (phoSite %in% all_scores$sites) {
            return(all_scores$probabilities[which(all_scores$sites == phoSite)])
        } else {
            return(0)
        }
    })

    # Filter based on functional score
    FS_filtered_KIN <- KIN[which(KIN[['FS']] >= beta),]

    return(FS_filtered_KIN)
}


#' Compute Node Filters
#'
#' This function computes node filters based on the provided parameters and species.
#'
#' @param md Metadata required for the computation of FS filters.
#' @param KIN Kinase interaction network.
#' @param output.path The directory path where the output results will be saved. Default is 'results'.
#' @param output.id The identifier for the output files. Default is 'key'.
#' @param gamma The gamma parameter used in the computation of DIFF and combined node filters. Default is 1.0.
#' @param beta The beta parameter used in the computation of FS filters. Default is 0.4.
#' @param apply.DIFF A logical value specifying whether to apply the DIFF filter. Default is TRUE.
#' @param apply.FS A logical value specifying whether to apply the FS filter. Default is TRUE.
#' @param log2FC_colname The column name containing the log2 fold change values. Default is 'f'.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#' @param translate_uniprots A boolean value indicating if the Uniprot IDs of other species should be translated to human Uniprot IDs. Default is TRUE.
#'
#' @export
compute_node_filters <- function(
    md,
    KIN,
    output.path = 'results',
    output.id = 'key',
    gamma = 1.0,
    beta = 0.4,
    apply.DIFF = TRUE,
    apply.FS = TRUE,
    log2FC_colname = 'f',
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL,
    translate_uniprots = TRUE
) {

    dir.create(paste0(output.path, '/node_filtered_KIN/'), recursive = TRUE)

    if (apply.FS) {
        dir.create(paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/'), recursive = TRUE)
        FS_filtered_KIN <- compute_FS_filter(
            md = md,
            KIN = KIN,
            output.path = output.path,
            output.id = output.id,
            beta = beta
        )
        if (translate_uniprots) {
            readr::write_tsv(
                Kinference::add_phosphosite_kinase_edges(FS_filtered_KIN, serine_threonine_kinase_data, tyrosine_kinase_data), 
                paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/', output.id, '_KIN.tsv')
            )
        } else {
            readr::write_tsv(
                FS_filtered_KIN, 
                paste0(output.path, '/node_filtered_KIN/FS_filtered_KIN/', output.id, '_KIN.tsv')
            )
        }
    }

    if (apply.DIFF) {
        DIFF_filtered_KIN <- compute_DIFF_filter(
            KIN = KIN,
            log2FC_colname = log2FC_colname,
            output.path = output.path,
            output.id = output.id,
            gamma = gamma
        )
        dir.create(paste0(output.path, '/node_filtered_KIN/DIFF_filtered_KIN/'), recursive = TRUE)
        if (translate_uniprots) {
            readr::write_tsv(
                Kinference::add_phosphosite_kinase_edges(DIFF_filtered_KIN, serine_threonine_kinase_data, tyrosine_kinase_data), 
                paste0(output.path, '/node_filtered_KIN/DIFF_filtered_KIN/', output.id, '_KIN.tsv')
            )
        } else {
            readr::write_tsv(
                DIFF_filtered_KIN, 
                paste0(output.path, '/node_filtered_KIN/DIFF_filtered_KIN/', output.id, '_KIN.tsv')
            )
        }
    } 

    if (apply.DIFF & apply.FS) {
        KIN <- compute_DIFF_filter(
            KIN = FS_filtered_KIN,
            log2FC_colname = log2FC_colname,
            output.path = output.path,
            output.id = output.id,
            gamma = gamma
        )
    } else if (apply.DIFF) {
        KIN <- DIFF_filtered_KIN
    } else {
        KIN <- FS_filtered_KIN
    }

    if (translate_uniprots) {
        readr::write_tsv(
            Kinference::add_phosphosite_kinase_edges(KIN, serine_threonine_kinase_data, tyrosine_kinase_data), 
            paste0(output.path, '/node_filtered_KIN/', output.id, '_KIN.tsv')
        )
    } else {
        readr::write_tsv(
            KIN, 
            paste0(output.path, '/node_filtered_KIN/', output.id, '_KIN.tsv')
        )
    }
    

    return(KIN)
}