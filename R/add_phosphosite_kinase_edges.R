#' Add Phosphosite-Kinase Edges to Kinase Interaction Network
#'
#' This function adds phosphosite-kinase edges to a given kinase interaction network (KIN).
#' It identifies the sources and targets for phosphosite of a kinase -> kinase edges and 
#' adds substrate -> protein edges to the network.
#'
#' @param KIN A data frame representing the kinase interaction network. It should contain 
#' columns 'Source', 'Target', and 'Target_Uniprot'.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#'
#' @return A data frame with the updated kinase interaction network, including the added 
#' phosphosite-kinase edges.
#' @export
add_phosphosite_kinase_edges <- function(
    KIN,
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL
) {
   
    KIN$Edge_type <- 'KS'
    # Get sources and targets for phosphosite of a kinase -> kinase edges 
    substrate_edges_sources <- KIN[which(
        KIN$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings[["ACC#"]], tyrosine_kinase_data$kinase_name_mappings[["ACC#"]])
        ),]$Target
    substrate_edges_targets <- KIN[which(
        KIN$Target_Uniprot %in% c(serine_threonine_kinase_data$kinase_name_mappings[["ACC#"]], tyrosine_kinase_data$kinase_name_mappings[["ACC#"]])
        ),]$Target_Uniprot 
    # Add substrate -> protein edges
    other_cols <- setdiff(colnames(KIN), c('Source', 'Target'))
    substrate_edges <- unique(data.frame(
        Source = substrate_edges_sources,
        Target = substrate_edges_targets
    ))
    if (nrow(substrate_edges) == 0) {
        return(KIN)
    }
    for (col in other_cols) {
        substrate_edges[[col]] <- NA
    }
    substrate_edges$Edge_type <- 'SK'
    KIN <- rbind(KIN, substrate_edges)

    return(KIN)
}