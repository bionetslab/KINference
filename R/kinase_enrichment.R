#' Compute Kinase Enrichment
#'
#' This function computes the enrichment of serine/threonine and tyrosine kinases
#' based on the provided metadata and baseline kinase data. The results are saved
#' as TSV files in the specified output path.
#'
#' @param md A data frame containing the metadata for the analysis.
#' @param baseline_KINs A list containing baseline kinase data for 'serine_threonine' and 'tyrosine'.
#' @param output.path A character string specifying the directory where the results will be saved. Default is 'results'.
#' @param output.id A character string used as a prefix for the output file names. Default is 'key'.
#' @param gamma A numeric value used as a parameter in the enrichment computation. Default is 1.0.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#'
#' @return A list containing the enrichment results for 'serine_threonine' and 'tyrosine' kinases.
#' @export
compute_kinase_enrichment <- function(
    md, 
    baseline_KINs, 
    output.path = 'results', 
    output.id = 'key',
    gamma = 1.0,
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL
) {
    
    kinase_enrichments <- list()

    kinase_enrichments[['serine_threonine']] <- compute_enrichment(
        md = md, 
        baseline_KIN = baseline_KINs[['serine_threonine']],
        kinase_type = 'serine_threonine',
        gamma = gamma,
        serine_threonine_kinase_data = serine_threonine_kinase_data,
        tyrosine_kinase_data = tyrosine_kinase_data
    )

    kinase_enrichments[['tyrosine']] <- compute_enrichment(
        md = md, 
        baseline_KIN = baseline_KINs[['tyrosine']],
        kinase_type = 'tyrosine',
        gamma = gamma,
        serine_threonine_kinase_data = serine_threonine_kinase_data,
        tyrosine_kinase_data = tyrosine_kinase_data
    )

    # Save results to TSV files
    dir.create(paste0(output.path, '/enrichments/'), recursive = TRUE)
    readr::write_tsv(
        kinase_enrichments[['serine_threonine']], 
        paste0(output.path, '/enrichments/', output.id, '_serine_threonine_kinase_enrichment.tsv')
    )
    readr::write_tsv(
        kinase_enrichments[['tyrosine']], 
        paste0(output.path, '/enrichments/', output.id, '_tyrosine_kinase_enrichment.tsv')
    )

    return(kinase_enrichments)
}

# Compute enrichment
#' This function computes the enrichment of kinases in upregulated and downregulated phosphorylation sites.
#'
#' @param md A data frame containing phosphorylation site data with at least columns 'f' (log2 fold change) and 'Protein'.
#' @param baseline_KIN A data frame containing kinase-substrate interactions with at least columns 'Source' (kinase) and 'Target' (substrate).
#' @param kinase_type A string specifying the type of kinase. Must be either 'serine_threonine' or 'tyrosine'.
#' @param gamma A numeric value specifying the threshold for log2 fold change to classify upregulated and downregulated sites. Default is 1.0.
#' @param serine_threonine_kinase_data An R object containing serine/threonine kinase data. Default is NULL.
#' @param tyrosine_kinase_data An R object containing tyrosine kinase data. Default is NULL.
#'
#' @return Data.frame containing the enrichment results
#' @noRd
compute_enrichment <- function(
    md, 
    baseline_KIN,
    kinase_type,
    gamma = 1.0,
    serine_threonine_kinase_data = NULL,
    tyrosine_kinase_data = NULL
) {
  
    if (kinase_type == 'serine_threonine') {
        md <- md %>% dplyr::filter(AA %in% c('S', 'T'))
        kinase_data <- serine_threonine_kinase_data
    } else if (kinase_type == 'tyrosine') {
        md <- md %>% dplyr::filter(AA == 'Y')
        kinase_data <- tyrosine_kinase_data
    } else {
        stop("Invalid kinase type. Please choose either 'Serine_Threonine' or 'Tyrosine'.")
    } 

    kinases <- kinase_data$kinase_name_mappings$`ACC#`
  
    # Upregulated set := Phosphosites with log2FC > gamma
    # Downregulated set := Phosphosites with log2FC < - (gamma)
    # Background set := Phosphosites with -(gamma) <= log2FC <= gamma
    upregulated_phosphorylationSites <- md[which(md$f > gamma),]
    downregulated_phosphorylationSites <- md[which(md$f < -gamma),]
    background_phosphorylationSites <- md[which(abs(md$f) <= gamma),]

    # Getting all interactions in each set 
    upregulated_interactions <- baseline_KIN[which(baseline_KIN$Target %in% upregulated_phosphorylationSites$Protein),]
    downregulated_interactions <- baseline_KIN[which(baseline_KIN$Target %in% downregulated_phosphorylationSites$Protein),]
    background_interactions <- baseline_KIN[which(baseline_KIN$Target %in% background_phosphorylationSites$Protein),]

    # Counting the occurrences of each kinase in the three sets
    up_counts <- data.frame(table(upregulated_interactions$Source))
    if (nrow(up_counts) != length(kinases)) up_counts <- rbind(up_counts, data.frame(Var1=setdiff(kinases, up_counts$Var1), Freq=0))
    down_counts <- data.frame(table(downregulated_interactions$Source))
    if (nrow(down_counts) != length(kinases)) down_counts <- rbind(down_counts, data.frame(Var1=setdiff(kinases, down_counts$Var1), Freq=0))
    background_counts <- data.frame(table(background_interactions$Source))
    if (nrow(background_counts) != length(kinases)) background_counts <- rbind(background_counts, data.frame(Var1=setdiff(kinases, background_counts$Var1), Freq=0))

    # Enrichment scoring
    enrich_val_up <- c()
    enrich_val_down <- c()
    dominant_enrichment_value <- c()
    dominant_enrichment_log2_value <- c()
    p_val_up <- c()
    p_val_down <- c()
    dominant_p_value <- c()
    dominant_directions <- c()

    up_set_size <- nrow(upregulated_phosphorylationSites)
    down_set_size <- nrow(downregulated_phosphorylationSites)
    background_set_size <- nrow(background_phosphorylationSites)
    for (i in 1:nrow(up_counts)) {

        up_hits <- up_counts$Freq[i]
        down_hits <- down_counts$Freq[i]
        background_hits <- background_counts$Freq[i]
        
        #TODO: Fix Haldane correction? Adding 0.5 to every cell of the contigency table in R does not seem to work for fisher.test
        # Haldane correction up_hits (No kinase hits in the upregulated (or background) interactions)
        if (up_hits == 0 | background_hits == 0) {
            d_up <- data.frame(
                x = c(up_hits + 0.5, background_hits + 0.5),
                y = c(up_set_size + 0.5, background_set_size + 0.5)
            )
        } else {
          d_up <- data.frame(
              x = c(up_hits, background_hits),
              y = c(up_set_size, background_set_size)
          )
        }
      
        # Haldane correction down_hits
        if (down_hits == 0 | background_hits == 0) {
            d_down <- data.frame(
                x = c(down_hits + 0.5, background_hits + 0.5),
                y = c(down_set_size + 0.5, background_set_size + 0.5)
            )
        } else {
            d_down <- data.frame(
                x = c(down_hits, background_hits),
                y = c(down_set_size, background_set_size)
            )
        }
        
        # Calculating p-value and enrichment score ((up_hits/up_set_size)/(background_hits/background_set_size))
        # using one sided fisher exact z test (assumption: background set is bigger than up-,downregulatory set)
        # Dominant value: max(enrich_up, enrich_down)
        test_up <- stats::fisher.test(d_up, alternative = 'greater')
        test_down <- stats::fisher.test(d_down, alternative = 'greater')
        
        enrich_up <- unname(test_up[['estimate']])
        enrich_down <- unname(test_down[['estimate']])
        
        enrich_val_up <- c(enrich_val_up, enrich_up)
        enrich_val_down <- c(enrich_val_down, enrich_down)
        dominant_value <- ifelse(enrich_up > enrich_down, 1, 2)
        dominant_directions <- c(dominant_directions, c('upregulated set', 'downregulated set')[dominant_value])
        dominant_enrichment_log2_value <- c(dominant_enrichment_log2_value, ifelse(enrich_up > enrich_down, log2(enrich_up), -log2(enrich_down)))
        dominant_enrichment_value <- c(dominant_enrichment_value, ifelse(enrich_up > enrich_down, enrich_up, 2 ** -log2(enrich_down)))
        p_val_up <- c(p_val_up, test_up[['p.value']])
        p_val_down <- c(p_val_down, test_down[['p.value']])
        dominant_p_value <- c(dominant_p_value, c(test_up[['p.value']], test_down[['p.value']])[dominant_value])
    }

    # Adjusting p values using Benjamini Hochberg correction
    p_val_up_adjusted <- p.adjust(p_val_up, method = 'BH')
    p_val_down_adjusted <- p.adjust(p_val_down, method = 'BH')
    dominant_p_value_adjusted <- p.adjust(dominant_p_value, method = 'BH')

    enrichment <- data.frame(uniprotID = up_counts$Var1, 
        gene = sapply(up_counts$Var1, function(x) {kinase_data$kinase_name_mappings$GENE[which(kinase_data$kinase_name_mappings$`ACC#` == x)]}),
        upregulated_set_hits = up_counts$Freq, 
        upregulated_set_size = up_set_size, 
        downregulated_set_hits = down_counts$Freq, 
        down_regulated_set_size = down_set_size, 
        background_set_hits = background_counts$Freq, 
        background_set_size = background_set_size, 
        upregulated_enrichment_value = enrich_val_up, 
        upregulated_enrichment_value_log2 = log2(enrich_val_up), 
        upregulated_p_value = p_val_up, 
        upregulated_adjusted_p_value = p_val_up_adjusted, 
        upregulated_p_value_log10_abs = abs(log10(p_val_up)), 
        upregulated_adjusted_p_value_log10_abs = abs(log10(p_val_up_adjusted)), 
        downregulated_enrichment_value = enrich_val_down,
        downregulated_enrichment_value_log2 = log2(enrich_val_down), 
        downregulated_p_value = p_val_down,
        downregulated_adjusted_p_value = p_val_down_adjusted, 
        downregulated_p_value_log10_abs = abs(log10(p_val_down)), 
        downregulated_adjusted_p_value_log10_abs = abs(log10(p_val_down_adjusted)),
        dominant_enrichment_value = dominant_enrichment_value, 
        dominant_enrichment_value_log2 = dominant_enrichment_log2_value,
        dominant_p_value = dominant_p_value, 
        dominant_adjusted_p_value = dominant_p_value_adjusted,
        dominant_p_value_log10_abs = abs(log10(dominant_p_value)), 
        dominant_adjusted_p_value_log10_abs = abs(log10(dominant_p_value_adjusted)),
        dominant_direction = dominant_directions
    )

    return(enrichment)
}