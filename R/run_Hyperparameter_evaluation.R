#' Run KINference hyperparameter evaluation
#'
#' This function performs the KINference analysis pipeline, which includes data preparation, 
#' baseline kinase inference, kinase enrichment computation, and application of various filters.
#'
#' @param x0.path Path to the first intensity matrix file. Default is NA.
#' @param x1.path Path to the second intensity matrix file. Default is NA.
#' @param f.path Path to the intensity vector file. Default is NA.
#' @param output.path Directory where the results will be saved. Default is 'results'.
#' @param output.id Identifier for the output files. Default is 'key'.
#' @param species Species name, either 'Homo sapiens' or 'Mus musculus'. Default is 'Homo sapiens'.
#' @param translate_uniprots If set to false, the Uniprot IDs will not be translated to human Uniprot IDs. If set to FALSE, the PCST filter will be automatically disabled because the PCST filter needs a connected network and, therefore, needs all Uniprot IDs to be human Uniprot IDs. Default is TRUE.
#' @param paired.samples Logical indicating if the samples are paired. Default is TRUE.
#' @param apply.log2 Logical indicating if log2 transformation should be applied to the data. Default is FALSE.
#' @param n Number of top kinases to infer. Default is 15.
#' @param alpha Parameter for baseline kinase inference. Default is 0.9.
#' @param apply.DIFF Logical indicating if the DIFF filter should be applied. Default is TRUE.
#' @param apply.FS Logical indicating if the FS filter should be applied. Default is TRUE.
#' @param apply.CORR Logical indicating if the CORR filter should be applied. Default is FALSE.
#' @param apply.PCST Logical indicating if the PCST filter should be applied. Default is TRUE.
#' @param beta Parameter for node filter computation. Default is 0.4.
#' @param gamma Parameter for kinase enrichment computation. Default is 1.0.
#' @param delta Parameter for edge filter computation. Default is 0.8.
#' @param epsilon Parameter for edge filter computation. Default is 0.05.
#' @param m Parameter for edge filter computation. Default is 10.
#' @param multiple_testing_correction Method for multiple testing correction. Default is 'BH'.
#' @param custom_serine_threonine_kinase_data.path Path to a custom serine/threonine kinase data file. Default is NULL.
#' @param custom_tyrosine_kinase_data.path Path to a custom tyrosine kinase data file. Default is NULL.
#'
#' @return None. The results are saved to the specified output path.
#' @export
run_Hyperparameter_evaluation <- function(
    x0.path = NA,
    x1.path = NA,
    f.path = NA,
    output.path = 'results',
    output.id = 'key',
    species = 'Homo sapiens',
    paired.samples = F,
    apply.log2 = F,
    n = 15,
    alpha = 0.9,
    apply.DIFF = TRUE,
    apply.FS = TRUE,
    apply.CORR = FALSE,
    apply.PCST = TRUE,
    beta = 0.4,
    gamma = 1.0,
    delta = 0.8,
    epsilon = 0.05,
    m = 10,
    multiple_testing_correction = 'BH',
    custom_serine_threonine_kinase_data.path = NULL,
    custom_tyrosine_kinase_data.path = NULL
) {

    # Checking validity of input arguments
    if ((!(species %in% c('Homo sapiens', 'Mus musculus'))) & apply.FS) {
        message("Cannot apply FS filter as the species is not 'Homo sapiens' or 'Mus Musculus'!")
        apply.FS <- FALSE
    }

    if (apply.CORR & (is.na(x0.path) | is.na(x1.path))) {
        message("Cannot apply CORR filter to data that is not given in matrix form!")
        apply.CORR <- FALSE
    }

    if (!is.null(custom_serine_threonine_kinase_data.path) && !dir.exists(custom_serine_threonine_kinase_data.path)) {
        stop('Custom serine threonine kinase data path cannot be found!')
    } else if (!is.null(custom_serine_threonine_kinase_data.path)) {
        serine_threonine_kinase_data <- Kinference::load_kinase_data(custom_serine_threonine_kinase_data.path)
    } else {
        serine_threonine_kinase_data <- list(
            kinase_motifs = Kinference::serine_threonine_kinase_data_chunk1$kinase_motifs,
            kinase_aprior_distributions = c(Kinference::serine_threonine_kinase_data_chunk1$combined_kinase_aprior_distributions, Kinference::serine_threonine_kinase_data_chunk2$combined_kinase_aprior_distributions, Kinference::serine_threonine_kinase_data_chunk3$combined_kinase_aprior_distributions),
            kinase_name_mappings = Kinference::serine_threonine_kinase_data_chunk1$kinase_name_mappings
        )
    }
    if (!is.null(custom_tyrosine_kinase_data.path) && !dir.exists(custom_tyrosine_kinase_data.path)) {
        stop('Custom tyrosine kinase data path cannot be found!')
    } else if (!is.null(custom_tyrosine_kinase_data.path)) {
        tyrosine_kinase_data <- Kinference::load_kinase_data(custom_tyrosine_kinase_data.path)
    } else {
        tyrosine_kinase_data <- Kinference::tyrosine_kinase_data
    }

    message("Loading data!")
    # Prepare data
    if (!is.na(x0.path) & !is.na(x1.path)) {
        data <- Kinference::prepare_intensity_matrices(
            x0.path = x0.path, 
            x1.path = x1.path,
            output.path = output.path,
            output.id = output.id,
            species = species,
            translate_uniprots = translate_uniprots,
            paired.samples = paired.samples,
            apply.log2 = apply.log2
        )

        if (ncol(data[['x0']]) < m & apply.CORR) {
            message("Won't apply CORR filter because the number of measurements threshold for the correlation coefficient computation 'm' is bigger than the number of columns!")
            apply.CORR <- FALSE
        }
    } else if (!is.na(f.path)) {
        data <- Kinference::prepare_intensity_vector(
            f.path = f.path,
            output.path = output.path,
            output.id = output.id,
            species = species,
            translate_uniprots = translate_uniprots
        )
    } else {
        stop("'x0.path' and 'x1.path' or 'f.path' have to be given!")
    }

    for (n_instance in n) {
        message("Starting KINference analysis with n = ", n_instance)
        for (alpha_instance in alpha) {        
            message("Starting KINference analysis with alpha = ", alpha_instance)
            # Infer serine/threonine and tyrosine baseline KIN
            message("Computing baseline KIN!")
            baseline_KINs <- Kinference::infer_baseline_KIN(
                md = data[['md']],
                output.path = output.path,
                output.id = paste0(output.id, '_n', n_instance, '_alpha', alpha_instance),
                n = n_instance,
                alpha = alpha_instance,
                serine_threonine_kinase_data = serine_threonine_kinase_data,
                tyrosine_kinase_data = tyrosine_kinase_data,
                translate_uniprots = translate_uniprots
            )

            for (gamma_instance in gamma) {
                message("Starting KINference analysis with gamma = ", gamma_instance)
                # Compute kinase enrichments
                message("Computing kinase enrichments!")
                Kinference::compute_kinase_enrichment(
                    md = data[['md']],
                    baseline_KINs = baseline_KINs,
                    output.path = output.path,
                    output.id = paste0(output.id, '_n', n_instance, '_alpha', alpha_instance, '_gamma', gamma_instance),
                    gamma = gamma_instance,
                    serine_threonine_kinase_data = serine_threonine_kinase_data,
                    tyrosine_kinase_data = tyrosine_kinase_data
                )

                for (beta_instance in beta) {
                    message("Starting KINference analysis with beta = ", beta_instance)
                    # Compute node filters
                    if (apply.DIFF | apply.FS) {
                        message("Computing node filters!")
                        KIN <- Kinference::compute_node_filters(
                            md = data[['md']],
                            KIN = baseline_KINs[['combined']],
                            output.path = output.path,
                            output.id = paste0(output.id, '_n', n_instance, '_alpha', alpha_instance, '_gamma', gamma_instance, '_beta', beta_instance),
                            gamma = gamma_instance,
                            beta = beta_instance,
                            apply.DIFF = apply.DIFF,
                            apply.FS = apply.FS,
                            serine_threonine_kinase_data = serine_threonine_kinase_data,
                            tyrosine_kinase_data = tyrosine_kinase_data,
                            translate_uniprots = translate_uniprots
                        )
                    } else {
                        KIN <- baseline_KINs[['combined']]
                    }

                    # Compute edge filters
                    if (apply.CORR | apply.PCST) {
                        message("Computing edge filters!")


                        if (apply.CORR) {
                            apply.PCST = FALSE
                            for (delta_instance in delta) {
                                message("Starting KINference analysis with delta = ", delta_instance)
                                edge_filtered_KINs <- Kinference::compute_edge_filters(
                                    data = data,
                                    KIN = KIN,
                                    output.path = output.path,
                                    output.id = paste0(output.id, '_n', n_instance, '_alpha', alpha_instance, '_gamma', gamma_instance, '_beta', beta_instance, '_delta', delta_instance),
                                    apply.CORR = apply.CORR,
                                    apply.PCST = apply.PCST,
                                    gamma = gamma_instance,
                                    delta = delta_instance,
                                    epsilon = epsilon,
                                    m = m,
                                    multiple_testing_correction = multiple_testing_correction,
                                    serine_threonine_kinase_data = serine_threonine_kinase_data,
                                    tyrosine_kinase_data = tyrosine_kinase_data,
                                    translate_uniprots = translate_uniprots
                                )
                            }
                            apply.PCST = TRUE
                            apply.CORR = FALSE
                        }

                        if (apply.PCST) {
                            reticulate::virtualenv_create(
                                envname = paste0(.libPaths(), "/r-Kinference"), 
                                packages = c("numpy==1.26.4", "pandas", "pcst_fast")
                            )
                            reticulate::use_virtualenv(paste0(.libPaths(), "/r-Kinference"))
                        }
                        edge_filtered_KINs <- Kinference::compute_edge_filters(
                            data = data,
                            KIN = KIN,
                            output.path = output.path,
                            output.id = paste0(output.id, '_n', n_instance, '_alpha', alpha_instance, '_gamma', gamma_instance, '_beta', beta_instance),
                            apply.CORR = apply.CORR,
                            apply.PCST = apply.PCST,
                            gamma = gamma_instance,
                            delta = delta[1],
                            epsilon = epsilon,
                            m = m,
                            multiple_testing_correction = multiple_testing_correction,
                            serine_threonine_kinase_data = serine_threonine_kinase_data,
                            tyrosine_kinase_data = tyrosine_kinase_data,
                            translate_uniprots = translate_uniprots
                        )
                        
                    }  
                }
            }
        }
    }
    message("Done!")
}
