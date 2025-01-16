#' Run KINference Analysis
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
run_KINference <- function(
    x0.path = NA,
    x1.path = NA,
    f.path = NA,
    output.path = 'results',
    output.id = 'key',
    species = 'Homo sapiens',
    translate_uniprots = TRUE,
    paired.samples = TRUE,
    apply.log2 = FALSE,
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
    message(paste0(
        '============= Parameters: ============= \n',
        'x0.path: ', x0.path, '\n',
        'x1.path: ', x1.path, '\n',
        'f.path: ', f.path, '\n',
        'output.path: ', output.path, '\n',
        'output.id: ', output.id, '\n',
        'paired.samples: ', paired.samples, '\n',
        'apply.log2: ', apply.log2, '\n',
        'species: ', species, '\n',
        'translate_uniprots: ', translate_uniprots, '\n',
        'apply.DIFF: ', apply.DIFF, '\n',
        'apply.FS: ', apply.FS, '\n',
        'apply.CORR: ', apply.CORR, '\n',
        'apply.PCST: ', apply.PCST, '\n',
        'beta: ', beta, '\n',
        'gamma: ', gamma, '\n',
        'delta: ', delta, '\n',
        'epsilon: ', epsilon, '\n',
        'm: ', m, '\n',
        'multiple_testing_correction: ', multiple_testing_correction, '\n',
        'custom_serine_threonine_kinase_data.path: ', custom_serine_threonine_kinase_data.path, '\n',
        'custom_tyrosine_kinase_data.path: ', custom_tyrosine_kinase_data.path, '\n',
        '========================================'
    ))

    # Checking validity of input arguments
    if (n <= 0 | alpha <= 0 | epsilon <= 0 | m <= 0) {
        stop("'n', 'alpha', 'epsilon', and 'm' must all be greater than 0!")
    }

    # Deactivating filters if the corresponding parameter is 0
    if (gamma == 0) {
        message('Wont apply Diff filter because gamma is 0!')
        apply.DIFF = FALSE
    }
    if (beta == 0) {
        message('Wont apply FS filter because beta is 0!')
        apply.FS = FALSE
    }
    if (delta == 0) {
        message('Wont apply CORR filter because delta is 0!')
        apply.CORR = FALSE
    }

    if ((!(species %in% c('Homo sapiens', 'Mus musculus'))) & apply.FS) {
        message("Cannot apply FS filter because the species is not 'Homo sapiens' or 'Mus musculus'!")
        apply.FS <- FALSE
    }

    if (apply.CORR & (is.na(x0.path) | is.na(x1.path))) {
        message("Cannot apply CORR filter to data that is not given in matrix form!")
        apply.CORR <- FALSE
    }

    if (!(species %in% c('Homo sapiens', 'Mus musculus'))) {
        message("Selected species does not exist! This may cause the script to break. Please select either 'Homo sapiens' or 'Mus musculus'")
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

    if (species == 'Homo sapiens' && !translate_uniprots) {
        message("Setting 'translate_uniprot' to TRUE because the species is set to 'Homo sapiens'.")
        translate_uniprots <- TRUE
    } 

    if (apply.PCST && !translate_uniprots && species != 'Homo sapiens') {
        message("Disabling PCST filter because the Uniprot IDs are not translated to human Uniprot IDs!")
        apply.PCST <- FALSE
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
        if ((ncol(data[['x0']]) < m) & apply.CORR) {
            message("Won't apply CORR filter because the number of measurements threshold for the correlation coefficient computation 'm' is bigger than the number of columns!")
            apply.CORR <- FALSE
        }
        if ((ncol(data[['x0']]) != ncol(data[['x1']])) & apply.CORR) {
            message("Won't apply CORR filter because the number of columns in the two matrices are not equal!")
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

    message("Computing baseline KIN!")
    # Infer serine/threonine and tyrosine baseline KIN
    baseline_KINs <- Kinference::infer_baseline_KIN(
        md = data[['md']],
        output.path = output.path,
        output.id = output.id,
        n = n,
        alpha = alpha,
        serine_threonine_kinase_data = serine_threonine_kinase_data,
        tyrosine_kinase_data = tyrosine_kinase_data,
        translate_uniprots = translate_uniprots
    )

    if (gamma > 0) {
        message("Computing kinase enrichments!")
        # Compute kinase enrichments
        Kinference::compute_kinase_enrichment(
            md = data[['md']],
            baseline_KINs = baseline_KINs,
            output.path = output.path,
            output.id = output.id,
            gamma = gamma,
            serine_threonine_kinase_data = serine_threonine_kinase_data,
            tyrosine_kinase_data = tyrosine_kinase_data
        )
    } else {
        message("Will not compute kinase enrichments because gamma is 0!")
    }

    # Compute node filters
    if (apply.DIFF | apply.FS) {
        message("Computing node filters!")
        KIN <- Kinference::compute_node_filters(
            md = data[['md']],
            KIN = baseline_KINs[['combined']],
            output.path = output.path,
            output.id = output.id,
            gamma = gamma,
            beta = beta,
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
            output.id = output.id,
            apply.CORR = apply.CORR,
            apply.PCST = apply.PCST,
            gamma = gamma,
            delta = delta,
            epsilon = epsilon,
            m = m,
            multiple_testing_correction = multiple_testing_correction,
            serine_threonine_kinase_data = serine_threonine_kinase_data,
            tyrosine_kinase_data = tyrosine_kinase_data,
            translate_uniprots = translate_uniprots
        )
    }
    message("Done!")
}
