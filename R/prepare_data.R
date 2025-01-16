#' Prepare Intensity Matrices
#'
#' This function reads in two data files containing protein intensity measurements, processes the data, and returns a list containing metadata and the processed intensity matrices.
#'
#' @param x0.path A string specifying the file path to the first data file.
#' @param x1.path A string specifying the file path to the second data file.
#' @param output.path A string specifying the directory where the output should be saved. Default is 'results'.
#' @param output.id A string specifying an identifier for the output. Default is 'key'.
#' @param species A string specifying the species of the data. Default is 'Homo sapiens'. If 'Mus musculus', Uniprot IDs will be converted to human Uniprot IDs.
#' @param translate_uniprots A boolean value indicating if the Uniprot IDs of other species should be translated to human Uniprot IDs. Default is TRUE.
#' @param paired.samples A logical value indicating whether the samples are paired. Default is FALSE.
#' @param apply.log2 A logical value indicating whether to apply log2 transformation to the data. Default is FALSE.
#'
#' @return A list containing:
#' \item{md}{A tibble with processed metadata information for the proteins.}
#' \item{x0}{A matrix with processed intensity values from the first data file.}
#' \item{x1}{A matrix with processed intensity values from the second data file.}
#' @export
prepare_intensity_matrices <- function(
    x0.path,
    x1.path,
    output.path = 'results',
    output.id = 'key',
    species = 'Homo sapiens',
    translate_uniprots = TRUE,
    paired.samples = FALSE,
    apply.log2 = FALSE
) {

    # Reading in data
    x0 <- readr::read_tsv(x0.path, show_col_types = FALSE)
    if (!('Protein' %in% colnames(x0))) {
        stop("x0 does not have a column named 'Protein'!")
    }
    x0 <- tibble::column_to_rownames(x0, 'Protein')
    x0 <- x0[, order(colnames(x0))]

    x1 <- readr::read_tsv(x1.path, show_col_types = FALSE)
    if (!('Protein' %in% colnames(x1))) {
        stop("x1 does not have a column named 'Protein'!")
    }
    x1 <- tibble::column_to_rownames(x1, 'Protein')
    x1 <- x1[, order(colnames(x1))]

    # Remove all proteins that are only measured in one condition
    x1 <- x1[which(rownames(x1) %in% rownames(x0)), ]
    x0 <- x0[which(rownames(x0) %in% rownames(x1)), ]

    # Order row names
    x0 <- x0[order(rownames(x0)), ]
    x1 <- x1[order(rownames(x1)), ]

    if (paired.samples) {
        # Adding NAs of one matrix to the other one to enforce paired samples
        x0[which(is.na(x1), arr.ind = T)] <- NA
        x1[which(is.na(x0), arr.ind = T)] <- NA
    }

    # Filter out all Proteins with no measurements (only NAs in the row)
    x0 <- x0[rowSums(is.na(x0)) != ncol(x0), ]
    x1 <- x1[rowSums(is.na(x1)) != ncol(x1), ]

    # Apply log2 transformation
    if (apply.log2) {
        x0 <- log2(x0)
        x1 <- log2(x1)
    }

    # Creating data object
    md <- dplyr::tibble(
      Protein = rownames(x1),
      Uniprot = sapply(strsplit(Protein, '_'), function(x) { x[1] }),
      AAPos = sapply(strsplit(Protein, '_'), function(x) { x[2] }),
      AA = sapply(strsplit(AAPos, '_'), function(x) { stringr::str_sub(x, 1, 1) }),
      Pos = sapply(strsplit(AAPos, '_'), function(x) { as.numeric(stringr::str_sub(x, 2)) }),
      Uniprot_Pos = paste0(Uniprot, '_', Pos)
    )

    # Convert mouse uniprot IDs to human uniprot IDs 
    if (species == 'Mus musculus' && translate_uniprots) {
        md$Protein_musMusculus <- md$Protein
        md <- OmnipathR::orthology_translate_column(
            data = md, 
            column = 'Uniprot', 
            target_organism = 'human', 
            source_organism = 'mouse',
            replace = TRUE
        )
        md$Protein <- paste0(md$Uniprot, '_', md$AA, md$Pos)
        md$Uniprot_Pos <- paste0(md$Uniprot, '_', md$Pos)
        
        x0 <- x0[rownames(x0) %in% md$Protein_musMusculus, ]
        x1 <- x1[rownames(x1) %in% md$Protein_musMusculus, ]
        # Changing rownames
        rownames(x0) <- sapply(rownames(x0), function(x) { md[md$Uniprot_musMusculus == x, 'Protein_musMusclus'] })
        rownames(x1) <- sapply(rownames(x1), function(x) { md[md$Uniprot_musMusculus == x, 'Protein_musMusclus'] })
    }

    # remove all proteins that are not phosphorylated at serine (S), threonine (T) or tyrosine (Y)
    ind_to_keep <- which(md$AA %in% c('S', 'T', 'Y'))
    x1 <- x1[ind_to_keep, ]
    x0 <- x0[ind_to_keep, ]
    md <- md[ind_to_keep, ]

    # Compute log2FC
    if (!paired.samples) {
        md$f <- rowMeans(x1, na.rm = T) - rowMeans(x0, na.rm = T)
    } else {
        md$f <- rowMeans(x1 - x0, na.rm = T)
    }

    return(list(md=md, x0=x0, x1=x1))
}

#' Prepare Intensity Vector
#'
#' This function reads in one data file containing the Protein names, phosphorylated amino acids, phosphorylation positions and the log2FCs.
#'
#' @param f.path A string representing the file path to the input TSV file.
#' @param output.path A string representing the output directory path. Default is 'results'.
#' @param output.id A string representing the output identifier. Default is 'key'.
#' @param species A string representing the species. Default is 'Homo sapiens'. If 'Mus Musculus', mouse Uniprot IDs will be converted to human Uniprot IDs.
#' @param translate_uniprots A boolean value indicating if the Uniprot IDs of other species should be translated to human Uniprot IDs. Default is TRUE.
#' @return Tibble with preprocessed data
#' @export
prepare_intensity_vector <- function(
    f.path,
    output.path = 'results',
    output.id = 'key',
    species = 'Homo sapiens',
    translate_uniprots = TRUE
) {

    # Reading in data
    f <- readr::read_tsv(f.path, show_col_types = FALSE)
    if (!('Protein' %in% colnames(f))) {
        stop("f does not have a column named 'Protein'!")
    }
    f <- tibble::column_to_rownames(f, 'Protein')

    # Creating meta data
    md <- dplyr::tibble(
        Protein = rownames(f),
        Uniprot = sapply(strsplit(Protein, '_'), function(x) { x[1] }),
        AAPos = sapply(strsplit(Protein, '_'), function(x) { x[2] }),
        AA = sapply(strsplit(AAPos, '_'), function(x) { stringr::str_sub(x, 1, 1) }),
        Pos = sapply(strsplit(AAPos, '_'), function(x) { as.numeric(stringr::str_sub(x, 2)) }),
        Uniprot_Pos = paste0(Uniprot, '_', Pos),
        f = f[[1]]
    )

    # Convert mouse uniprot IDs to human uniprot IDs 
    if (species == 'Mus musculus' && translate_uniprots) {
        md$Protein_musMusculus <- md$Protein
        md <- OmnipathR::orthology_translate_column(
            data = md, 
            column = 'Uniprot', 
            target_organism = 'human', 
            source_organism = 'mouse',
            replace = TRUE
        )
        md$Protein <- paste0(md$Uniprot, '_', md$AA, md$Pos)
        md$Uniprot_Pos <- paste0(md$Uniprot, '_', md$Pos)
    }

    # remove all proteins that are not phosphorylated at serine (S), threonine (T) or tyrosine (Y)
    md <- md[which(md$AA %in% c('S', 'T', 'Y')), ]

    return (list(md=md))
}

#' Load Kinase Data
#'
#' This function loads kinase data from a specified file path.
#'
#' @param kinase_data.path A string representing the file path to the kinase data.
#'
#' @return The loaded kinase data.
#'
#' @export
load_kinase_data <- function(
    kinase_data.path
) {
    kinase_motifs.path <- paste0(kinase_data.path, '/kinase_motifs.csv')
    kinase_name_mappings.path <- paste0(kinase_data.path, '/kinase_name_mappings.tsv')
    kinase_aprior_distributions.path <- paste0(kinase_data.path, '/apriori_distributions/')

    # reading in kinase motif data
    kinase_motifs <- readr::read_csv(kinase_motifs.path, show_col_types = FALSE)
    # reading in mapping of kinases to uniprot id file
    kinase_names_to_uniprotmapping <- readr::read_tsv(kinase_name_mappings.path, show_col_types = FALSE)

    kinases <- kinase_names_to_uniprotmapping[['ACC#']]

    # Computing a-priori distributions of kinases based on log2 score database
    kinase_aprior_distributions <- list()
    for (i in seq_len(nrow(kinase_motifs))) {
        row <- kinase_motifs[i, ]
        kin <- row$KINASE
        kin_aprior_file_name <- paste(kinase_aprior_distributions.path, kinase_motifs$KINASE[i], '.txt', sep = '')
        aprior_distribution <- stats::ecdf(readr::read_tsv(kin_aprior_file_name, show_col_types = FALSE)[[1]])
    kinase_aprior_distributions <- append(kinase_aprior_distributions, aprior_distribution)
    }
    names(kinase_aprior_distributions) <- kinases

    return(list(kinase_motifs = as.data.frame(kinase_motifs), kinase_aprior_distributions = kinase_aprior_distributions, kinase_name_mappings = as.data.frame(kinase_names_to_uniprotmapping)))
}