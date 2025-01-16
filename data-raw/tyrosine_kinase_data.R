load_kinase_data <- function(
    kinase_motifs.path,
    kinase_aprior_distributions.path,
    kinase_name_mappings.path
) {
  
    # reading in kinase motif data
    kinase_motifs <- data.table::fread(kinase_motifs.path)
    # reading in mapping of kinases to uniprot id file
    kinase_names_to_uniprotmapping <- data.table::fread(kinase_name_mappings.path)

    kinases <- kinase_names_to_uniprotmapping[['ACC#']]

    # Computing a-priori distributions of kinases based on log2 score database
    kinase_aprior_distributions <- list()
    for (i in seq_len(nrow(kinase_motifs))) {
        row <- kinase_motifs[i, ]
        kin <- row$KINASE
        kin_aprior_file_name <- paste(kinase_aprior_distributions.path, kinase_motifs$KINASE[i], '.txt', sep = '')
        aprior_distribution <- stats::ecdf(data.table::fread(kin_aprior_file_name)[[1]])
        kinase_aprior_distributions <- append(kinase_aprior_distributions, aprior_distribution)
    }
    names(kinase_aprior_distributions) <- kinases

    return(list(kinase_motifs = kinase_motifs, kinase_aprior_distributions = kinase_aprior_distributions, kinase_name_mappings = kinase_names_to_uniprotmapping))
}

tyrosine_kinase_data <- load_kinase_data(
    kinase_motifs.path = 'kinase_data/tyrosine_kinases/kinase_motifs.csv',
    kinase_aprior_distributions.path = 'kinase_data/tyrosine_kinases/apriori_distributions/',
    kinase_name_mappings.path = 'kinase_data/tyrosine_kinases/kinase_name_mappings.tsv'
)

usethis::use_data(tyrosine_kinase_data, overwrite = TRUE)  
