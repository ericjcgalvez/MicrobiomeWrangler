#' Summarises Phyloseq object to taxlevel and highest n abundance
#'
#'
#' @param phylo phyloseq object
#' @param top_n single numeric, prunes taxa with highest n abundance
#' @param taxlev character. Taxonomic level.
#'
#' @return phyloseq object
#' @example
#' data("mice_B6_N")
#' tax_subset(mice_B6_N, 10, "Familiy")
tax_subset <- function(phylo, top_n, taxlev) {
  phylo_sum_tax <- phyloseq::tax_glom(phylo, taxlev)
  phylo_sum_abund <- sort(phyloseq::taxa_sums(phylo_sum_tax), T)
  phylo_reduced <- phyloseq::prune_taxa(names(phylo_sum_abund)[1:top_n], phylo_sum_tax)

  return(phylo_reduced)
}

#' Calculate relative abundance
#'
#'
#' @param phylo Phyloseq object
#'
#' @return phyloseq object
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' calc_rel_abund(mice_B6_N)
calc_rel_abund <- function(phylo, n = 1) {
  phylo_rel <- transform_sample_counts(phylo, function(x) n * x / sum(x))
  return(phylo_rel)
}


#' Export Phyloseq as tidy data frame
#'
#' @param phylo Phyloseq object
#'
#' It is slow, so only use when necessary
#'
#' @return tidy data frame
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' phyloseq_to_tidy(mice_B6_N)
phyloseq_to_tidy <- function(phylo) {
  sam_data <- data.frame(phylo@sam_data) %>%
    dplyr::mutate(SampleID = as.character(SampleID))

  otu_table_tidy <- data.frame(otu_table(phylo)) %>%
    tibble::rownames_to_column(var = "OTU") %>%
    tidyr::gather(SampleID, OTU_abundance,-OTU) %>%
    dplyr::mutate(SampleID = as.character(SampleID)) %>%
    dplyr::mutate(SampleID = substr(SampleID,2,length(SampleID)))

  tax_table_tidy <- data.frame(phylo@tax_table) %>%
    tibble::rownames_to_column(var = "OTU")

  phylo_tidy <- dplyr::left_join(sam_data, otu_table_tidy, by = "SampleID")
  phylo_tidy <- dplyr::left_join(phylo_tidy, tax_table_tidy, by = "OTU")
  return(phylo_tidy)
}
