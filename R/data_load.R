#' Read 16S amplicon biom file into phyloseq object
#'
#' This function allows you to load the data recived from amplicon sequencing
#' from MIKI-lab.
#' @param otu character string. Otu_table file path. In biom file oldformat
#' from qiime 1.8.
#' @param map character string. Mapping file path. In biom file oldformat from
#' qiime 1.8.
#' @param tree character string. tree file path. In nexus format from qiime 1.8.
#' @param min_otu numeric. Filter OTUs that do not reach this count number
#' @param min_sample numeric. Filter samples that do not reach this count number
#' @param verbose boolean. If TRUE cat()s step progress
#' @keywords load amplicon
#' @export
#'
#' @importFrom phyloseq import_biom phy_tree taxa_names sample_data tax_table
#' @importFrom phyloseq sample_sums taxa_sums
#'
#' @examples
#' phylo <- MIKIbiomeR::data_load(otu = otu_path, map = map_path,
#' tree = tree_path)

data_load <- function(otu, map, tree, min_otu = 0.2, min_sample = 999,
                      verbose = FALSE) {
  # import modelue, to be update with JSON biom
  if (is.null(otu)) {
    stop("otutable file not supplied\n")
  }
  if (is.null(map)) {
    stop("maping_file not supplied\n")
  }

  if (verbose) cat("Loading data and making a phyloseq object.\n")

  # load data
  phylo_import <- import_biom(otu,
                              treefilename = NULL,
                              refseqfilename = NULL,
                              refseqFunction = readDNAStringSet,
                              refseqArgs = NULL,
                              parallel = FALSE, version = 1.0)

  otu_tbl <- phyloseq::otu_table(phylo_import)
  tax_tbl <- phyloseq::tax_table(phylo_import)
  #meta_tbl <- phyloseq::import_qiime_sample_data(map). # not working in new R version
  meta_tbl <- phyloseq::sample_data(read.delim(file = map))
  rownames(meta_tbl) <- meta_tbl$X.SampleID

  tree_nexus <- phyloseq::phy_tree(phyloseq::import_qiime(treefilename = tree))

  # create phyloseq object
  phylo <- phyloseq::phyloseq(otu_tbl, tax_tbl, meta_tbl, tree_nexus)

  phyloseq::phy_tree(phylo) <- ape::root(phy_tree(phylo),
                                         sample(taxa_names(phylo), 1),
                                         resolve.root = TRUE)
  ape::is.rooted(phy_tree(phylo))
  colnames(phyloseq::tax_table(phylo)) <- c("Kingdom",
                                  "Phylum",
                                  "Class",
                                  "Order",
                                  "Family",
                                  "Genus",
                                  "Species")

  # !!! this needs to be replaced by something better
  # sometimes X. gets introduced
  if (any(grepl("^X\\.", colnames(phyloseq::sample_data(phylo))))) {
    warning("'X.' found in sample data. Will be removed\n")
    colnames(phyloseq::sample_data(phylo)) <- sub("^X\\.", "",
                                        colnames(sample_data(phylo)))
  }
  phyloseq::tax_table(phylo) <- gsub("^(.+)__", "", tax_table(phylo))

  if (verbose) cat("Loading successful\n")

  ## consider to add the Decotamination step wit "decontam"

  # Filtering module
  if (verbose) cat("Performing OTU rel. abundance filtering OTUs >",
                  min_otu,
                  "\nPerforming sample filtering  >",
                  min_sample,
                  "\n")

  ps <- phyloseq::prune_taxa(taxa_sums(phylo) > 0, phylo)
  otu <- phyloseq::otu_table(ps)
  otu_rel <- prop.table(otu, margin = 2) * 100
  keep <- apply(otu_rel, 1, max) >= min_otu
  phylo_clean <- phyloseq::prune_taxa(keep, ps)

  # Taxa cleaning and prunning
  phyloseq::tax_table(phylo_clean)[is.na(
    tax_table(phylo_clean))] <- "Ambiguous_taxa"
  phyloseq::tax_table(phylo_clean) <- sub("uncultured bacterium",
                                          "Ambiguous_taxa",
                                          tax_table(phylo_clean))
  phyloseq::tax_table(phylo_clean) <- sub("uncultured", "Ambiguous_taxa",
                                          tax_table(phylo_clean))
  clean_ps <- phyloseq::subset_taxa(phylo_clean, Kingdom == "Bacteria" &
    Family != "Mitochondria" &
    Class != "Chloroplast" &
    Order != "Ambiguous_taxa")

  phylo_clean <- phyloseq::prune_samples(sample_sums(phylo_clean) >= min_sample,
                                         phylo_clean)
  phylo_clean <- phyloseq::prune_taxa(taxa_sums(phylo_clean) > 0, phylo_clean)

  phylo_clean
}
