#' Alpha diversity
#'
#' Calulates alpha diversity distances
#'
#' @param phylo reads a phyloseq object
#' @param measures character vector. Type of alpha diversity indexes. Possible
#' values: c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson",
#' "BF_ratio"). "BF_ratio" calculates the ratio of Bacteroidetes and Firmicutes.
#' @param  min_count samples that do not reach this count number
#'
#' @keywords Alpha diversity..
#' @export
#' @examples
#' data("mice_B6_N")
#' alpha_diversity_plus_dt(mice_B6_N, 2000)
#' @importFrom magrittr %>%
#' @importFrom phyloseq prune_samples
#' @importFrom phyloseq sample_sums
#' @importFrom phyloseq rarefy_even_depth
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq estimate_richness
#' @importFrom phyloseq tax_glom
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq subset_taxa

alpha_div_calc <- function(phylo,
                           min_count = 2000,
                           measures = c("Observed", "Chao1", "ACE", "Shannon",
                                        "Simpson", "InvSimpson", "BF_ratio")) {

  if (class(phylo)[1] != "phyloseq") stop("Input is not a phyloseq object!")

  ## Calculate alpha diversity
  # this steps is because phyloseq further complains for the rarefaction step
  phylo_pruned <- prune_samples(sample_sums(phylo) >= min_count, phylo)
  phylo_pruned_rarefied <- rarefy_even_depth(phylo_pruned,
                                             sample.size = min_count)

  # TO DO: For some reason SampleIDs are not seen as valid rownames,
  # thats why an additional X gets concatonated at the beginning of an ID
  # BEFORE: SampleID: 123
  # AFTER: SampleID: X123
  diversityindex <- measures %>%
    estimate_richness(phylo_pruned_rarefied, measures = .) %>%
    as.matrix(., row.names = 1)

  if (length(grep("^X", row.names(diversityindex))) > 0) {
    warning("X in the first position of Sample Identifier is getting deleted\n")
    row.names(diversityindex) <- sub("^X", "", row.names(diversityindex))
  }


  # Add alpha diversity to metadata

  phyloseq::sample_data(phylo) <- cbind(
    sample_data(phylo_pruned_rarefied),
    diversityindex[rownames(sample_data(phylo_pruned_rarefied)), ]
  )

  if ("BF_ratio" %in% measures) phylo <- bacfir_ratio_calc(phylo)

  return(phylo)
}

#' Alpha diversity plot
#'
#' Plots boxplot of alpha-diversity distance together with jitter plot.
#'
#' @param phylo phyloseq object with alpha-diversity distances in the
#' sample data
#' @param x_axis Determines x variable in plots. Column name in sample data
#' from phyloseq object as character
#' @param color_var Variable in sample data from phyloseq object
#' for color aesthetics
#' @param ... More Paremeters
#' @param jitter_params Additional parameters for ggplot2 geom_jitter
#' @param boxplot_params Additional parameters for ggplot2 geom_boxplot
#' @param measure Alpha-diversity distance. Character string of column name of
#' sample data from the phyloseq object. Normally calculated with
#' alpha_div_calc(). Examples: "Observed"; "Simpson";  "Shannon"; "BF_ratio"
#'
#' @keywords alpha-diversity plot
#' @export
#'
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#' @examples
#' data("mice_B6_N")
#' phylo <- phyloseq::subset_samples(mice_B6_N, Description == "IN")
#' phylo <- phyloseq::prune_taxa(phyloseq::taxa_sums(mice_B6_N) > 0, phylo)
#' phylo_a_dist <- alpha_div_calc(phylo, 2000)
#' alpha_div_plot(phylo_a_dist, "Vendor", "Vendor", "Simpson")
alpha_div_plot <- function(phylo,
                            x_axis,
                            color_var,
                            alpha_dist,
                            jitter_params = list(),
                            boxplot_params = list(),
                            ...) {
  alpha_tbl <- data.frame(sample_data(phylo))

  # sanity check
  if (is.null(alpha_tbl[[alpha_dist]])) {
    stop(paste(alpha_dist, "is not a variable of", deparse(substitute(phylo))))
  }

  # handle custome plot parameter options
  params <- list(...)
  jitter_params <- modifyList(params, jitter_params)
  boxplot_params <- modifyList(params, boxplot_params)

  jitter <- do.call("geom_point",
                    modifyList(
                      list(position = ggplot2::position_jitterdodge(dodge.width=0.9),
                           alpha = 0.5),
                      jitter_params
                    ))
  boxplot <- do.call("geom_boxplot",
                     modifyList(
                       list(position = ggplot2::position_dodge(width=0.9),
                            fill = "white"),
                       boxplot_params
                     ))

  alpha_plot <- ggplot(alpha_tbl, aes(
    x = .data[[x_axis]],
    y = .data[[alpha_dist]],
    color = .data[[color_var]],
    fill = .data[[color_var]]
  )) +
    boxplot +
    jitter +
    ggplot2::xlab(x_axis) +
    ggplot2::ylab(alpha_dist) +
    ggplot2::labs(color = color_var) +
    ggplot2::guides(fill = FALSE)

  alpha_plot
}


#' Alpha diversity plotting
#'
#' Plots boxplots of alpha-diversity distances with jitter plot.
#'
#' @param phylo phyloseq object with alpha-diversity distances in the
#' sample data
#' @param x_axis Determines x variable in plots. Column name in sample data
#' from phyloseq object as character
#' @param color_var Variable in sample data from phyloseq object
#' for color aesthetics
#' @param ... More Paremeters
#' @param jitter_params Additional parameters for ggplot2 geom_jitter
#' @param boxplot_params Additional parameters for ggplot2 geom_boxplot
#' @param measures Alpha-diversity distances. Atomic character vector with
#' column names of sample data from the phyloseq object. Normally calculated
#' with alpha_div_calc()
#'
#' @keywords alpha-diversity plot
#' @export
#' @examples
#' data("mice_B6_N")
#' phylo <- phyloseq::subset_samples(mice_B6_N, Description == "IN")
#' phylo <- phyloseq::prune_taxa(phyloseq::taxa_sums(mice_B6_N) > 0, phylo)
#' phylo_a_dist <- alpha_div_calc(phylo, 2000)
#' alpha_div_plots(phylo_a_dist, "Vendor", "Vendor")
alpha_div_plots <- function(phylo,
                            x_axis,
                            color_var,
                            measures = c("Observed",
                                         "Simpson",
                                         "Shannon",
                                         "BF_ratio"),
                            jitter_params = list(),
                            boxplot_params = list(),
                            ...) {
  ## for to do: Add the contrast groups for stat_compare_means




# Plots
  alpha_plot_list <- measures %>% purrr::map(
    function(measure)
      alpha_div_plot(phylo, x_axis, color_var, measure,
                     jitter_params, boxplot_params, ...)
  )
  alpha_plot_list
}



#' Calculate bacteroides to firmicutes ratio
#'
#' @param phylo phyloseq object
#'
#'
#' @return phyloseq object with additional variable in sample_data
#'
#'
#' @examples
#' data("mice_B6_N")
#' bacfir_ratio(mice_B6_N)

bacfir_ratio_calc <- function(physeq = phylo) {

  # Collapse on phyla
  phyla_psq <- tax_glom(physeq = physeq, taxrank = "Phylum")

  # Keep B/F taxa
  tax_table(phyla_psq)
  phyla_rel_bact <- otu_table(subset_taxa(phyla_psq, Phylum == "Bacteroidetes"))
  phyla_rel_firm <- suppressWarnings(
    otu_table(subset_taxa(phyla_psq, Phylum == "Firmicutes")))

  # Calculate ratio
  bf_ratio <- phyla_rel_bact / (phyla_rel_bact + phyla_rel_firm)

  # Add to sample metadata
  phyloseq::sample_data(physeq)$BF_ratio <- as.numeric(bf_ratio)

  # Return phyloseq object
  return(physeq)
}
