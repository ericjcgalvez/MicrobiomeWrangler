beta_ordination2 <- function(phylo, distance_method, myf, ordination_type,
                             color, shape) {

  dist <- phyloseq::distance(phylo, distance_method)
  phylo_sample_data <- as(sample_data(phylo), "data.frame")

  # Adonis Test
  # adonis_formula <- paste("vegan::adonis(d ~", myf,
  #                     ",phylo_sample_data, permutations = 999, parallel = 2 )")
  # adonis_res <- eval(parse(text = adonis_formula))
  myadonis <- paste("vegan::adonis(dist ~", myf, ",phylo_sample_data, permutations =999, parallel = 2 )")
  adonis <- eval(parse(text = myadonis))
  # myf <- as.formula(paste("dist ~ ", deparse(substitute(myf))))
  # adonis_res <- vegan::adonis(eval( myf),
  #                             data = phylo_sample_data,
  #                             permutations = permutations,
  #                             parallel = parallel)

  # Ordination
  phylo_ord <- phyloseq::ordinate(pyhlo = phylo, ordination_type, distance_method)

  ## Ordination Plot
  phylo_ord_plot <- phyloseq::plot_ordination(pyhlo, phylo_ord,
                                              type = "samples",
                                              color = color,
                                              shape = shape)
  phylo_ord_plot <- phylo_ord_plot +
    ggplot2::geom_point(size = 5, alpha = 0.70) +
    ggthemes::scale_colour_gdocs() +
    ggpubr::stat_stars() +
    ggpubr::stat_mean(geom = "point", size = 4) +
    ggplot2::scale_shape_manual(values = c(16, 15, 17, 18:30))
  ggplot2::theme(legend.position = "bottom")

  phylo_ord_plot
}

facet_barplot <- function(phylo,
                          taxlev = "Family",
                          orderbar = "",
                          facetby = "Description",
                          ntaxa = 14) {
  ###
  # First part need to be put in another function
  ###

  # fetch user input to be able to unquote them later
  taxlev <- rlang::enexpr(taxlev)
  facetby <- rlang::ensym(facetby)
  if(orderbar != "") orderbar_sym <- rlang::ensym(orderbar)

  # fetch faulty user input
  if (ntaxa <= 0) stop("ntaxa must be higher than 0")
  if (as.integer(ntaxa) != ntaxa) stop("ntaxa must be an integer")

  # Subset OTUs with highest abundance for plotting
  glomed_psq <- phyloseq::tax_glom(phylo, taxlev)
  glomed_top14_psq <- phyloseq::prune_taxa(names(sort(phyloseq::taxa_sums(glomed_psq), T))[1:ntaxa], glomed_psq)

  ## Adding taxon Others in the phyloseq object
  # calculating "others" values
  others <- abs(sample_sums(glomed_psq) - sample_sums(glomed_top14_psq))

  # create "Others" taxa
  others.otu <- rbind(others, glomed_top14_psq@otu_table)
  others.otu <- otu_table(others.otu, taxa_are_rows = TRUE)

  # adding "Others" taxa and taxonomy to the provided phyloseq object
  glomed_top14_psq@otu_table <- others.otu
  others <- rep("Others", ncol(tax_table(glomed_top14_psq)))
  glomed_top14_psq@tax_table <- tax_table(rbind(others, tax_table(glomed_top14_psq)))

  # make phyloseq object tidy because ggplot2 requires it
  melted.dt <- psmelt(glomed_top14_psq)

  # Plotting -----------------------------------------
  # order ranks
  colour_count <- length(unique(eval(rlang::expr(`$`(melted.dt, !!taxlev)))))
  taxlev <- rlang::ensym(taxlev)
  melted.dt <- dplyr::arrange(melted.dt, !!taxlev)
  melted.dt[[taxlev]] <- forcats::fct_relevel(melted.dt[[taxlev]], "Others", after = Inf)

  p <- ggplot2::ggplot(melted.dt)
  if (orderbar == "") {
    p <- p + ggplot2::aes(x = SampleID)
  } else {
    p <- p + ggplot2::aes(x = reorder(SampleID, !!orderbar_sym))
  }
  p <- p + ggplot2::aes(
    y = Abundance,
    fill = !!taxlev)
  p <- p + ggplot2::geom_bar(stat = "identity")
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.9, size = 8)) +
    ggplot2::ylab("Abundance") +
    ggplot2::facet_wrap(rlang::expr(~ !!facetby),
                        scales = "free_x",
                        labeller = ggplot2::label_wrap_gen(multi_line = FALSE)
    ) +
    ggplot2::scale_fill_manual(values = ggpubr::get_palette(palette = "Paired", k = colour_count)) +
    ggplot2::xlab("Sample")
  # ggplot2::theme(axis.text.x = ggplot2::element_blank())
  return(p)
}
