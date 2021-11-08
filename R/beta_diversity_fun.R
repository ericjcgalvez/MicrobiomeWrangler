#'  Plot beta ordination function
#'
#' This function performs ordination plots and adonis test
#'
#' @param phylo Phyloseq object
#' @param formula stats for adonis i.e. c("Genotype + Experiment + Cage.ID")
#' @param color metadata category color
#' @param shape metadata category shape
#' @param distance_method Similarity or distance measurement i.e. "bray",
#' "jaccard", "uunifrac", "wunifrac"
#' @param ordination_type Type of ordination that should be performed
#' @param geom_star boolean. If TRUE create a star plot by drawing segments from
#' group centroid to each points.
#' @param permutations a list of control values for the permutations as
#' returned by the function how, or the number of permutations required,
#' or a permutation matrix where each row gives the permuted indices.
#' @param parallel Parallel computations
#' @param tagtype Label type of the plots. i.e. "A" "1" "I".
#' @param verbose If TRUE calculation results will be printed. If FALSE not.
#' @param layout List of patch_area objects. Control the plot layout.
#'
#' @importFrom magrittr %>%
#' @importFrom patchwork area
#'
#' @keywords Beta diversity
#' @export
#' @examples
#' data("mice_B6_N")
#' beta_ordination(mice_B6_N, mycolor = "Genotype", myshape = "Experiment",
#' myformula = c("Genotype + Experiment + Cage.ID"))

beta_ordination <- function(phylo,
                            formula = "",
                            color,
                            distance_method = "bray",
                            ordination_type = "NMDS",
                            geom_star = TRUE,
                            permutations = 999,
                            parallel = 4,
                            shape = NULL,
                            tagtype = "A",
                            verbose = TRUE,
                            layout = c(area(1, 1), area(2, 1, 5))
                            ) {
  # sanity checks
  if (!is.numeric(permutations)) stop("permutations needs to be numeric")
  permutations <- as.character(permutations)

  if (!is.numeric(parallel)) stop("parallel needs to be numeric")
  permutations <- as.character(parallel)

  color <- c(color)
  shape <- c(shape)





  # Adonis test
  if (formula != "") {
    dist <- phyloseq::distance(phylo, distance_method)
    df <- as(phyloseq::sample_data(phylo), "data.frame")

    # because Adonis outputs the EXACT input parameters this string trick was used
    # to make a pleasent looking output
    myadonis <- paste("vegan::adonis(dist ~", formula,
                      ",df, permutations =", permutations,
                      ", parallel = ", parallel, ")")
    adonis_res <- eval(parse(text = myadonis))

    # prepare data for Adonis plot
    adonis_res_df <- data.frame(adonis_res$aov.tab)
    adonis_res_df[["model..id"]] <- row.names(adonis_res_df)
    adonis_res_df <- adonis_res_df[-NROW(adonis_res_df),]
    limits <- row.names(adonis_res_df[-NROW(adonis_res_df),])

    adonis_plot <- adonis_plot(adonis_res = adonis_res_df,
                               x = model..id,
                               value = R2,
                               ylab = "R2",
                               scale_x_discrete_limits = limits)
  }
  # Ordination plot
  phylo_ord <- phyloseq::ordinate(phylo, ordination_type, distance_method)

  ord_plot <- ordination_plot(phylo = phylo,
                              ordination = phylo_ord,
                              type = "samples",
                              color = color,
                              shape = shape,
                              geom_star = geom_star)

  if (formula != ""){
  # Plots adonis and ordination results together
  plots <- patchwork::wrap_plots(adonis_plot, ord_plot, design = layout) +
    patchwork::plot_annotation(tag_levels = tagtype)
  } else {
    plots <- ord_plot
  }

  plots
}

#' Adonis Plot
#'
#' Barplot for vegan::adonis results.
#'
#' @param adonis_res adonis object
#' @param x List of variables to map
#' @param value Result values to map
#' @param ylab label for y axis
#' @param xlab label for x axis
#' @param scale_x_discrete_limits
#'
#' @return ggplot
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @examples adonis_plot(adonis_res, ...)
adonis_plot <- function(adonis_res, x, value, ylab, xlab = "",
                        scale_x_discrete_limits) {
  x <- rlang::ensym(x)
  value <- rlang::ensym(value)

  adonis_plot <- ggplot(adonis_res, aes(x = !!x, y = !!value)) +
    ggplot2::geom_col(stat = "identity") +
    ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::scale_x_discrete(limits = scale_x_discrete_limits) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                                position = "right", limits = c(0.0, 1.0))
  adonis_plot
}

#' Plot Ordination
#'
#' @param phylo phyloseq object
#' @param ordination ordination object
#' @param type Plot type. Supported options are
#' c("samples", "sites", "species", "taxa", "biplot", "split", "scree")
#' @param color Map colors based on character string
#' of a sample variable or taxonomic rank
#' @param shape Map shapes based on character string of a sample variable
#' or taxonomic rank
#' @param geom_star boolean. If TRUE create a star plot by drawing segments from
#' group centroid to each points.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' ordination <- phyloseq::ordinate(mice_B6_N, "NMDS", "bray")
#' ordination_plot(mice_B6_N, ordination)
ordination_plot <- function(phylo, ordination,
                            type = "samples",
                            color = NULL,
                            shape = NULL,
                            geom_star = TRUE) {
  ggplot2::theme_set(ggplot2::theme_bw(16))
  ord_plot <- phyloseq::plot_ordination(phylo, ordination = ordination,
                                        type = type,
                                        color = color,
                                        shape = shape)
  ord_plot <- ord_plot + ggplot2::geom_point(size = 2)
    # + ggthemes::scale_colour_gdocs()
    if(geom_star) ord_plot <- ord_plot + ggpubr::stat_stars() + ggpubr::stat_mean(geom = "point", size = 2)
    ord_plot <- ord_plot +
    ggplot2::scale_shape_manual(values = c(16, 15, 17, 18:30))
  ggplot2::theme(legend.position = "bottom")

  ord_plot
}


#' Beta diversity distances inter/intra samples
#'
#' This function calculate beta diversity distances inter/intra samples and
#' returns a data frame
#' @param phylo phyloseq object
#' @param var variable in sample data of the phyloseq object, which should be
#' compared
#' @param dist Filter samples that do not reach this count number
#' @param sample_ID Sample identifier column
#' @keywords beta diversity
#' @export
#' @return data.frame
#' @examples
#' data("mice_B6_N")
#' get_beta_distances(mice_B6_N, var = "Genotype", dist = "unifrac")
#'
beta_dist_calc <- function(phylo, var, dist, sample_ID = "SampleID") {
  # calc distances
  phylo_dist <- phyloseq::distance(phylo, dist) # !!!second time distance is calculated. Get rid of it
  phylo_dist_df <- reshape2::melt(as.matrix(phylo_dist))

  # remove self-comparisons
  phylo_dist_df <- phylo_dist_df %>%
    dplyr::transmute(sample_a = as.character(Var1),
                     sample_b = as.character(Var2),
                     value = value) %>%
    dplyr::filter(sample_a != sample_b)

  # get sample data
  phylo_sample <- sample_data(phylo) %>%
    as_tibble() %>%
    dplyr::select(all_of(sample_ID), all_of(var))
  # Type conversion in S4 raises warnings, but gives expected results. Supressed
  phylo_sample <- suppressWarnings(
    dplyr::mutate_if(phylo_sample, is.factor, as.character)
  )

  # combined distances with sample data
  colnames(phylo_sample) <- c("sample_a", "type_a")
  phylo_dist_data <- dplyr::left_join(phylo_dist_df,
                                      phylo_sample, by = "sample_a")

  colnames(phylo_sample) <- c("sample_b", "type_b")
  phylo_dist_data <- phylo_dist_data %>%
    dplyr::left_join(phylo_sample, by = "sample_b") %>%
    dplyr::mutate(comparison = paste0(type_a, "_vs_", type_b))

  phylo_dist_data
}


#' Plot Beta diversity distances inter/intra samples
#'
#' This function
#' @param beta_distance data frame from beta_dist_calc
#' @param x x axis. Descriptive data.
#' @param y y axis. Value data.
#' @param jitter if TRUE adds addition geom_jitter
#' @param alpha only for jitter, regulates alpha value. Possible values from
#' 0 to 1
#' @keywords Beta diversity plot
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter
#'
#' @examples
#' data("mice_B6_N")
#' beta_dist <- beta_dist_calc(mice_B6_N, var = "Genotype", dist = "unifrac")
#' get_beta_distances_plot(beta_dist)

beta_dist_plot <- function(beta_distance, x = "comparison", y = "value",
                                    jitter = FALSE, alpha = 0.5) {
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  p <- ggplot(beta_distance, aes(x = !!x, y = !!y)) +
    geom_boxplot(aes(color = ifelse(type_a == type_b, "#80b1d3", "#fb8072")))
  if(jitter) {
    p <- p + geom_jitter(position = ggplot2::position_jitter(width = .1),
                         alpha = alpha,
                         aes(color = ifelse(type_a == type_b,
                                            "#80b1d3", "#fb8072")))
  }
  p <- p + ggplot2::scale_color_identity() +
    ggplot2::facet_wrap(~ type_a, scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Beta diversity") +
    ggplot2::ylab("distance")
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      legend.text = ggplot2::element_text(size = 22),
      strip.text.x = ggplot2::element_text(size = 22),
      axis.text.y = ggplot2::element_text(size = 22)
    )
  p
}

## To do add stat_compare_means
##  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

#' Calculate Abundance
#'
#' Calculate abundance in regards to taxonomic level.
#'
#' @param phylo phyloseq object
#' @param taxlev character. taxonomic level
#' @param ntaxa numeric. Number of most abundant taxa, which should be
#' displayed. Remaining taxa will be summarised in "Others".
#'
#' @return data.frame
#' @export
#'
#' @importFrom phyloseq prune_taxa taxa_sums otu_table
#'
#' @examples
#' data("mice_B6_N")
#' abundance_calc(mice_B6_N, "Family", 14)
abundance_calc <- function(phylo, taxlev, ntaxa) {
  # Subset OTUs with highest abundance for plotting
  phylo_glom <- phyloseq::tax_glom(phylo, taxlev)

  phylo_glom_names <- phylo_glom %>%
    phyloseq::taxa_sums() %>%
    sort(TRUE) %>%
    names()
  phylo_glom_top14 <- prune_taxa(phylo_glom_names[1:ntaxa], phylo_glom)

  # Adding taxon Others in the phyloseq object
  # calculating "others" values
  others_num <- abs(sample_sums(phylo_glom) - sample_sums(phylo_glom_top14))

  # create "Others" taxa
  others_otu_tbl <- rbind(others_num, otu_table(phylo_glom_top14))
  others_otu_tbl <- otu_table(others_otu_tbl, taxa_are_rows = TRUE)

  # adding "Others" taxa and taxonomy to the provided phyloseq object
  phylo_glom_top14@otu_table <- others_otu_tbl
  others_num <- rep("Others", ncol(tax_table(phylo_glom_top14)))

  phylo_glom_top14@tax_table <- phylo_glom_top14 %>%
    tax_table() %>%
    rbind(others_num, .) %>%
    tax_table()

  # make phyloseq object tidy because ggplot2 requires it
  melted_df <- phyloseq::psmelt(phylo_glom_top14)

  # order ranks
  taxlev <- rlang::ensym(taxlev)
  melted_df <- dplyr::arrange(melted_df, !!taxlev)
  melted_df[[taxlev]] <- forcats::fct_relevel(melted_df[[taxlev]],
                                              "Others",
                                              after = Inf)

  melted_df
}

#' Abundance Plot
#'
#' @param phylo Phyloseq object
#' @param x String Variable for x axis
#' @param taxlev String Taxonomic Level
#' @param orderbar String. Numeric Variable of sample_data. The Bars will be
#' ordered in respect to the values
#' @param facetby String Variable the plot should be faceted by
#' @param ntaxa Integer number of taxa which should be displayed
#' @param bar.params Additional parameters for the geom_bar plot
#' @param facet.params Additional parameters for facet_wrap
#' @param ...
#'
#' @return ggplot
#' @export
#'
#' @importFrom ggplot2 geom_bar facet_wrap aes ggplot
#'
#' @examples
#' data("mice_B6_N")
#' abundance_plot(mice_B6_N, "Microbiota")
abundance_plot <- function(phylo,
                            x,
                            taxlev = "Family",
                            orderbar = "",
                            facetby = "Description",
                            ntaxa = 14,
                            bar.params = list(),
                            facet.params = list(),
                           ...) {


  # fetch faulty user input
  if (ntaxa <= 0) stop("ntaxa must be higher than 0")
  if (as.integer(ntaxa) != ntaxa) stop("ntaxa must be an integer")

  # Calculate abundace for Plot
  melted_df <- abundance_calc(phylo, taxlev, ntaxa)


  # fetch user input to be able to unquote them later
  taxlev <- rlang::enexpr(taxlev)
  if(orderbar != "") orderbar_sym <- rlang::ensym(orderbar)

  colour_count <- length(unique(melted_df[[taxlev]]))
  taxlev <- rlang::ensym(taxlev)
  melted_df <- dplyr::arrange(melted_df, !!taxlev)
  melted_df[[taxlev]] <- forcats::fct_relevel(melted_df[[taxlev]], "Others", after = Inf)

  # handle custome plot options
  params <- list(...)
  bar.params <- modifyList(params, bar.params)
  facet.params <- modifyList(params, facet.params)

  bar <- do.call("geom_bar", modifyList(
    list(stat = "identity"), bar.params)
  )

  facet <- do.call("facet_wrap", modifyList(
    list(as.formula( paste("~", facetby)),
         scales = "free_x",
         labeller = ggplot2::label_wrap_gen(multi_line = FALSE)),
    facet.params)
  )

  x <- rlang::ensym(x)

  p <- ggplot(melted_df)
  if (orderbar == "") {
    p <- p + aes(x = !!x)
  } else {
    p <- p + aes(x = reorder(!!x, !!orderbar_sym))
  }
  p <- p + aes(y = Abundance, fill = !!taxlev) +
    ggplot2::geom_bar(stat = "identity") +
    facet +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 0.9,
                                                       size = 8)) +
    ggplot2::scale_fill_manual(values = ggpubr::get_palette(palette = "Paired",
                                                            k = colour_count)) +
    ggplot2::ylab("Abundance") +
    ggplot2::xlab(x)

  p
}
