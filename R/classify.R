#' Random Forest Model
#'
#'Creates a classification model based on Random Forest
#'
#' @param phylo phyloseq object
#' @param classify_by variable with which the data should be classified by
#' @param minlib normalizing to minimum reads
#' @param prunescale minimum abundance in percantage an OTU must have to be considered
#' @param seed set seed for reproducable research
#' @param ntree number of trees to grow
#'
#' @return randomForest object
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' classify_randomForest(mice_B6_N, "treatment_group")
classify_randomForest <- function(phylo, classify_by, minlib = 1500, prunescale = 0.01, seed = 2, ntree = 500 ) {

  # Prune rare OTUs by mean relative abundance
  tax_mean <- phyloseq::taxa_sums(phylo)/nsamples(phylo)
  sites_prune <- phyloseq::prune_taxa(tax_mean > prunescale*minlib, phylo)

  # Make training data with OTUs as column and sample as row
  predictors <- t(phyloseq::otu_table(sites_prune))
  response <- as.factor(phyloseq::sample_data(sites_prune)[[classify_by]])
  rf_data <- data.frame(response,predictors)
  set.seed(seed)
  phylo_classfied <- randomForest::randomForest(response~., data = rf_data, ntree = ntree)

  return(phylo_classfied)
}

#' Plot Random Forest Results
#'
#' @param randomForest_obj randomForest object
#' @param top number of top OTUs, which should be plotted
#'
#' @return plot
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' prediction_model <- classify_randomForest(mice_B6_N, "treatment_group")
#' plot_randomForest_results(prediction_model, 10)
plot_randomForest_results <- function(randomForest_obj, phylo, classified_by, top = 20) {

  # Make a data frame with predictor names and their importance
  imp <- randomForest::importance(randomForest_obj)
  imp <- data.frame(OTU = rownames(imp), imp)

  # Order the predictor levels by importance
  imp_sort <- dplyr::arrange(imp, desc(MeanDecreaseGini))
  imp_sort$predictors <- factor(imp_sort$OTU, levels = imp_sort$OTU)

  imp_top20 <- imp_sort[1:top, ]

  # Find the top predictors for every variable
  phylo_tidy <- phyloseq_to_tidy(phylo)
  phylo_top20 <- dplyr::right_join(phylo_tidy, imp_top20, by = "OTU")

  p_simple <- ggplot2::ggplot(phylo_top20, ggplot2::aes(x = predictors, y = MeanDecreaseGini)) +
    ggplot2::geom_bar(stat = "identity", fill = "indianred") +
    ggplot2::coord_flip() +
    ggplot2::ggtitle("Most important OTUs for classification")

  p_complex <- ggplot2::ggplot(phylo_top20, ggplot2::aes(x = MeanDecreaseGini, y = OTU_abundance)) +
    ggplot2::geom_boxplot(ggplot2::aes(color = factor(OTU))) +
    ggplot2::facet_wrap(as.formula(paste("~",classified_by)),
                        scales = "free_y") +
    ggplot2::labs(colour = "OTU") +
    ggplot2::scale_color_discrete(labels = ggplot2::waiver())
  return(p_complex)
}
