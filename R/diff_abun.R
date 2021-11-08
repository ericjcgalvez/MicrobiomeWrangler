#' Differential expression analysis
#'
#' Test for differential expression based on
#' negative binomal distribution
#'
#' @param phylo Phyloseq Object with OTU table and sample data
#' @param var Variable of sample data to analyze
#' @param test significance test used in DESeq2. Either "Wald" or "LRT"
#' @param fitType type of fitting of dispersions to the mean intensity
#' in DESeq2. Either "parametric", "local" or "mean"
#'
#' @return DESeqDataSet object including data set and test results
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' desq2_analysis(mice_B6_N, "Vendor")
deseq2_analysis <- function(phylo, var, alpha = 0.01, test = "Wald", fitType = "parametric") {

  diagdds <-phyloseq::phyloseq_to_deseq2(phylo, as.formula(paste("~",var)))

  # Performing Tests
  diagdds <- DESeq2::DESeq(diagdds, test = test, fitType = fitType)

  # Formatting results
  res = DESeq2::results(diagdds, cooksCutoff = FALSE)
  # sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(res, "data.frame"), as(tax_table(phylo)[rownames(res), ], "matrix"))

  return(sigtab)
}

#' Plot DESeq2 analysis
#'
#' Plots results from DESeq2 analysis
#'
#' @param deseq2_result deseq2 object
#' @param title Title of the plot
#' @param label_taxlev Taxonomic label for the top hits
#' @param top_hits Top hits which should be labled
#' @param fdr single numeric. false discovery rate
#' @param fc single numeric. fold change
#' @param size single numeric
#' @param alpha single numeric.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' deseq2_res <- deseq2_analysis(mice_B6_N, "Mircobiota")
#' deseq2_plot(deseq2_res, "Deseq", "Family")
deseq2_plot <- function(deseq2_result, title, label_taxlev, top_hits = 10,
                        fdr = 0.05, fc = 2, size = 2, alpha = 0.3) {
  p <- ggpubr::ggmaplot(deseq2_result, main = title,
                        fdr = fdr, fc = fc, size = size, alpha = alpha,
                        palette = c("#B31B21", "#3366cc", "darkgray"),
                        genenames = as.vector(deseq2_result[[label_taxlev]]),
                        legend = "top", top = top_hits,
                        label.rectangle = F,
                        font.label = c("bold", 14),
                        font.legend = "bold",
                        font.main = "bold",
                        ggtheme = ggplot2::theme_classic(20))
  return(p)
}

#' Join ALDEx2 results with taxonomic information
#'
#' @param phylo phyloseq object
#' @param aldex2_res aldex2 object
#'
#' @return tibble
#' @export
#' @example data("mice_B6_N")
#' aldex2_res <- ALDEx2::aldex(mice_B6_N@@otu_table,
#' as.character(sample_data(phylo)$Vendor))
#' aldex2_res_sig(mice_B6_N, aldex2_res)
#'
aldex2_join_tax <- function(phylo, aldex2_res) {
  aldex2_res <- tibble::rownames_to_column(aldex2_res)
  phylo_tax <- tibble::rownames_to_column(data.frame(phylo@tax_table))
  sig_tax <- dplyr::left_join(phylo_tax, aldex2_res)
  sig_tax <- tibble::remove_rownames(sig_tax)
  sig_tax <- tibble::column_to_rownames(sig_tax, var = "rowname")

  return(sig_tax)
}

#' ALDEx2 Significant OTUs
#'
#' Format ALDEx2 results returning taxonomic
#' information for significant results
#'
#' @param phylo phyloseq object
#' @param aldex2_res aldex2 object
#' @param sig single numeric. Significance threshold
#' @param all_sig identify which values are siginifcant in both t-test and
#' glm test
#'
#' @return data.frame
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' aldex2_res <- ALDEx2::aldex(mice_B6_N@@otu_table,
#' as.character(sample_data(phylo)$Vendor))
#' aldex2_res_sig(mice_B6_N, aldex2_res)
aldex2_res_sig <- function(phylo, aldex2_res, sig = 0.05, all_sig = TRUE) {
  sig_tax <- aldex2_join_tax(phylo, aldex2_res)

  if(all_sig == TRUE) {
    # identify which values are siginifcant in both t-test and glm test
    aldex2_sig <- sig_tax[which(sig_tax$we.eBH < sig & sig_tax$wi.eBH < sig),]
  } else if (all_sig == FALSE) {
    # identify which values are significant in at least one test
    aldex2_sig <- sig_tax[which(sig_tax$we.eBH < sig | sig_tax$wi.eBH < sig),]
  } else {
    message("The parameter all_sig must be a boolean")
  }

  return(aldex2_sig)
}

#' Compare Differential Abundance Analysis
#'
#' Returns matching or missmatching rows from the
#' Differential Abundance Analysis of different methods
#'
#' @param x data.frame with significant observations as rownames
#' @param y data.frame with significant observations as rownames
#' @param by identifier in both tables
#'
#' @return list
#' @export
#'
#' @examples
#' data("mice_B6_N")
#' deseq2_res <- deseq2_analysis(mice_B6_N, "Mircobiota")
#' aldex2_res <- ALDEx2::aldex(mice_B6_N@@otu_table,
#' as.character(sample_data(phylo)$Vendor))
#' comp_diff_abund_analysis(deseq2_res, aldex2_res, type = "match")
comp_diff_abund_analysis <- function(x,y, type, by = "OTU") {
  x_tidy <- tibble::rownames_to_column(x, var = by)
  y_tidy <- tibble::rownames_to_column(y, var = by)

  if(type == "match"){
    shared_sig <- dplyr::semi_join(x_tidy, y_tidy, by=by)
    return(shared_sig)
  } else if (type == "mismatch"){
    diff_sig <- x_tidy[!x_tidy[[by]] %in% y_tidy[[by]],]
    diff_sig2 <- y_tidy[!y_tidy[[by]] %in% x_tidy[[by]],]

    return(list(in_x_not_y = diff_sig, in_y_not_x = diff_sig2))
  } else {
    stop("type must be either 'match' or 'mismatch'")
  }
}
