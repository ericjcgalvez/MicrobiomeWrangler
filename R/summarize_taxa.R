summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    plyr::alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    BiocGenerics::tapply(counts, taxonomy, sum)
  }
}


phyloseq_summarize_taxa <- function(psdata, taxonomic.ranks = rank_names(psdata)[2,6]) {
  if(length(taxonomic.ranks) > 1) {
    names(taxonomic.ranks) <- taxonomic.ranks
    plyr::llply(taxonomic.ranks, phyloseq_summarize_taxa, psdata = psdata)
  } else {
    taxa <- as(tax_table(psdata)[, taxonomic.ranks], 'character')
    sum_tax_table <- summarize_taxa(as(otu_table(psdata), 'matrix'), taxa)
    phyloseq(otu_table(sum_tax_table, taxa_are_rows = TRUE), sample_data(psdata, FALSE))

  }
}

# summarize_taxa Generic function
## modified from https://github.com/joey711/phyloseq/issues/418
