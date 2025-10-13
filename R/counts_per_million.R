#' Generating counts per million from a SingleCellExperiment object
#'
#' @param sce A SingleCellExperiment object
#' @param input_file A character vector of counts Anndata file path
#' @param output_file A character vector of CPM Anndata file path
#' @return A directory stores counts per million Anndata
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay assays assays<-
#' @importFrom purrr map
#' @noRd
get_counts_per_million <- function(sce, input_file, output_file) {

  # Save SCE to the cache directory counts folder
  sce |> zellkonverter::writeH5AD(input_file, compression = "gzip")
  
  # Avoid completely empty cells
  col_sums <- colSums(as.matrix(assay(sce)))
  selected_cols <- which(col_sums >0 & col_sums < Inf)
  
  assay_name <- assays(sce) |> names()
  sce <- SingleCellExperiment(list(cpm = scuttle::calculateCPM(sce[,selected_cols ,drop=FALSE ], assay.type = assay_name)))
  rownames(sce) <- rownames(sce[,selected_cols  ])
  colnames(sce) <- colnames(sce[,selected_cols  ])
  
  # Avoid scaling zeros
  sce_zero <- SingleCellExperiment(list(cpm = assay(sce)[, !selected_cols ,drop=FALSE ]))
  rownames(sce_zero) <- rownames(sce[, !selected_cols  ])
  colnames(sce_zero) <- colnames(sce[, !selected_cols ])
  
  sce <- sce |> cbind(sce_zero)
  
  sce <- sce[,colnames(sce)]
  
  # Check if there is a memory issue 
  assays(sce) <- assays(sce) |> map(DelayedArray::realize)
  
  sce |> zellkonverter::writeH5AD(output_file, compression = "gzip")
} 

