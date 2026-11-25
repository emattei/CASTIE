.pheno_cache <- new.env(parent = emptyenv()) # read_pheno_h5, options(CASTIE.cache = TRUE)

#' Read a phenotype HDF5 file produced by `convert_phenoFile_to_h5.R`.
#'
#' Reads metadata columns named in `checkColList` from /meta and one gene
#' column (named in `phenoCol`) from the CSC sparse matrix under /counts.
#' Returns a data.frame with one row per cell.
#'
#' Schema:
#'   /meta/<col>           - one dataset per metadata column
#'   /counts/data          - non-zero values
#'   /counts/indices       - row indices of non-zeros (0-based)
#'   /counts/indptr        - column pointers (length ngenes+1, double)
#'   /counts/gene_names    - gene names (one per column)
read_pheno_h5 <- function(phenoFile, phenoCol, sampleIDColinphenoFile, checkColList) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r package is required to read HDF5 phenotype files. Install with: install.packages('hdf5r')")
  }
  h5 <- hdf5r::H5File$new(phenoFile, mode = "r")
  on.exit(h5$close(), add = TRUE)

  metaNames <- names(h5[["meta"]])
  metaCols <- intersect(checkColList, metaNames)
  # options(CASTIE.cache = TRUE) (step12_batch_qtl.R): the metadata, gene names and indptr are the same for
  # every gene of a file, so they are read once per session; only the gene's counts are read each time.
  key <- paste(normalizePath(phenoFile), file.mtime(phenoFile), paste(metaCols, collapse = ","))
  if (isTRUE(getOption("CASTIE.cache")) && identical(.pheno_cache$key, key)) {
    data <- .pheno_cache$data
    geneNames <- .pheno_cache$geneNames
    indptr_raw <- .pheno_cache$indptr
  } else {
    data <- data.frame(lapply(metaCols, function(col) h5[["meta"]][[col]]$read()),
      stringsAsFactors = FALSE
    )
    colnames(data) <- metaCols
    geneNames <- h5[["counts"]][["gene_names"]]$read()
    ## indptr is stored as double to support nnz > 2^31; keep as numeric for
    ## indexing (R doubles give 53-bit integers, ample for any realistic dataset)
    indptr_raw <- h5[["counts"]][["indptr"]]$read()
    if (isTRUE(getOption("CASTIE.cache"))) {
      list2env(list(key = key, data = data, geneNames = geneNames, indptr = indptr_raw), envir = .pheno_cache)
    }
  }

  geneIdx <- match(phenoCol, geneNames)
  if (is.na(geneIdx)) {
    stop("ERROR! phenoCol '", phenoCol, "' not found in HDF5 counts/gene_names\n")
  }

  ptr_start <- indptr_raw[geneIdx] + 1 # 0-based -> 1-based
  ptr_end <- indptr_raw[geneIdx + 1L]
  nrows <- nrow(data)
  geneVec <- rep(0L, nrows)
  if (ptr_end >= ptr_start) {
    ds_idx <- h5[["counts"]][["indices"]]
    ds_val <- h5[["counts"]][["data"]]
    idx <- ds_idx[ptr_start:ptr_end]
    val <- ds_val[ptr_start:ptr_end]
    geneVec[idx + 1L] <- as.integer(val)
  }
  data[[phenoCol]] <- geneVec

  if (sampleIDColinphenoFile %in% colnames(data)) {
    data[[sampleIDColinphenoFile]] <- as.character(data[[sampleIDColinphenoFile]])
  }
  cat(sprintf(
    "Read phenotype from HDF5: %d rows, gene '%s' (%d non-zeros)\n",
    nrows, phenoCol, sum(geneVec != 0)
  ))
  data
}
