.model_cache <- new.env(parent = emptyenv()) # options(CASTIE.cache = TRUE): the model last saved by castie_save_model

.model_key <- function(file) {
  i <- file.info(file, extra_cols = FALSE)
  paste(normalizePath(file), i$size, as.numeric(i$mtime))
}

#' Save a null model as `save(modglmm, file = file, compress = FALSE)`. With options(CASTIE.cache = TRUE) the
#' object is also kept, so castie_load_model() can return it without reading the file back.
#' @param modglmm the null model object.
#' @param file the .rda path.
castie_save_model <- function(modglmm, file) {
  save(modglmm, file = file, compress = FALSE)
  if (isTRUE(getOption("CASTIE.cache"))) list2env(list(key = .model_key(file), model = modglmm), envir = .model_cache)
  invisible(file)
}

#' The `modglmm` stored in an .rda file: the kept object when it is what castie_save_model() last wrote to that
#' file (same path, size and modification time), otherwise read with load().
#' @param file the .rda path.
castie_load_model <- function(file) {
  if (isTRUE(getOption("CASTIE.cache")) && identical(.model_cache$key, .model_key(file))) {
    return(.model_cache$model)
  }
  e <- new.env()
  load(file, envir = e)
  e$modglmm
}
