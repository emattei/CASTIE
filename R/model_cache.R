# The null model store: only these functions save, load, check, rename or remove null models, keyed by path.
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

#' Whether a null model is stored at `file`. The empty file castie_reserve_model() leaves is not one.
#' @param file the .rda path.
#' @param reserved TRUE to also count the empty file castie_reserve_model() leaves.
castie_model_exists <- function(file, reserved = FALSE) file.exists(file) && (reserved || file.size(file) > 0)

#' Create the empty file a null model will be saved to, so an unwritable path warns before the fit.
#' @param file the .rda path.
castie_reserve_model <- function(file) file.create(file, showWarnings = TRUE)

#' Move the null model at `from`, or the empty file reserved there, to `to`. Nothing happens if there is neither.
#' @param from,to the .rda paths.
castie_rename_model <- function(from, to) if (file.exists(from)) file.rename(from, to)

#' Remove the null model at `file`, or the empty file reserved there.
#' @param file the .rda path.
castie_remove_model <- function(file) if (file.exists(file)) file.remove(file)
