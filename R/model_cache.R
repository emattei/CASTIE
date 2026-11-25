# The null model store: only these functions save, load, check, rename or remove null models, keyed by path.
# With options(CASTIE.keepModels = FALSE) (step12_batch_qtl.R --keepModels=FALSE) they act on .model_store, not on
# files: no null model file is created, renamed or removed. Paths not in .model_store are still read from file.
.model_cache <- new.env(parent = emptyenv()) # options(CASTIE.cache = TRUE): the model last saved by castie_save_model
.model_store <- new.env(parent = emptyenv()) # keepModels = FALSE: path -> null model, NULL if reserved, FALSE removed

.in_memory <- function() isFALSE(getOption("CASTIE.keepModels"))
.store_key <- function(file) file.path(normalizePath(dirname(file), mustWork = FALSE), basename(file)) # no file to normalise
.in_store <- function(file) exists(.store_key(file), envir = .model_store, inherits = FALSE)

.model_key <- function(file) {
  i <- file.info(file, extra_cols = FALSE)
  paste(normalizePath(file), i$size, as.numeric(i$mtime))
}

#' Save a null model as `save(modglmm, file = file, compress = FALSE)`. With options(CASTIE.cache = TRUE) the
#' object is also kept, so castie_load_model() can return it without reading the file back.
#' @param modglmm the null model object.
#' @param file the .rda path.
castie_save_model <- function(modglmm, file) {
  if (.in_memory()) return(invisible(assign(.store_key(file), modglmm, envir = .model_store)))
  save(modglmm, file = file, compress = FALSE)
  if (isTRUE(getOption("CASTIE.cache"))) list2env(list(key = .model_key(file), model = modglmm), envir = .model_cache)
  invisible(file)
}

#' The `modglmm` stored in an .rda file: the kept object when it is what castie_save_model() last wrote to that
#' file (same path, size and modification time), otherwise read with load().
#' @param file the .rda path.
castie_load_model <- function(file) {
  if (.in_memory() && isFALSE(.model_store[[.store_key(file)]])) stop("null model ", file, " was removed")
  if (.in_memory() && !is.null(m <- .model_store[[.store_key(file)]])) return(m) # else the file: MODEL_FROM, low-memory LOCO
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
castie_model_exists <- function(file, reserved = FALSE) {
  if (.in_memory() && .in_store(file)) {
    m <- .model_store[[.store_key(file)]]
    return(!isFALSE(m) && (reserved || !is.null(m)))
  }
  file.exists(file) && (reserved || file.size(file) > 0)
}

#' Create the empty file a null model will be saved to, so an unwritable path warns before the fit.
#' @param file the .rda path.
castie_reserve_model <- function(file) {
  if (.in_memory()) return(invisible(assign(.store_key(file), NULL, envir = .model_store)))
  file.create(file, showWarnings = TRUE)
}

#' Move the null model at `from`, or the empty file reserved there, to `to`. Nothing happens if there is neither.
#' @param from,to the .rda paths.
castie_rename_model <- function(from, to) {
  if (!.in_memory()) return(if (file.exists(from)) file.rename(from, to))
  if (.in_store(from) && castie_model_exists(from, reserved = TRUE)) {
    assign(.store_key(to), .model_store[[.store_key(from)]], envir = .model_store)
    rm(list = .store_key(from), envir = .model_store)
  }
}

#' Remove the null model at `file`, or the empty file reserved there.
#' @param file the .rda path.
castie_remove_model <- function(file) {
  if (!.in_memory()) return(if (file.exists(file)) file.remove(file))
  assign(.store_key(file), FALSE, envir = .model_store) # so a stale file left at this path is not read
}

#' Forget every null model kept in memory (CASTIE.keepModels = FALSE): step12_batch_qtl.R calls it at each new feature.
castie_clear_models <- function() rm(list = ls(.model_store, all.names = TRUE), envir = .model_store)
