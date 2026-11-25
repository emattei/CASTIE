# gc() unless options(CASTIE.gc = FALSE): step12_batch_qtl.R turns the explicit collections off, which cost
# ~1.5 s per gene in step 1 + step 2 and change no result.
castie_gc <- function() if (isTRUE(getOption("CASTIE.gc", TRUE))) invisible(gc())
