#!/usr/bin/env Rscript
# Run step 1 and step 2 for many genes in one R session, so R starts and loads CASTIE once per batch
# instead of twice per gene. Each line of --batch is
#   <gene> TAB <step: 1 or 2> TAB <log file> TAB <arguments separated by \x1f>
# and step1_fitNULLGLMM_qtl.R / step2_tests_qtl.R (found on PATH) are sourced with those arguments,
# exactly as on the command line. One line per run goes to --results:
#   <gene> TAB <step> TAB <seconds> TAB <"ok" or the error message>
# Step 2 of a gene whose step 1 failed is skipped.
suppressPackageStartupMessages({
  library(CASTIE)
  library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--batch", type = "character", help = "Required. The runs, one per line (see above)"),
  make_option("--results", type = "character", help = "Required. Where to append one line per run")
)))
script <- c("1" = unname(Sys.which("step1_fitNULLGLMM_qtl.R")), "2" = unname(Sys.which("step2_tests_qtl.R")))
if (!all(nzchar(script))) stop("step1_fitNULLGLMM_qtl.R and step2_tests_qtl.R must be on PATH")
options(CASTIE.gc = FALSE, CASTIE.cache = TRUE) # no explicit gc(); read the phenotype metadata and transform the covariates once

failed <- character(0)
for (ln in readLines(opt$batch)) {
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  gene <- f[1]
  step <- f[2]
  if (step == "2" && gene %in% failed) next
  con <- file(f[3], open = "a")
  sink(con)
  sink(con, type = "message")
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    {
      options(CASTIE.args = strsplit(f[4], "\x1f", fixed = TRUE)[[1]])
      source(script[[step]], local = new.env(parent = globalenv()))
      "ok"
    },
    error = function(e) {
      message("Error: ", conditionMessage(e)) # into the gene's log
      gsub("[\t\n]", " ", conditionMessage(e))
    }
  )
  options(CASTIE.args = NULL)
  sink(type = "message")
  sink()
  close(con)
  if (res != "ok") failed <- c(failed, gene)
  cat(gene, step, round(proc.time()[["elapsed"]] - t0), res, sep = "\t", file = opt$results, append = TRUE)
  cat("\n", file = opt$results, append = TRUE)
}
