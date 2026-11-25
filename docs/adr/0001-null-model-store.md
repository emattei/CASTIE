---
status: accepted
---

# Null models go through one store; batch mode can keep them in memory

The null model `.rda` was both the hand-off from step 1 to step 2 and step 1's working state: the fallback stages
rename and remove stage models by path, and the theta check and report look for the file. In one-session batch
mode every feature still wrote ~59 MB that was deleted unread at the end of the session. We decided that every
save, load, existence check, rename and removal of a null model goes through `R/model_cache.R`, keyed by path,
and that writing to disk is a setting of that store: `step12_batch_qtl.R --keepModels=FALSE` keeps models in
memory only, and the runner empties the store at each new feature.

## Considered Options

- Pass the model object from step 1 to step 2 in `step12_batch_qtl.R`: changes both scripts' interfaces and
  leaves step 1's fallback stages working on files.
- Keep writing files and delete each one right after its step 2: saves disk space held but not bytes written.

## Consequences

- Nothing outside `R/model_cache.R` touches a null model file (`file.exists`, `file.rename`, `load`, ...), except the
  LOCO saves in `fitNULLGLMM_multiV`, which write their own files; a path the store does not hold is read from file.
- `MODEL_FROM`, re-running step 2, or inspecting a model needs `--keepModels=TRUE` (the default) or standalone step 1.
- Shrinking the saved model (no `XVX_inv_XV`; `X` stored instead of `XV` and `XXVX_inv`, rebuilt by `readInGLMM`)
  means a step 2 from an older CASTIE cannot read new models; a newer step 2 still reads old ones.
