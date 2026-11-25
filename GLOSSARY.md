# CASTIE

Single-cell genetic association testing: a per-cell count model per feature, with donors as the unit of genotype and contexts along which genetic effects may vary.

## Language

**Feature**:
The molecular unit whose per-cell counts are tested: a gene for eQTLs, a peak for caQTLs.
_Avoid_: phenotype, trait (when meaning the unit itself)

**Peak**:
A feature defined as a genomic interval of open chromatin, whose per-cell value is its fragment count.
_Avoid_: region, bin

**Cis window**:
The genomic range around a feature whose variants are tested against it; it changes which variants are tested, never the feature's counts.
_Avoid_: region (alone), flank

**Donor**:
The genetically distinct individual that cells come from; the unit that carries a genotype.
_Avoid_: sample, individual

**Library**:
A batch of cells processed and sequenced together, containing many donors; a cell-level technical covariate.
_Avoid_: pool, sample, PKR

**Context**:
A per-cell numeric value along which a genetic effect is allowed to vary; the GxC test asks whether the effect changes with it.
_Avoid_: dynamic covariate, condition, state (unless naming what the context encodes)

**Reference context**:
The category left out when a categorical context is encoded as 0/1 columns; GxC effects are relative to it.
_Avoid_: baseline

**Program**:
A gene module whose per-cell activity could serve as a context, but is not one until chosen as such.
_Avoid_: context (when no choice has been made), factor, topic

**Detection rate**:
The fraction of cells (or of donors) with a non-zero count for a feature.
_Avoid_: sparsity, coverage

**Null model**:
The per-feature count model fitted without genotypes; each variant in the feature's cis window is tested against it.
_Avoid_: step 1 model, GLMM (alone), .rda
