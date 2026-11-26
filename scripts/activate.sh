#!/usr/bin/env bash

# Make CASTIE's source-tree command-line scripts available in `pixi shell`
# and to commands launched with `pixi run`.
export PATH="${PIXI_PROJECT_ROOT}/extdata:${PATH}"

# Skip the personal R library (a CRAN R's arm64 packages there crash on load); R ignores a missing path.
export R_LIBS_USER="${PIXI_PROJECT_ROOT}/.pixi/no-user-library"
