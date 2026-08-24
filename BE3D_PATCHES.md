# Portability patches applied to BE3D for these runs

We ran BE3D (`beclust3d`) on Windows / Python 3.14 / pandas 3.x. The upstream code needed a few small,
**non-scientific** portability fixes to run to completion. None change the LFC3D math, the null model, the
significance calculation, or the clustering. All are also captured with context in
[BE3D_IMPROVEMENTS.md](BE3D_IMPROVEMENTS.md). Upstream: https://github.com/broadinstitute/BE3D

1. **`beclust3d/aggregate/nonaggregate.py`** and **`.../metaaggregate.py`** — `to_csv(filename, "\t", ...)`
   passes `sep` positionally, which pandas 3.x rejected. Fix: `to_csv(filename, sep="\t", ...)`.
2. **`beclust3d/lfc3d/clustering_plot.py`** — significance labels like `p<0.05` were interpolated into output
   *filenames*; `<` and `>=` are illegal in Windows paths. Fix: sanitize only the path components
   (`<`→`lt`, `>=`→`ge`); the in-data labels are unchanged.
3. **`examples/be3d_local.py`** (`find_union`) — positional `Series[0]` indexing broke under label-based
   pandas. Fix: `list(series)[0]`.

### DSSP
`sequence_structural_features` shells out to `mkdssp`, which we did not have installed. DSSP only feeds the
final *characterization* step (not LFC3D / significance / clustering), so we supply a placeholder
`user_dssp` generated from the structure (`harness/gen_dssp.py`) to skip the `mkdssp` call. The
recommendation to upstream is to wrap that call in try/except or ship a pure-Python fallback (see the log).

### Recommended upstream fixes (beyond portability)
See [BE3D_IMPROVEMENTS.md](BE3D_IMPROVEMENTS.md) — highest-value: an `assay_direction`/`invert_score` flag
(BE-QA silently assumes a dropout sign, inverting activity/enrichment screens); FDR/q-values (the base rate
reaches 37–45% on some targets); pLDDT propagation into significance; a config-nameable control category;
and PPI/complex support for a structural-only (untiled) partner chain.
