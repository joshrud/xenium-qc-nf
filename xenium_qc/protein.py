"""Protein QC for Xenium Protein runs: detection, RNA-protein concordance, lineage exclusivity."""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

from .common import Bundle, StepOutput, load_markers, setup_logging, timed


def run(bundle: str, outdir: str, cells_qc: str, markers: str | None = None, n_sample: int = 100_000):
    log = setup_logging("protein", outdir)
    b = Bundle(bundle)
    out = StepOutput(outdir, "protein", log)
    M = out.metrics
    cfg = load_markers(markers)
    raw = sc.read_10x_h5(b.file("cell_feature_matrix.h5"), gex_only=False)
    ftype = raw.var["feature_types"]
    if not (ftype == "Protein Expression").any():
        log.info("no Protein Expression features; protein QC skipped")
        M["skipped"] = True
        out.write_metrics()
        return
    qc = pd.read_parquet(cells_qc).set_index("cell_id")
    keep = qc.index[qc.pass_qc]
    prot = raw[raw.obs_names.isin(keep), ftype == "Protein Expression"]
    rna = raw[raw.obs_names.isin(keep), ftype == "Gene Expression"].copy()
    rna.var_names_make_unique()
    P = pd.DataFrame(prot.X.toarray() if sp.issparse(prot.X) else prot.X, index=prot.obs_names, columns=prot.var_names)
    log.info("%d protein features on %d QC-passing cells", P.shape[1], P.shape[0])

    with timed(log, "protein detection"):
        ps = pd.DataFrame({"frac_positive": (P > 0).mean(), "median_if_positive": P.where(P > 0).median(),
                           "p99": P.quantile(0.99)}).sort_values("frac_positive")
        out.table(ps.rename_axis("protein").reset_index(), "qc_protein_summary")
        M["protein_frac_positive"] = ps.frac_positive.round(3).to_dict()

    with timed(log, "RNA-protein concordance"):
        pairs = {p: g for p, g in cfg.get("protein_gene_pairs", {}).items() if p in P.columns and g in rna.var_names}
        idx = P.sample(min(n_sample, len(P)), random_state=0).index
        R = rna[idx, list(pairs.values())].to_df()
        conc = pd.DataFrame([{"protein": p, "gene": g, "spearman_rho": stats.spearmanr(P.loc[idx, p], R[g]).statistic,
                              "frac_rna_pos": (R[g] > 0).mean(), "frac_prot_pos": (P.loc[idx, p] > 0).mean()}
                             for p, g in pairs.items()]).sort_values("spearman_rho", ascending=False)
        out.table(conc, "qc_rna_protein_concordance")
        M["weak_concordance"] = conc[conc.spearman_rho < 0.15].protein.tolist()
        fig, ax = plt.subplots(figsize=(6, 0.25 * len(conc) + 1))
        ax.barh(conc.protein + " / " + conc.gene, conc.spearman_rho)
        ax.axvline(0.15, color="k", ls="--", lw=1); ax.invert_yaxis(); ax.set_xlabel("Spearman ρ (per cell)")
        ax.set_title("RNA–protein concordance")
        out.figure(fig, "qc_rna_protein_concordance")

    with timed(log, "protein lineage exclusivity"):
        rows = []
        for a, c in [tuple(x) for x in cfg.get("protein_exclusive", []) if all(p in P.columns for p in x)]:
            ta, tc = P[a][P[a] > 0].quantile(0.9), P[c][P[c] > 0].quantile(0.9)
            obs_rate, bg = (P.loc[P[a] >= ta, c] >= tc).mean(), (P[c] >= tc).mean()
            rows.append({"high_in": a, "also_high": c, "observed": obs_rate, "background_rate": bg,
                         "enrichment": obs_rate / bg if bg else None})
        if rows:
            out.table(pd.DataFrame(rows), "qc_protein_exclusivity")
    out.write_metrics()
