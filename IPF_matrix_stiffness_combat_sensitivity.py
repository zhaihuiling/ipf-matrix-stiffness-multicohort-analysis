from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from pycombat import Combat


DATA = Path(r"C:\Users\ASUS\Documents\Codex\2026-10-13\zhe\work\data")
OUT = Path(r"C:\Users\ASUS\Documents\Codex\2026-10-21\12-combat-4-6-truncated-nfatc4\outputs")
UP = [
    "KIF20A", "IL7R", "HMMR", "KCNK2", "DLGAP5", "PRR11", "SERPINB7", "CDKN3",
    "ASPM", "CPA4", "VEPH1", "CENPI", "KNL1", "ANLN", "MKI67", "DEPDC1", "SHCBP1",
    "TTK", "DHCR24", "APOBEC3B",
]
DOWN = [
    "MMP10", "IL24", "PTGS2", "FOS", "MMP3", "SLC16A6", "NFATC2", "BMP2", "FOSB",
    "TNFAIP6", "CXCL8", "PLIN2", "SLC6A15", "ATP6V0D2", "ATP2A3", "ATF3", "FBXO32",
    "PDK4", "GPNMB", "KLHL24",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_gse150910() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(DATA / "GSE150910_counts.csv.gz")
    expr = raw.groupby("symbol", as_index=False).sum(numeric_only=True).set_index("symbol")
    lib = expr.sum(axis=0)
    expr = np.log2(expr.div(lib, axis=1) * 1_000_000 + 1)
    meta = pd.DataFrame({"sample": expr.columns})
    meta["group"] = meta["sample"].str.split("_", n=1).str[0].str.upper().replace({"CONTROL": "Control"})
    meta["batch"] = "GSE150910"
    meta["patient"] = meta["sample"]
    return expr, meta.set_index("sample")


def load_gse134692() -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(DATA / "GSE134692_log2cpm.txt.gz", sep="\t").set_index("ensembl")
    annotation = pd.read_csv(DATA / "GSE134692_gene_annotation.txt.gz", sep="\t")
    mapping = annotation.dropna(subset=["GeneName"]).drop_duplicates("ensembl").set_index("ensembl")["GeneName"]
    expr = expr.join(mapping, how="inner").groupby("GeneName").mean()
    design = pd.read_csv(DATA / "GSE134692_design.txt.gz", sep="\t").set_index("sample_id")
    design = design.loc[expr.columns].copy()
    donor = design["Sample.property"].astype(str)
    expr = expr.T.assign(_donor=donor).groupby("_donor").mean(numeric_only=True).T
    meta = design.reset_index().drop_duplicates("Sample.property").set_index("Sample.property").loc[expr.columns]
    meta = meta.rename_axis("sample").rename(columns={"DiseaseStatus": "group"})
    meta["batch"] = "GSE134692"
    meta["patient"] = meta.index
    return expr, meta[["group", "batch", "patient"]]


def load_gse124685() -> tuple[pd.DataFrame, pd.DataFrame]:
    import tarfile
    import xml.etree.ElementTree as ET

    raw = pd.read_csv(DATA / "GSE124685_fpkm.txt.gz", sep="\t")
    sample_cols = [c for c in raw.columns if c.startswith("GSM")]
    expr = np.log2(raw.set_index("gene_short_name")[sample_cols].groupby(level=0).mean() + 1)
    with tarfile.open(DATA / "GSE124685_family.xml.tgz") as archive:
        member = archive.getmembers()[0]
        root = ET.parse(archive.extractfile(member)).getroot()
    ns = {"m": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
    rows = []
    for sample in root.findall("m:Sample", ns):
        gsm = sample.attrib["iid"]
        title = sample.findtext("m:Title", namespaces=ns).strip()
        characteristics = {
            c.attrib.get("tag", ""): (c.text or "").strip()
            for c in sample.findall(".//m:Characteristics", ns)
        }
        group = characteristics.get("disease state", "")
        group = {"IPF": "IPF", "healthy donor": "Control"}.get(group, group)
        rows.append({"sample": gsm, "group": group, "patient": title.split(".", 1)[0]})
    meta = pd.DataFrame(rows).set_index("sample").loc[sample_cols]
    meta["batch"] = "GSE124685"
    return expr, meta[["group", "batch", "patient"]]


def gene_z_score(matrix: pd.DataFrame) -> pd.DataFrame:
    sd = matrix.std(axis=0, ddof=0).replace(0, np.nan)
    return matrix.subtract(matrix.mean(axis=0), axis=1).divide(sd, axis=1)


def score(matrix: pd.DataFrame, genes_up: list[str], genes_down: list[str]) -> pd.Series:
    z = gene_z_score(matrix)
    return z[genes_up].mean(axis=1) - z[genes_down].mean(axis=1)


def summary_scores(scores: pd.Series, meta: pd.DataFrame) -> pd.DataFrame:
    x = meta.copy()
    x["score"] = scores.reindex(x.index)
    rows = []
    comparisons = {
        "GSE150910_IPF_vs_Control": ("GSE150910", "IPF", "Control", "sample"),
        "GSE150910_CHP_vs_Control": ("GSE150910", "CHP", "Control", "sample"),
        "GSE134692_IPF_vs_Normal": ("GSE134692", "IPF", "Normal", "sample"),
        "GSE124685_IPF_vs_Control": ("GSE124685", "IPF", "Control", "patient"),
    }
    for label, (batch, case, control, unit) in comparisons.items():
        sub = x[x["batch"] == batch].copy()
        if unit == "patient":
            sub = sub.groupby(["patient", "group"], as_index=False)["score"].mean()
        a = sub.loc[sub["group"] == case, "score"]
        b = sub.loc[sub["group"] == control, "score"]
        rows.append({
            "comparison": label,
            "n_case": int(len(a)),
            "n_control": int(len(b)),
            "case_mean": float(a.mean()),
            "control_mean": float(b.mean()),
            "difference_case_minus_control": float(a.mean() - b.mean()),
            "unit": unit,
        })
    return pd.DataFrame(rows)


def component_concordance(
    matrix: pd.DataFrame,
    meta: pd.DataFrame,
    batch: str,
    case: str,
    control: str,
    unit: str,
    genes_up: list[str],
    genes_down: list[str],
) -> tuple[float, float]:
    keep = meta["batch"] == batch
    m = matrix.loc[keep, genes_up + genes_down]
    mm = meta.loc[keep].copy()
    if unit == "patient":
        m = m.groupby(mm["patient"]).mean()
        mm = mm.groupby("patient", as_index=True).first()
    z = gene_z_score(m)
    a = z.loc[mm["group"] == case].mean(axis=0)
    b = z.loc[mm["group"] == control].mean(axis=0)
    diff = a - b
    return float((diff.loc[genes_up] > 0).mean()), float((diff.loc[genes_down] < 0).mean())


def main() -> None:
    datasets = [load_gse150910(), load_gse134692(), load_gse124685()]
    common = set(UP + DOWN)
    for expr, _ in datasets:
        common &= set(expr.index)
    genes_up = [g for g in UP if g in common]
    genes_down = [g for g in DOWN if g in common]
    if len(genes_up) < 15 or len(genes_down) < 15:
        raise RuntimeError(f"Too few common signature genes: {len(genes_up)} UP, {len(genes_down)} DOWN")

    matrices = []
    metas = []
    input_files = [
        DATA / "GSE150910_counts.csv.gz",
        DATA / "GSE134692_log2cpm.txt.gz",
        DATA / "GSE134692_gene_annotation.txt.gz",
        DATA / "GSE134692_design.txt.gz",
        DATA / "GSE124685_fpkm.txt.gz",
        DATA / "GSE124685_family.xml.tgz",
    ]
    for expr, meta in datasets:
        matrices.append(expr.loc[sorted(common)])
        metas.append(meta)
    expr_all = pd.concat(matrices, axis=1).T
    meta_all = pd.concat(metas, axis=0).loc[expr_all.index]

    # Preserve broad biological group indicators while estimating batch effects.
    # Control is the reference category; CHP and ALI are retained as separate context indicators.
    design = pd.get_dummies(meta_all["group"], columns=["group"], dtype=float)
    design = design.reindex(columns=["IPF", "CHP", "ALI"], fill_value=0.0)
    batch = meta_all["batch"].to_numpy()

    before = summary_scores(score(expr_all, genes_up, genes_down), meta_all)
    combat = Combat(mode="p")
    corrected = combat.fit_transform(expr_all.to_numpy(dtype=float), batch, X=design.to_numpy(dtype=float))
    corrected_df = pd.DataFrame(corrected, index=expr_all.index, columns=expr_all.columns)
    after = summary_scores(score(corrected_df, genes_up, genes_down), meta_all)
    comparison = before.merge(after, on="comparison", suffixes=("_before", "_after"))
    comparison["direction_preserved"] = np.sign(comparison["difference_case_minus_control_before"]) == np.sign(comparison["difference_case_minus_control_after"])
    comparison["absolute_change_in_difference"] = (comparison["difference_case_minus_control_after"] - comparison["difference_case_minus_control_before"]).abs()
    concordance_specs = [
        ("GSE150910_IPF_vs_Control", "GSE150910", "IPF", "Control", "sample"),
        ("GSE150910_CHP_vs_Control", "GSE150910", "CHP", "Control", "sample"),
        ("GSE134692_IPF_vs_Normal", "GSE134692", "IPF", "Normal", "sample"),
        ("GSE124685_IPF_vs_Control", "GSE124685", "IPF", "Control", "patient"),
    ]
    concordance_rows = []
    for label, batch_name, case, control, unit in concordance_specs:
        up_b, down_b = component_concordance(expr_all, meta_all, batch_name, case, control, unit, genes_up, genes_down)
        up_a, down_a = component_concordance(corrected_df, meta_all, batch_name, case, control, unit, genes_up, genes_down)
        concordance_rows.append({
            "comparison": label,
            "UP_concordance_before": up_b,
            "DOWN_concordance_before": down_b,
            "UP_concordance_after": up_a,
            "DOWN_concordance_after": down_a,
        })
    concordance = pd.DataFrame(concordance_rows)
    comparison = comparison.merge(concordance, on="comparison", how="left")

    OUT.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(OUT / "ComBat_sensitivity_score_comparison.csv", index=False)
    corrected_df.to_csv(OUT / "ComBat_sensitivity_corrected_common_signature_matrix.csv")
    meta_all.reset_index(names="sample").to_csv(OUT / "ComBat_sensitivity_sample_manifest.csv", index=False)
    library_manifest = pd.read_csv(DATA / "GSE134692_design.txt.gz", sep="\t")
    library_manifest.to_csv(OUT / "GSE134692_library_to_donor_manifest.csv", index=False)
    meta_all.loc[meta_all["batch"] == "GSE124685"].reset_index(names="sample").to_csv(OUT / "GSE124685_sample_manifest.csv", index=False)
    result = {
        "analysis": "Exploratory ComBat sensitivity analysis across GSE150910, GSE134692, and GSE124685",
        "not_included": ["GSE22011: different in-vitro system", "GSE47460: no independent expression input available", "GSE135893: single-cell metadata rather than bulk expression"],
        "genes": {"common_total": len(common), "UP": genes_up, "DOWN": genes_down},
        "sample_counts": {"|".join(map(str, k)): int(v) for k, v in meta_all.groupby(["batch", "group"]).size().items()},
        "design_columns_preserved": list(design.columns),
        "direction_preserved_all_comparisons": bool(comparison["direction_preserved"].all()),
        "comparison": comparison.to_dict(orient="records"),
        "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "pycombat": "0.20"},
        "input_sha256": {str(p.name): sha256(p) for p in input_files},
        "interpretation": "Exploratory sensitivity analysis only; not independent validation. Batch and biological group are partially confounded, and GSE124685 contains repeated regions per patient.",
    }
    (OUT / "ComBat_sensitivity_audit.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
