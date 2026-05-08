"""Loaders for the supplementary datasets that ship in
``<ZENODO>/DATASETS/DatasetS{1..8}.xlsx``.

Each dataset has a "Description" sheet plus one or more data sheets.
We read the XLSX files directly (rather than the matching .txt exports,
which mix encodings) and apply a small per-sheet layout hint so the
result is a clean DataFrame even when the XLSX has banner rows or
nested per-cluster blocks.

Layout hints
------------
* ``flat``     — row 0 is a banner sentence, row 1 is a header, then data.
* ``no_banner``— row 0 is the header, then data.
* ``list``     — a single column of IDs (no banner, no header).
* ``raw``      — leave as-is (used for the per-cluster nested blocks
                 in S2 A/E and S8 A — these don't fit a single rectangular
                 schema; the user gets the raw, banner-stripped sheet).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st


# ---------------------------------------------------------------------------
# Registry: which sheets in which file, and how to lay them out.
# ---------------------------------------------------------------------------

@dataclass
class SheetSpec:
    label: str       # nice name shown in the UI
    sheet: str       # exact XLSX sheet name
    layout: str      # 'flat' | 'no_banner' | 'list' | 'raw'
    note: str = ""   # short caption shown above the table


@dataclass
class DatasetSpec:
    dataset_id: str    # e.g. 'S1'
    file: str          # 'DatasetS1.xlsx'
    title: str
    sheets: List[SheetSpec]


REGISTRY: Dict[str, DatasetSpec] = {
    "S1": DatasetSpec(
        dataset_id="S1",
        file="DatasetS1.xlsx",
        title="S1 — Features used for evolutionary signatures",
        sheets=[
            SheetSpec("A) Zarin et al.", "A) Zarin et al.", "flat",
                      "51 features overlapping with Zarin et al. eLife 2019."),
            SheetSpec("B) Additional features", "B) Additional features", "flat",
                      "93 additional human-IDR-ome features."),
        ],
    ),
    "S2": DatasetSpec(
        dataset_id="S2",
        file="DatasetS2.xlsx",
        title="S2 — Overrepresented GO terms & features in clusters",
        sheets=[
            SheetSpec("A) Selected clusters", "A) Selected clusters", "raw",
                      "Manually selected clusters across thresholds — nested per-cluster GO + feature blocks."),
            SheetSpec("B) IDRs in selected clusters", "B) IDRs in selected clusters", "flat",
                      "IDR → cluster annotations for the manually selected clusters in Tab A."),
            SheetSpec("C) GO terms classification", "C) GO terms classification", "flat",
                      "Higher-level functional categories assigned to selected clusters."),
            SheetSpec("D) Features over categories", "D) Features over categories", "flat",
                      "Qualitative feature → function classification."),
            SheetSpec("E) Automatic clusters (0.7)", "E) Automatic clusters", "raw",
                      "Automatic clusters at threshold 0.7 — nested per-cluster blocks."),
            SheetSpec("F) IDRs in automatic clusters", "F) IDRs in automatic clusters", "flat",
                      "IDR → automatic-cluster annotations at threshold 0.7."),
            SheetSpec("G) Automatic cluster GO O/E", "G) Automatic cluster GO O_E", "flat",
                      "Overrepresentation ratios for GO terms across the 393 automatic clusters."),
            SheetSpec("H) Auto clusters GO coverage", "H) Auto clusters GO coverage", "flat",
                      "Per-threshold counts and % of clusters with overrepresented GO terms."),
        ],
    ),
    "S3": DatasetSpec(
        dataset_id="S3",
        file="DatasetS3.xlsx",
        title="S3 — IDRs in clusters featured in Figure 1",
        sheets=[
            SheetSpec("A) Cluster 2133 / 0.4 — histone & chromatin binding",
                      "A) Cluster 2133 0.4", "list"),
            SheetSpec("B) Cluster 497 / 0.5 — high Glu, low Ser",
                      "B) Cluster 497 0.5", "list"),
            SheetSpec("C) Cluster 2403 / 0.5 — Cys-rich, S/G repeats",
                      "C) Cluster 2403 0.5", "list"),
            SheetSpec("D) Cluster 3046 / 0.5 — Gly/Arg motifs, SH3-binding",
                      "D) Cluster 3046 0.5", "list"),
            SheetSpec("E) Cluster 277 / 0.7 — nuclear pore, signal-seq binding",
                      "E) Cluster 277 0.7", "list"),
            SheetSpec("F) Cluster 326 / 0.7 — mRNA processing, RNA splicing",
                      "F) Cluster 326 0.7", "list"),
        ],
    ),
    "S4": DatasetSpec(
        dataset_id="S4",
        file="DatasetS4.xlsx",
        title="S4 — IDPs in the proteome and on the IDR-ome map",
        sheets=[
            SheetSpec("All IDPs (≥95% disordered)", "IDPs", "list",
                      "718 UniProt IDs predicted >=95% disordered by SPOT-Disorder v1."),
            SheetSpec("Clustered IDPs", "clustered IDPs", "list",
                      "492 IDPs that land in clusters of the IDR-ome map."),
            SheetSpec("Clustered IDPs in GO-overrepresented clusters",
                      "clustered IDPs with GO terms", "list"),
        ],
    ),
    "S5": DatasetSpec(
        dataset_id="S5",
        file="DatasetS5.xlsx",
        title="S5 — Proteins of unknown function on the IDR-ome map",
        sheets=[
            SheetSpec("Unknown function (neXtProt)", "Unknown_function", "list",
                      "1,512 UniProt IDs from the neXtProt 'Functional Proteome' set."),
            SheetSpec("Unknown function with IDRs", "Unknown_function_with_IDRs", "list"),
            SheetSpec("Clustered (unknown-function with IDRs)", "Clustered", "list"),
            SheetSpec("Clustered & GO-overrepresented", "With GO terms", "list"),
        ],
    ),
    "S6": DatasetSpec(
        dataset_id="S6",
        file="DatasetS6.xlsx",
        title="S6 — FAIDR GO-prediction model performance",
        sheets=[
            SheetSpec("A) FAIDR Testing", "A) FAIDR Testing", "flat",
                      "Per-GO-term test-set classification performance (601 models)."),
            SheetSpec("B) Test-set probabilities (top 146)", "B) Test Set Probabilities", "flat",
                      "Per-protein FAIDR probabilities for the 146 best models."),
            SheetSpec("C) IDR prediction statistics", "C) IDR Prediction Statistics", "flat"),
            SheetSpec("D) IDR assignment matrix", "D) IDR Assignment Matrix", "matrix",
                      "Per-IDR × GO-term assignment fractions across 100 FAIDR repetitions. "
                      "Use the directional queries below to find IDRs assigned to a given "
                      "GO term, or GO terms assigned to a given IDR, above a chosen threshold."),
            SheetSpec("E) Consensus-filtered assignment stats", "E) IDR Consistent Assign Stats", "flat"),
            SheetSpec("F) Representative GO terms", "F) Representative IDR Assign", "flat"),
        ],
    ),
    "S7": DatasetSpec(
        dataset_id="S7",
        file="DatasetS7.xlsx",
        title="S7 — Disease-risk genes & biomolecular condensates",
        sheets=[
            SheetSpec("A) Disease-risk in selected clusters", "A) Disease-risk select", "flat"),
            SheetSpec("B) Disease-risk in automatic clusters", "B) Disease-risk auto", "flat"),
            SheetSpec("C) Condensate proteins in selected clusters", "C) Biomol condens select", "flat"),
            SheetSpec("D) Condensate proteins in automatic clusters", "D) Biomol condens auto", "flat"),
            SheetSpec("E) FAIDR top 10% — biomolecular condensates", "E) FAIDR Biomol condens", "flat"),
            SheetSpec("F) FAIDR top 10% — ASD-risk", "F) FAIDR ASD", "flat"),
        ],
    ),
    "S8": DatasetSpec(
        dataset_id="S8",
        file="DatasetS8.xlsx",
        title="S8 — SLiM-based clustering & interaction analyses",
        sheets=[
            SheetSpec("A) SLiM clusters", "A) SLiM_CLUSTERS", "raw",
                      "SLiM-dominated clusters with their member proteins (nested blocks)."),
            SheetSpec("B) SLiM ↔ binding-domain interactors",
                      "B) SLiM_DOMAINS_INTERACTORS", "no_banner",
                      "Enrichment of canonical SLiM-binding domains among interactors."),
        ],
    ),
}


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def list_dataset_ids() -> List[str]:
    return list(REGISTRY.keys())


def get_spec(dataset_id: str) -> DatasetSpec:
    return REGISTRY[dataset_id]


def datasets_dir(zenodo_root: str) -> str:
    return os.path.join(zenodo_root, "DATASETS")


@st.cache_data(show_spinner=False)
def load_description(zenodo_root: str, dataset_id: str) -> str:
    spec = REGISTRY[dataset_id]
    fp = os.path.join(datasets_dir(zenodo_root), spec.file)
    try:
        xl = pd.ExcelFile(fp)
    except FileNotFoundError:
        return f"_File not found: {spec.file}_"
    desc_sheet = next(
        (s for s in xl.sheet_names if s.lower().startswith("desc")), None
    )
    if desc_sheet is None:
        return ""
    raw = pd.read_excel(xl, sheet_name=desc_sheet, header=None)
    parts: List[str] = []
    for _, row in raw.iterrows():
        for v in row.tolist():
            if isinstance(v, str) and v.strip():
                parts.append(v.strip())
                break
    return "\n\n".join(parts).strip()


@st.cache_data(show_spinner="Loading sheet…")
def load_sheet(
    zenodo_root: str, dataset_id: str, sheet_label: str
) -> Tuple[pd.DataFrame, str]:
    """Load and tidy one sheet. Returns (df, layout)."""
    spec = REGISTRY[dataset_id]
    fp = os.path.join(datasets_dir(zenodo_root), spec.file)
    sh = next((s for s in spec.sheets if s.label == sheet_label), None)
    if sh is None:
        raise KeyError(f"{dataset_id} has no sheet labelled {sheet_label!r}")
    raw = pd.read_excel(fp, sheet_name=sh.sheet, header=None)
    raw = raw.dropna(how="all").reset_index(drop=True)

    if sh.layout == "list":
        # Take the first non-empty cell of each row in column 0.
        col0 = raw.iloc[:, 0].dropna().astype(str).str.strip()
        # Drop legend strings ("n = 718", description sentences) by keeping
        # only ID-like tokens.
        ids = [v for v in col0 if re.match(r"^[A-Z0-9][A-Z0-9_]{2,}$", v)]
        return pd.DataFrame({"UniProt_ID": ids}), sh.layout

    if sh.layout == "no_banner":
        df = _df_with_header(raw, header_row=0)
        return df, sh.layout

    if sh.layout == "flat":
        # Drop leading banner rows (long sentences in any of the first cols)
        offset = 0
        for i in range(min(3, len(raw))):
            if _is_banner_row(raw.iloc[i]):
                offset = i + 1
            else:
                break
        df = _df_with_header(raw.iloc[offset:].reset_index(drop=True), header_row=0)
        return df, sh.layout

    # 'raw' — just drop a banner row if present and keep generic columns.
    if len(raw) > 0 and _is_banner_row(raw.iloc[0]):
        raw = raw.iloc[1:].reset_index(drop=True)
    raw.columns = [f"col_{i}" for i in range(raw.shape[1])]
    return raw, sh.layout


def _is_banner_row(row: pd.Series) -> bool:
    """A row is a banner if it has *very few* populated cells (so it isn't
    a real header row) and any of those cells holds a long sentence
    starting with "Tab " or is just long descriptive text."""
    populated: List[str] = []
    for v in row.tolist():
        if isinstance(v, str):
            s = v.strip()
            if s and s.lower() != "nan":
                populated.append(s)
        elif v is not None and not (isinstance(v, float) and np.isnan(v)):
            populated.append(str(v))
    if not populated or len(populated) > 2:
        # Real header rows always populate most/all columns.
        return False
    for s in populated:
        if s.startswith("Tab ") or s.startswith('"Tab') or len(s) > 70:
            return True
    return False


@st.cache_data(show_spinner="Loading FAIDR assignment matrix…")
def load_assignment_matrix(zenodo_root: str) -> pd.DataFrame:
    """Load Dataset S6 Tab D as a (n_idrs × n_go_terms) DataFrame with the
    IDR ID as the index and GO term identifiers as columns."""
    fp = os.path.join(datasets_dir(zenodo_root), "DatasetS6.xlsx")
    df = pd.read_excel(fp, sheet_name="D) IDR Assignment Matrix", header=1)
    df = df.dropna(how="all").reset_index(drop=True)
    idr_col = df.columns[0]
    df = df.rename(columns={idr_col: "IDR"})
    df = df.dropna(subset=["IDR"]).reset_index(drop=True)
    df["IDR"] = df["IDR"].astype(str).str.strip()
    df = df.set_index("IDR")
    # Coerce all data columns to floats; non-numeric becomes NaN.
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _df_with_header(raw: pd.DataFrame, header_row: int) -> pd.DataFrame:
    if raw.empty or header_row >= len(raw):
        return raw
    header = raw.iloc[header_row].tolist()
    seen: Dict[str, int] = {}
    cols: List[str] = []
    for i, v in enumerate(header):
        name = str(v) if pd.notna(v) and str(v).strip() else f"col_{i}"
        seen[name] = seen.get(name, 0) + 1
        if seen[name] > 1:
            name = f"{name}_{seen[name] - 1}"
        cols.append(name)
    df = raw.iloc[header_row + 1 :].copy()
    df.columns = cols
    df = df.reset_index(drop=True)
    # Drop trailing fully-empty columns (common in xlsx exports).
    while df.shape[1] > 1 and df.iloc[:, -1].isna().all():
        df = df.iloc[:, :-1]
    return df
