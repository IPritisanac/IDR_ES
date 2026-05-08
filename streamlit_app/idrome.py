"""Data loaders for the IDR-ome Streamlit explorer.

The app reads directly from a Zenodo-style data folder with the following
layout (same as the one shipped with the paper):

    <ZENODO>/
        MAP/HUMAN_ES.txt          # full IDR-ome Z-score map (CDT format)
        MAP/HUMAN_ES.gtr          # hierarchical tree from Cluster3.0
        CLUSTERS_AUTO/CLUSTERS_0p4.zip ... CLUSTERS_0p8.zip
        CLUSTERS_AUTO/AUTO_GO_FEATS.xlsx
        IDROME_SEQUENCES/UP000005640_9606_SPOTD_MIN_30AA.fasta
        IDR_ALN/<UniProt>_<ENSG>_ENSEMBL_ORTHOLOGUES_ALN_IDR_<start>_<end>.fa

Every loader that touches disk is wrapped in ``st.cache_data`` so the user
pays the parsing cost at most once per session/threshold.
"""

from __future__ import annotations

import os
import re
import zipfile
from dataclasses import dataclass
from glob import glob
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st


THRESHOLDS: Tuple[float, ...] = (0.4, 0.5, 0.6, 0.7, 0.8)


def _thr_label(thr: float) -> str:
    """0.5 -> '0p5'."""
    return f"0p{int(round(thr * 10))}"


def _thr_dot(thr: float) -> str:
    """0.5 -> '0.5' (used in cluster CDT filenames)."""
    return f"{thr:.1f}"


# ---------------------------------------------------------------------------
# Main IDR-ome map (MAP/HUMAN_ES.txt)
# ---------------------------------------------------------------------------

@dataclass
class IdromeMap:
    features: List[str]           # 148 feature column names
    idrids: List[str]             # row order in CDT (post-cluster)
    zscores: np.ndarray           # (n_idrs, n_feats)
    gid_to_idrid: Dict[str, str]  # GENE####X -> IDRID
    idrid_to_gid: Dict[str, str]
    uniprot_to_idrids: Dict[str, List[str]]


def _resolve_path(zenodo_root: str, *candidates: str) -> str:
    """Return the first existing path among ``zenodo_root/<candidate>``.
    Lets the app accept either the manually-curated ``MAP/`` layout or
    the verbatim Zenodo archive layout ``ES_MAP/``, etc."""
    for c in candidates:
        full = os.path.join(zenodo_root, *c.split("/"))
        if os.path.exists(full):
            return full
    return os.path.join(zenodo_root, *candidates[0].split("/"))


@st.cache_data(show_spinner="Loading IDR-ome map…")
def load_idrome_map(zenodo_root: str) -> IdromeMap:
    path = _resolve_path(
        zenodo_root,
        "MAP/HUMAN_ES.txt",
        "ES_MAP/HUMAN_ES.txt",
    )
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Cannot find HUMAN_ES.txt under {zenodo_root}/MAP or "
            f"{zenodo_root}/ES_MAP. If you downloaded from Zenodo, the "
            "ES_MAP.zip archive should have been unpacked here."
        )
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    # CDT columns: GID, IDRID, NAME, GWEIGHT, then features
    feats = header[4:]
    df = pd.read_csv(
        path, sep="\t", skiprows=2, header=None, names=header,
        na_values=["nan", "NA", ""], low_memory=False,
    )
    # Drop EWEIGHT or any stray rows without a real IDRID
    df = df[df["IDRID"].astype(str).str.contains("_IDR_", na=False)].reset_index(drop=True)
    gid = df["GID"].astype(str).tolist()
    idrid = df["IDRID"].astype(str).tolist()
    zs = df[feats].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    uniprot_to_idrids: Dict[str, List[str]] = {}
    for iid in idrid:
        parsed = parse_idrid(iid)
        if parsed is None:
            continue
        uniprot_to_idrids.setdefault(parsed["uniprot"], []).append(iid)

    return IdromeMap(
        features=feats,
        idrids=idrid,
        zscores=zs,
        gid_to_idrid=dict(zip(gid, idrid)),
        idrid_to_gid=dict(zip(idrid, gid)),
        uniprot_to_idrids=uniprot_to_idrids,
    )


def parse_idrid(idrid: str) -> Optional[Dict[str, object]]:
    """e.g. 'Q2M385_IDR_1_34' -> {'uniprot': 'Q2M385', 'start': 1, 'end': 34}."""
    m = re.match(r"^([A-Za-z0-9]+)_IDR_(\d+)_(\d+)$", idrid)
    if not m:
        return None
    return {"uniprot": m.group(1), "start": int(m.group(2)), "end": int(m.group(3))}


# ---------------------------------------------------------------------------
# Cluster assignments (CLUSTERS_AUTO/CLUSTERS_0p?.zip)
# ---------------------------------------------------------------------------

@dataclass
class ClusterIndex:
    threshold: float
    clusters: Dict[int, List[str]]       # cluster_id -> [member IDRIDs]
    member_to_cluster: Dict[str, int]    # member IDRID -> cluster_id


@st.cache_data(show_spinner="Indexing clusters…")
def load_cluster_index(zenodo_root: str, threshold: float) -> ClusterIndex:
    zip_path = os.path.join(
        zenodo_root, "CLUSTERS_AUTO", f"CLUSTERS_{_thr_label(threshold)}.zip"
    )
    if not os.path.isfile(zip_path):
        raise FileNotFoundError(f"Cannot find cluster archive: {zip_path}")
    clusters: Dict[int, List[str]] = {}
    member_to_cluster: Dict[str, int] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            name = info.filename
            m = re.search(r"CLUSTER_(\d+)_[0-9p.]+\.out\.cdt$", name)
            if not m:
                continue
            cid = int(m.group(1))
            with zf.open(info) as f:
                raw = f.read().decode("utf-8", errors="replace").splitlines()
            members: List[str] = []
            # header line 0, EWEIGHT line 1, data from line 2
            for ln in raw[2:]:
                cells = ln.split("\t")
                if len(cells) >= 2 and "_IDR_" in cells[1]:
                    members.append(cells[1])
            if members:
                clusters[cid] = members
                for mm in members:
                    member_to_cluster[mm] = cid
    return ClusterIndex(threshold=threshold,
                       clusters=clusters,
                       member_to_cluster=member_to_cluster)


@st.cache_data(show_spinner="Loading cluster…")
def load_cluster_cdt(
    zenodo_root: str, threshold: float, cluster_id: int
) -> pd.DataFrame:
    """Return a DataFrame of cluster members × features (Z-scores)."""
    zip_path = os.path.join(
        zenodo_root, "CLUSTERS_AUTO", f"CLUSTERS_{_thr_label(threshold)}.zip"
    )
    target = f"CLUSTERS_{_thr_label(threshold)}/CLUSTER_{cluster_id}_{_thr_dot(threshold)}.out.cdt"
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(target) as f:
            raw = f.read().decode("utf-8", errors="replace").splitlines()
    header = raw[0].split("\t")
    feats = header[4:]
    rows: List[List[float]] = []
    idrids: List[str] = []
    for ln in raw[2:]:
        cells = ln.split("\t")
        if len(cells) < 5 or "_IDR_" not in cells[1]:
            continue
        idrids.append(cells[1])
        vals: List[float] = []
        for v in cells[4:4 + len(feats)]:
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                vals.append(np.nan)
        while len(vals) < len(feats):
            vals.append(np.nan)
        rows.append(vals)
    return pd.DataFrame(np.asarray(rows, dtype=float), index=idrids, columns=feats)


# ---------------------------------------------------------------------------
# GO + feature enrichments (CLUSTERS_AUTO/AUTO_GO_FEATS.xlsx)
# ---------------------------------------------------------------------------

@dataclass
class EnrichmentData:
    stats: pd.DataFrame                           # per-threshold summary
    go: Dict[float, Dict[str, pd.DataFrame]]      # thr -> {cluster_name: GO table}
    feats: Dict[float, Dict[str, pd.DataFrame]]   # thr -> {cluster_name: features table}


_GO_COLS = ["GO_term", "GO_description", "GO_category", "fold_enrichment",
            "p_value", "n_idrome", "n_cluster"]
_FEAT_COLS = ["direction", "feature", "fold_enrichment", "p_value"]


@st.cache_data(show_spinner="Loading GO enrichments…")
def load_enrichments(zenodo_root: str) -> EnrichmentData:
    path = os.path.join(zenodo_root, "CLUSTERS_AUTO", "AUTO_GO_FEATS.xlsx")
    xl = pd.ExcelFile(path)
    stats = pd.read_excel(xl, "AUTO_CLUSTERS_STATS", header=1)
    # First col is 'THRESHOLD' even though the xlsx puts a 'Table 1' banner
    stats = stats.rename(columns={stats.columns[0]: "THRESHOLD"})
    stats = stats.dropna(subset=["THRESHOLD"]).reset_index(drop=True)

    go_all: Dict[float, Dict[str, pd.DataFrame]] = {}
    feats_all: Dict[float, Dict[str, pd.DataFrame]] = {}

    for thr in THRESHOLDS:
        sh = f"CLUSTER_DISTANCE_{_thr_label(thr)}"
        if sh not in xl.sheet_names:
            continue
        raw = pd.read_excel(xl, sh, header=None)

        per_go: Dict[str, List[Dict[str, object]]] = {}
        per_feat: Dict[str, List[Dict[str, object]]] = {}

        current_cluster: Optional[str] = None
        mode: Optional[str] = None  # 'go' | 'feat'

        for i in range(len(raw)):
            row = raw.iloc[i].tolist()
            c0 = row[0]
            s0 = "" if (isinstance(c0, float) and np.isnan(c0)) else str(c0)

            if s0.startswith("CLUSTER_"):
                current_cluster = s0
                per_go.setdefault(current_cluster, [])
                per_feat.setdefault(current_cluster, [])
                mode = "go"
                continue
            if s0.upper().startswith("SIGNIFICANT FEATURES"):
                mode = "feat"
                continue
            if current_cluster is None:
                continue

            if mode == "go" and "#GO:" in s0:
                # col0 is 'name#GO:id'; col1-6 are description, category, fold, p, ratios
                name_go = s0
                if "#" in name_go:
                    _, go_term = name_go.split("#", 1)
                else:
                    go_term = ""
                per_go[current_cluster].append({
                    "GO_term": go_term,
                    "GO_description": row[1] if len(row) > 1 else "",
                    "GO_category": row[2] if len(row) > 2 else "",
                    "fold_enrichment": _to_float(row[3] if len(row) > 3 else None),
                    "p_value": _to_float(row[4] if len(row) > 4 else None),
                    "n_idrome": row[5] if len(row) > 5 else "",
                    "n_cluster": row[6] if len(row) > 6 else "",
                })
            elif mode == "feat" and s0 and s0 not in {"(+)", "(-)"}:
                direction, feature = _split_feat(s0)
                if feature is None:
                    continue
                per_feat[current_cluster].append({
                    "direction": direction,
                    "feature": feature,
                    "fold_enrichment": _to_float(row[1] if len(row) > 1 else None),
                    "p_value": _to_float(row[2] if len(row) > 2 else None),
                })

        go_all[thr] = {
            k: pd.DataFrame(v, columns=_GO_COLS)
            for k, v in per_go.items() if v
        }
        feats_all[thr] = {
            k: pd.DataFrame(v, columns=_FEAT_COLS)
            for k, v in per_feat.items() if v
        }

    return EnrichmentData(stats=stats, go=go_all, feats=feats_all)


def _to_float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _split_feat(s: str) -> Tuple[str, Optional[str]]:
    """'(+) id=REP_G2 G{2,2}' -> ('+', 'id=REP_G2 G{2,2}')."""
    m = re.match(r"^\s*\(([+\-])\)\s*(.+?)\s*$", s)
    if not m:
        return ("", None)
    return (m.group(1), m.group(2))


def cluster_name(cluster_id: int, threshold: float) -> str:
    """Match the keys used in AUTO_GO_FEATS.xlsx, e.g. CLUSTER_3048_0.5."""
    return f"CLUSTER_{cluster_id}_{_thr_dot(threshold)}"


# ---------------------------------------------------------------------------
# IDR sequences (IDROME fasta + per-IDR alignment files)
# ---------------------------------------------------------------------------

@dataclass
class IdromeSequences:
    uniprot_to_gene: Dict[str, str]            # P04637 -> 'P53' (UniProt entry name stem)
    gene_to_uniprots: Dict[str, List[str]]     # 'P53' (upper) -> [P04637, ...]


@st.cache_data(show_spinner="Indexing IDR sequences…")
def load_idrome_index(zenodo_root: str) -> IdromeSequences:
    """Build a UniProt↔entry-name index from the IDROME fasta. This is
    a "gene-like" name (e.g. ``P53`` for TP53, ``PECR`` for PECR, ``NUD4B``
    for NUD4B) — not always the HGNC symbol. Per-IDR sequences are
    fetched on demand from the IDR_ALN folder."""
    fasta = os.path.join(
        zenodo_root, "IDROME_SEQUENCES", "UP000005640_9606_SPOTD_MIN_30AA.fasta"
    )
    uni_to_gene: Dict[str, str] = {}
    if os.path.isfile(fasta):
        with open(fasta, "r", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if not ln.startswith(">"):
                    continue
                # '>sp|Q9H8V3|ECT2_HUMAN.diso:1' or '>tr|A0A024R1R8|A0A024R1R8_HUMAN.diso:1'
                h = ln[1:].strip().split("|")
                if len(h) >= 3 and h[0] == "sp":
                    uni = h[1]
                    gene_tok = h[2].split("_HUMAN")[0]
                    if gene_tok and gene_tok != uni:
                        uni_to_gene.setdefault(uni, gene_tok)
    gene_to_unis: Dict[str, List[str]] = {}
    for u, g in uni_to_gene.items():
        gene_to_unis.setdefault(g.upper(), []).append(u)
    return IdromeSequences(uniprot_to_gene=uni_to_gene,
                           gene_to_uniprots=gene_to_unis)


@st.cache_data(show_spinner=False)
def load_idr_sequence(zenodo_root: str, idrid: str) -> Optional[str]:
    """Pull the HUMAN reference sequence for an IDR out of its alignment file."""
    parsed = parse_idrid(idrid)
    if parsed is None:
        return None
    pattern = os.path.join(
        zenodo_root, "IDR_ALN",
        f"{parsed['uniprot']}_*_ALN_IDR_{parsed['start']}_{parsed['end']}.fa",
    )
    matches = glob(pattern)
    if not matches:
        return None
    with open(matches[0], "r", encoding="utf-8", errors="replace") as f:
        seq_parts: List[str] = []
        in_human = False
        for ln in f:
            if ln.startswith(">"):
                if in_human:
                    break
                in_human = "HUMAN" in ln.upper()
                continue
            if in_human:
                seq_parts.append(ln.strip())
    seq = "".join(seq_parts).replace("-", "")
    return seq or None


# ---------------------------------------------------------------------------
# Small helpers used by the UI
# ---------------------------------------------------------------------------

def resolve_query(
    query: str,
    idrome: IdromeMap,
    idx: IdromeSequences,
) -> pd.DataFrame:
    """Return a small DataFrame (IDRID, UniProt, gene, start, end) of
    matches for a free-text query. Tries, in order: exact IDRID,
    UniProt accession, UniProt entry name (e.g. ``P53``), common
    gene-name suffix/prefix variants (``TP53``→``P53``), and finally
    substring fallbacks. Empty DataFrame if no match."""
    q = (query or "").strip()
    if not q:
        return pd.DataFrame(columns=["IDRID", "UniProt", "Gene", "start", "end"])
    q_up = q.upper()

    # 1. Exact IDRID match
    candidates: List[str] = []
    if q in idrome.idrid_to_gid:
        candidates.append(q)
    # 2. UniProt accession
    if not candidates and q_up in idrome.uniprot_to_idrids:
        candidates.extend(idrome.uniprot_to_idrids[q_up])
    if not candidates:
        for u, iids in idrome.uniprot_to_idrids.items():
            if u.upper() == q_up:
                candidates.extend(iids)
                break
    # 3. UniProt entry name -> UniProt(s)
    if not candidates and q_up in idx.gene_to_uniprots:
        for u in idx.gene_to_uniprots[q_up]:
            candidates.extend(idrome.uniprot_to_idrids.get(u, []))
    # 4. HGNC -> entry-name heuristics (TP53 <-> P53, EGFR<->EGFR, etc.)
    if not candidates:
        variants = {q_up}
        if q_up.startswith("TP") and q_up[2:].isalnum():
            variants.add("P" + q_up[2:])
        if q_up.startswith("P") and q_up[1:].isdigit():
            variants.add("TP" + q_up[1:])
        for v in variants:
            if v in idx.gene_to_uniprots:
                for u in idx.gene_to_uniprots[v]:
                    candidates.extend(idrome.uniprot_to_idrids.get(u, []))
                if candidates:
                    break
    # 5. Substring fallback: first in IDRIDs, then in entry-name keys
    if not candidates:
        candidates = [i for i in idrome.idrids if q_up in i.upper()][:50]
    if not candidates:
        entry_hits = [g for g in idx.gene_to_uniprots.keys() if q_up in g][:20]
        for g in entry_hits:
            for u in idx.gene_to_uniprots[g]:
                candidates.extend(idrome.uniprot_to_idrids.get(u, []))

    rows: List[Dict[str, object]] = []
    for iid in candidates:
        p = parse_idrid(iid)
        if p is None:
            continue
        rows.append({
            "IDRID": iid,
            "UniProt": p["uniprot"],
            "Gene": idx.uniprot_to_gene.get(p["uniprot"], ""),
            "start": p["start"],
            "end": p["end"],
        })
    return pd.DataFrame(rows, columns=["IDRID", "UniProt", "Gene", "start", "end"])


def idr_zscores(idrome: IdromeMap, idrid: str) -> Optional[pd.Series]:
    try:
        row = idrome.idrids.index(idrid)
    except ValueError:
        return None
    return pd.Series(idrome.zscores[row], index=idrome.features, name=idrid)


# ---------------------------------------------------------------------------
# Unified cluster collection (auto OR manually selected)
# ---------------------------------------------------------------------------

@dataclass
class ClusterCollection:
    """A grouping of IDR clusters that the Browse tab can render. Both the
    automatic clusters at a given correlation threshold and the manually
    selected clusters from Dataset S2 are exposed through this single
    structure so the UI doesn't need to branch."""

    name: str                                   # display name, e.g. "Automatic @ 0.7"
    kind: str                                   # 'auto' | 'selected'
    threshold: Optional[float]                  # only set for 'auto'
    clusters: Dict[str, List[str]]              # canonical_id -> [IDR ids]
    member_to_cluster: Dict[str, str]           # IDR id -> canonical_id
    go: Dict[str, pd.DataFrame]                 # canonical_id -> GO table
    feats: Dict[str, pd.DataFrame]              # canonical_id -> features table

    def display_id(self, canonical_id: str) -> str:
        """Strip the 'CLUSTER_' prefix and any trailing _<threshold>."""
        s = canonical_id
        if s.startswith("CLUSTER_"):
            s = s[len("CLUSTER_"):]
        if self.kind == "auto":
            # 'CLUSTER_1977_0.7' -> '1977'
            parts = s.rsplit("_", 1)
            if len(parts) == 2 and re.match(r"^[0-9]+\.[0-9]+$", parts[1]):
                s = parts[0]
        return s


@st.cache_data(show_spinner="Loading automatic clusters…")
def load_auto_collection(zenodo_root: str, threshold: float) -> ClusterCollection:
    ci = load_cluster_index(zenodo_root, threshold)
    enr = load_enrichments(zenodo_root)
    canonical = {
        cluster_name(cid, threshold): members
        for cid, members in ci.clusters.items()
    }
    m2c = {
        iid: cluster_name(cid, threshold)
        for iid, cid in ci.member_to_cluster.items()
    }
    return ClusterCollection(
        name=f"Automatic @ {threshold}",
        kind="auto",
        threshold=threshold,
        clusters=canonical,
        member_to_cluster=m2c,
        go=enr.go.get(threshold, {}),
        feats=enr.feats.get(threshold, {}),
    )


_DATASET_S2_PATH = ("DATASETS", "DatasetS2.xlsx")


def _normalize_cluster_id(s: str) -> str:
    """For matching Tab A names like 'CLUSTER_27A' against Tab B's '27_A'."""
    if s.startswith("CLUSTER_"):
        s = s[len("CLUSTER_"):]
    return s.replace("_", "").upper()


@st.cache_data(show_spinner="Loading selected clusters…")
def load_selected_collection(zenodo_root: str) -> ClusterCollection:
    """Parse Dataset S2 Tab A (per-cluster GO + feature blocks) and Tab B
    (IDR ↔ cluster annotations) into a ClusterCollection."""
    fp = os.path.join(zenodo_root, *_DATASET_S2_PATH)
    if not os.path.isfile(fp):
        raise FileNotFoundError(f"Cannot find {fp}")

    # Tab B → cluster membership
    b = pd.read_excel(fp, sheet_name="B) IDRs in selected clusters", header=1)
    b = b.dropna(subset=["IDR", "CLUSTER"])
    clusters: Dict[str, List[str]] = {}
    m2c: Dict[str, str] = {}
    for iid, cid in zip(b["IDR"].astype(str), b["CLUSTER"]):
        cid_str = str(cid).strip()
        if cid_str.endswith(".0"):  # ints round-tripped via float
            cid_str = cid_str[:-2]
        canonical = f"CLUSTER_{cid_str}"
        clusters.setdefault(canonical, []).append(str(iid).strip())
        m2c[str(iid).strip()] = canonical

    # Tab A → nested GO + features per cluster
    raw = pd.read_excel(fp, sheet_name="A) Selected clusters", header=None)
    go_per: Dict[str, List[dict]] = {}
    feats_per: Dict[str, List[dict]] = {}
    current: Optional[str] = None
    mode: Optional[str] = None
    for i in range(len(raw)):
        c0 = raw.iloc[i, 0]
        s0 = "" if pd.isna(c0) else str(c0).strip()
        if s0.startswith("CLUSTER_"):
            current = s0
            go_per.setdefault(current, [])
            feats_per.setdefault(current, [])
            mode = "go"
            continue
        if "SIGNIFICANT FEATURES" in s0.upper():
            mode = "feat"
            continue
        if current is None:
            continue
        if mode == "go" and "#GO:" in s0:
            go_per[current].append({
                "GO_term": s0.split("#", 1)[1] if "#" in s0 else "",
                "GO_description": raw.iloc[i, 1] if raw.shape[1] > 1 else "",
                "GO_category": raw.iloc[i, 2] if raw.shape[1] > 2 else "",
                "fold_enrichment": _to_float(raw.iloc[i, 3] if raw.shape[1] > 3 else None),
                "p_value": _to_float(raw.iloc[i, 4] if raw.shape[1] > 4 else None),
                "n_idrome": raw.iloc[i, 5] if raw.shape[1] > 5 else "",
                "n_cluster": raw.iloc[i, 6] if raw.shape[1] > 6 else "",
            })
        elif mode == "feat" and (s0.startswith("(+) ") or s0.startswith("(-) ")):
            d, f = _split_feat(s0)
            if f is not None:
                feats_per[current].append({
                    "direction": d,
                    "feature": f,
                    "fold_enrichment": _to_float(raw.iloc[i, 1] if raw.shape[1] > 1 else None),
                    "p_value": _to_float(raw.iloc[i, 2] if raw.shape[1] > 2 else None),
                })

    # Best-effort name matching between Tab A and Tab B
    by_norm = {_normalize_cluster_id(k): k for k in clusters.keys()}
    go: Dict[str, pd.DataFrame] = {}
    feats: Dict[str, pd.DataFrame] = {}
    for tab_a_name, rows in go_per.items():
        canonical = by_norm.get(_normalize_cluster_id(tab_a_name), tab_a_name)
        if rows:
            go[canonical] = pd.DataFrame(rows, columns=_GO_COLS)
    for tab_a_name, rows in feats_per.items():
        canonical = by_norm.get(_normalize_cluster_id(tab_a_name), tab_a_name)
        if rows:
            feats[canonical] = pd.DataFrame(rows, columns=_FEAT_COLS)

    return ClusterCollection(
        name="Selected (Dataset S2)",
        kind="selected",
        threshold=None,
        clusters=clusters,
        member_to_cluster=m2c,
        go=go,
        feats=feats,
    )


def member_zscores_df(idrome: IdromeMap, idrids: List[str]) -> pd.DataFrame:
    """Pull the Z-score matrix for an arbitrary list of IDR IDs from the
    main IDR-ome map. Used for selected clusters (no per-cluster CDT
    file exists) and as a unified fallback."""
    pos = {iid: i for i, iid in enumerate(idrome.idrids)}
    rows: List[np.ndarray] = []
    keep: List[str] = []
    for iid in idrids:
        i = pos.get(iid)
        if i is None:
            continue
        rows.append(idrome.zscores[i])
        keep.append(iid)
    if not rows:
        return pd.DataFrame(columns=idrome.features)
    return pd.DataFrame(np.vstack(rows), index=keep, columns=idrome.features)


@st.cache_data(show_spinner="Building autocomplete list…")
def build_search_options(
    _idrome_idrids: Tuple[str, ...],
    _uniprot_to_gene_items: Tuple[Tuple[str, str], ...],
) -> Tuple[List[str], Dict[str, str]]:
    """Build a sorted list of searchable strings for a Streamlit selectbox.

    Each option is either ``"GENE — UNIPROT"`` (when a gene name is known) or
    just ``"UNIPROT"``. Returned alongside a label→UniProt lookup dict.

    Hash-friendly args (tuples) so ``@st.cache_data`` can memoise.
    """
    uniprot_to_gene = dict(_uniprot_to_gene_items)
    # Collect the UniProts that actually appear in the IDR-ome map.
    uniprots: set = set()
    for iid in _idrome_idrids:
        m = re.match(r"^([A-Za-z0-9]+)_IDR_", iid)
        if m:
            uniprots.add(m.group(1))

    options: List[str] = []
    label_to_uniprot: Dict[str, str] = {}
    for u in sorted(uniprots):
        gene = uniprot_to_gene.get(u, "")
        if gene:
            label = f"{gene} — {u}"
        else:
            label = u
        options.append(label)
        label_to_uniprot[label] = u
    # Sort with gene-named entries first (alphabetical), then bare UniProts.
    named = sorted(o for o in options if " — " in o)
    bare = sorted(o for o in options if " — " not in o)
    options = named + bare
    return options, label_to_uniprot


def feature_short_name(feat: str) -> str:
    """Compact the long Cluster3.0 feature names for plot axes."""
    if feat.startswith("id=REP_"):
        tok = feat.split()[0]
        return tok.replace("id=", "")
    return feat
