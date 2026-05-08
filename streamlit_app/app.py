"""Interactive explorer for the human IDR-ome ES map.

Run with::

    streamlit run streamlit_app/app.py

Point the sidebar at the Zenodo data folder (default path below can be
edited or overridden via the ``IDR_ES_ZENODO`` environment variable).
"""

from __future__ import annotations

import os
from typing import List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from idrome import (
    ClusterCollection,
    ClusterIndex,
    EnrichmentData,
    IdromeMap,
    IdromeSequences,
    THRESHOLDS,
    build_search_options,
    cluster_name,
    feature_short_name,
    idr_zscores,
    load_auto_collection,
    load_cluster_cdt,
    load_cluster_index,
    load_enrichments,
    load_idr_sequence,
    load_idrome_index,
    load_idrome_map,
    load_selected_collection,
    member_zscores_df,
    parse_idrid,
    resolve_query,
)
import datasets as supp_datasets


# Default Zenodo data folder. Resolution order:
#   1. IDR_ES_ZENODO environment variable
#   2. <repo_root>/ZENODO  (the layout produced by download_zenodo.py)
#   3. <repo_root>/../ZENODO
# Whichever exists first wins; otherwise falls back to a placeholder the
# user must edit in the sidebar.
def _default_zenodo() -> str:
    env = os.environ.get("IDR_ES_ZENODO")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    for candidate in (
        os.path.join(here, "..", "ZENODO"),
        os.path.join(here, "ZENODO"),
        os.path.join(os.path.expanduser("~"), "ZENODO"),
    ):
        if os.path.isdir(candidate):
            return os.path.abspath(candidate)
    return os.path.abspath(os.path.join(here, "..", "ZENODO"))


DEFAULT_ZENODO = _default_zenodo()


# ---------------------------------------------------------------------------
# App shell
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="IDR-ome ES explorer",
    page_icon=":microscope:",
    layout="wide",
)

st.title("Human IDR-ome — evolutionary signature explorer")
st.caption(
    "Interactive companion to the PNAS article. Browse Z-score clusters of "
    "human intrinsically disordered regions (IDRs) or look up a protein of "
    "interest to see where it lands on the map."
)


# ---------------------------------------------------------------------------
# Sidebar — data source + global controls
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Data")
    zenodo_root = st.text_input(
        "Zenodo data folder",
        value=DEFAULT_ZENODO,
        help="Folder containing MAP/, CLUSTERS_AUTO/, IDROME_SEQUENCES/, IDR_ALN/.",
    )
    if not os.path.isdir(zenodo_root):
        st.error("That path does not exist on this machine.")
        st.stop()

    st.header("Cluster source")
    cluster_source = st.radio(
        "Use which cluster set?",
        options=["Automatic", "Selected (Dataset S2)"],
        index=0,
        help=(
            "Automatic: clusters derived by Cluster3.0 from the IDR-ome at "
            "the chosen correlation threshold (0.4 – 0.8). "
            "Selected: the manually curated clusters listed in Dataset S2 "
            "Tab A/B (highlighted in the paper figures)."
        ),
    )
    if cluster_source == "Automatic":
        threshold = st.select_slider(
            "Correlation cutoff",
            options=list(THRESHOLDS),
            value=0.7,
            help=(
                "Correlation-based linkage cutoff used to partition the dendrogram "
                "in Cluster3.0. Lower values give more, tighter clusters "
                "(e.g. 0.4 ≈ 7,800 clusters); higher values give fewer, coarser "
                "groupings (e.g. 0.7 ≈ 400 clusters, 0.8 ≈ 230)."
            ),
        )
    else:
        threshold = None


# ---------------------------------------------------------------------------
# Load data (cached)
# ---------------------------------------------------------------------------

try:
    idrome: IdromeMap = load_idrome_map(zenodo_root)
    seq_idx: IdromeSequences = load_idrome_index(zenodo_root)
    if cluster_source == "Automatic":
        coll: ClusterCollection = load_auto_collection(zenodo_root, threshold)
    else:
        coll = load_selected_collection(zenodo_root)
except FileNotFoundError as e:
    st.error(f"Missing data file: {e}")
    st.stop()

with st.sidebar:
    st.markdown("---")
    st.metric("IDRs in map", f"{len(idrome.idrids):,}")
    st.metric("Features", f"{len(idrome.features)}")
    st.metric(
        f"Clusters ({coll.name})",
        f"{len(coll.clusters):,}",
    )


# ---------------------------------------------------------------------------
# Cross-tab helpers
# ---------------------------------------------------------------------------

def _make_feature_df(z: pd.Series) -> pd.DataFrame:
    return pd.DataFrame({
        "feature": [feature_short_name(f) for f in z.index],
        "full_feature": z.index,
        "z": z.values,
    })


def _feature_bar(z: pd.Series, title: str) -> go.Figure:
    df = _make_feature_df(z)
    df = df.sort_values("z", ascending=False).reset_index(drop=True)
    fig = px.bar(
        df, x="feature", y="z",
        hover_data={"full_feature": True, "feature": False, "z": ":.2f"},
        color="z",
        color_continuous_scale="RdBu_r",
        range_color=[-max(abs(df["z"].min()), abs(df["z"].max()), 1e-6),
                     max(abs(df["z"].min()), abs(df["z"].max()), 1e-6)],
    )
    fig.update_layout(
        title=title,
        xaxis_title="",
        yaxis_title="Z-score",
        xaxis={"tickangle": -75, "tickfont": {"size": 9}},
        margin=dict(l=40, r=20, t=50, b=120),
        height=420,
        coloraxis_showscale=False,
    )
    return fig


def _cluster_heatmap(
    members: pd.DataFrame,
    gene_lookup: dict,
    max_rows: int = 120,
) -> go.Figure:
    # Order columns by |mean Z|, rows by their projection on that order
    mean_abs = members.abs().mean(axis=0).fillna(0.0)
    col_order = mean_abs.sort_values(ascending=False).index.tolist()
    M = members[col_order].fillna(0.0)
    # Row order by first principal direction (simple: by mean of top features)
    row_order = M.iloc[:, : min(10, M.shape[1])].mean(axis=1).sort_values(ascending=False).index.tolist()
    M = M.loc[row_order]

    truncated = len(M) > max_rows
    if truncated:
        # Keep the top and bottom rows by projection so both extremes are visible.
        top = M.iloc[: max_rows // 2]
        bot = M.iloc[-max_rows // 2 :]
        M = pd.concat([top, bot])

    row_labels = []
    for iid in M.index:
        p = parse_idrid(iid) or {}
        gene = gene_lookup.get(p.get("uniprot", ""), "")
        row_labels.append(f"{gene} · {iid}" if gene else iid)

    col_labels = [feature_short_name(c) for c in M.columns]
    # Fixed diverging palette: pale blues for negative Z (down to -15),
    # pale yellows for positive Z (up to +15), white near zero. Values
    # outside [-15, +15] saturate at the endpoint colors.
    Z_RANGE = 15.0
    colorscale = [
        [0.00, "#5b8db8"],  # Z = -15 : pale steel blue
        [0.25, "#a8c5dd"],  # Z = -7.5
        [0.45, "#e9eff6"],
        [0.50, "#ffffff"],  # Z =   0 : white
        [0.55, "#fff7d9"],
        [0.75, "#f5d878"],  # Z = +7.5
        [1.00, "#d4a017"],  # Z = +15 : pale goldenrod
    ]
    fig = go.Figure(data=go.Heatmap(
        z=M.values,
        x=col_labels,
        y=row_labels,
        colorscale=colorscale,
        zmin=-Z_RANGE, zmax=Z_RANGE,
        colorbar=dict(title="Z-score", tickvals=[-15, -10, -5, 0, 5, 10, 15]),
        hovertemplate="<b>%{y}</b><br>%{x}<br>Z = %{z:.2f}<extra></extra>",
    ))

    # Make the chart explicitly wide: ~30 px per feature column so labels
    # are readable. The caller embeds the figure in a horizontally
    # scrollable container.
    px_per_col = 30
    label_pad = 240          # left margin for IDR labels
    cbar_pad = 110           # right margin reserved for the colorbar
    fig_width = label_pad + px_per_col * max(M.shape[1], 1) + cbar_pad
    fig_height = max(320, 20 * len(row_labels) + 160)

    fig.update_layout(
        width=fig_width,
        height=fig_height,
        xaxis=dict(tickangle=-75, tickfont=dict(size=10), automargin=True),
        yaxis=dict(tickfont=dict(size=10), autorange="reversed", automargin=True),
        margin=dict(l=label_pad, r=cbar_pad, t=30, b=140),
    )
    if truncated:
        fig.add_annotation(
            xref="paper", yref="paper", x=0.0, y=1.04, showarrow=False,
            text=f"Showing top+bottom {max_rows} of {len(members)} members — use the members table for the full list.",
            font=dict(size=11, color="gray"),
        )
    return fig


def _render_assignment_matrix(zenodo_root: str, ds_id: str, sh) -> None:
    """Interactive query UI for Dataset S6 Tab D — the FAIDR per-IDR × GO-term
    assignment matrix. Supports both directional queries (pick a GO term →
    list IDRs above threshold; pick an IDR → list GO terms above threshold)."""
    try:
        mat = supp_datasets.load_assignment_matrix(zenodo_root)
    except FileNotFoundError as e:
        st.error(f"{e}")
        return
    except Exception as e:
        st.error(f"Could not load assignment matrix: {e}")
        return

    n_idrs, n_go = mat.shape
    info_a, info_b, info_c = st.columns(3)
    info_a.metric("IDRs", f"{n_idrs:,}")
    info_b.metric("GO terms", f"{n_go:,}")
    info_c.metric("Cells ≥ 0.5",
                  f"{int((mat.values >= 0.5).sum()):,}")

    direction = st.radio(
        "Query",
        options=[
            "By GO term  →  list IDRs",
            "By IDR  →  list GO terms",
        ],
        index=0,
        horizontal=True,
        key=f"matrix_dir_{ds_id}",
    )

    threshold_q = st.slider(
        "Minimum assignment fraction",
        min_value=0.0, max_value=1.0, value=0.5, step=0.05,
        help=(
            "An IDR/GO pair passes if it was assigned in at least this "
            "fraction of the 100 FAIDR repetitions. 0.5 / 0.6 are the "
            "thresholds used in the paper."
        ),
        key=f"matrix_thr_{ds_id}",
    )

    if direction.startswith("By GO term"):
        go_choice = st.selectbox(
            "GO term",
            options=list(mat.columns),
            key=f"matrix_go_{ds_id}",
        )
        col = mat[go_choice]
        hits = col[col >= threshold_q].sort_values(ascending=False)
        rows = []
        for iid, score in hits.items():
            p = parse_idrid(str(iid)) or {}
            rows.append({
                "IDR": iid,
                "UniProt": p.get("uniprot", ""),
                "Gene": seq_idx.uniprot_to_gene.get(p.get("uniprot", ""), ""),
                "start": p.get("start", np.nan),
                "end": p.get("end", np.nan),
                "score": float(score),
            })
        result = pd.DataFrame(rows)
        st.caption(
            f"**{len(result):,}** IDRs assigned to *{go_choice}* "
            f"in ≥ {threshold_q:.2f} of 100 FAIDR runs."
        )
        st.dataframe(result, use_container_width=True, hide_index=True,
                     height=520)
        st.download_button(
            "Download CSV",
            data=result.to_csv(index=False).encode("utf-8"),
            file_name=f"S6_TabD_GO_{go_choice}_thr{threshold_q:.2f}.csv",
            mime="text/csv",
            key=f"matrix_dl_go_{ds_id}",
        )
    else:
        # Build IDR options once, with gene/UniProt context for autocomplete.
        idr_options, idr_label_to_id = _matrix_idr_options(
            tuple(mat.index.tolist()),
            tuple(seq_idx.uniprot_to_gene.items()),
        )
        idr_label = st.selectbox(
            "IDR",
            options=idr_options,
            index=None,
            placeholder="Start typing a gene, UniProt, or IDR ID…",
            key=f"matrix_idr_{ds_id}",
        )
        if not idr_label:
            st.info("Pick an IDR to list its GO-term assignments.")
            return
        idr_id = idr_label_to_id[idr_label]
        row = mat.loc[idr_id]
        hits = row[row >= threshold_q].sort_values(ascending=False)
        result = pd.DataFrame({
            "GO_term": hits.index,
            "score": hits.values.astype(float),
        })
        st.caption(
            f"**{len(result):,}** GO terms assigned to *{idr_id}* "
            f"in ≥ {threshold_q:.2f} of 100 FAIDR runs."
        )
        st.dataframe(result, use_container_width=True, hide_index=True,
                     height=520)
        st.download_button(
            "Download CSV",
            data=result.to_csv(index=False).encode("utf-8"),
            file_name=f"S6_TabD_IDR_{idr_id}_thr{threshold_q:.2f}.csv",
            mime="text/csv",
            key=f"matrix_dl_idr_{ds_id}",
        )


@st.cache_data(show_spinner="Building IDR autocomplete…")
def _matrix_idr_options(
    _idr_ids: tuple,
    _uniprot_to_gene_items: tuple,
) -> tuple:
    uniprot_to_gene = dict(_uniprot_to_gene_items)
    options: list = []
    label_to_id: dict = {}
    for iid in _idr_ids:
        p = parse_idrid(str(iid)) or {}
        uni = p.get("uniprot", "")
        gene = uniprot_to_gene.get(uni, "")
        if gene:
            label = f"{gene} · {iid}"
        else:
            label = str(iid)
        options.append(label)
        label_to_id[label] = iid
    options.sort()
    return options, label_to_id


def _show_scrollable_plotly(fig: go.Figure) -> None:
    """Render a Plotly figure in a horizontally-scrollable container so the
    chart can be wider than the Streamlit column."""
    width = int(fig.layout.width or 1200)
    height = int(fig.layout.height or 500)
    inner_html = fig.to_html(include_plotlyjs="cdn", full_html=False,
                             config={"displaylogo": False})
    wrapper = f"""
    <div style="overflow-x: auto; overflow-y: hidden; width: 100%;
                border: 1px solid #eee; border-radius: 4px;">
      <div style="width: {width}px;">{inner_html}</div>
    </div>
    """
    components.html(wrapper, height=height + 40, scrolling=False)


# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

PAGES = ("Find a protein", "Browse clusters", "Supplementary datasets")
if "page" not in st.session_state:
    st.session_state["page"] = PAGES[0]
st.radio(
    "View",
    PAGES,
    key="page",
    horizontal=True,
    label_visibility="collapsed",
)
page = st.session_state["page"]


def _jump_to_cluster(canonical_cid: str) -> None:
    """Button callback: switch to the Browse page with the cluster preselected.
    Runs BEFORE the next script rerun, so the radio picks up the new page value
    cleanly without Streamlit complaining about widget-state writes."""
    st.session_state["page"] = "Browse clusters"
    st.session_state["browse_cluster"] = canonical_cid


# ------------------------------ Tab 1 -------------------------------------

if page == "Find a protein":
    st.subheader("Look up an IDR")

    search_options, label_to_uniprot = build_search_options(
        tuple(idrome.idrids),
        tuple(seq_idx.uniprot_to_gene.items()),
    )

    pick = st.selectbox(
        "Gene name or UniProt ID",
        options=search_options,
        index=None,
        placeholder="Start typing a gene (e.g. TP53) or UniProt (e.g. Q9H8V3)…",
        key="find_pick",
    )

    with st.expander(
        "Or paste a free-form query (IDR ID, substring, …)", expanded=False
    ):
        free_q = st.text_input(
            "Free-form search",
            placeholder="e.g. Q9H8V3_IDR_832_914 or 'P04637' or 'kinase'",
            key="find_query_free",
            label_visibility="collapsed",
        )

    if pick:
        q = label_to_uniprot.get(pick, pick)
    else:
        q = free_q
    matches = resolve_query(q, idrome, seq_idx)
    if q and matches.empty:
        st.warning("No matches in the IDR-ome map.")
    elif not matches.empty:
        st.caption(f"{len(matches)} match(es)")
        st.dataframe(matches, use_container_width=True, hide_index=True)

        choice = st.selectbox(
            "Select an IDR to inspect",
            options=matches["IDRID"].tolist(),
            key="find_selected",
        )
        if choice:
            parsed = parse_idrid(choice) or {}
            uniprot = parsed.get("uniprot", "")
            gene = seq_idx.uniprot_to_gene.get(uniprot, "")
            cid = coll.member_to_cluster.get(choice)
            cid_display = coll.display_id(cid) if cid else None

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("UniProt", uniprot or "—")
            c2.metric("Gene", gene or "—")
            c3.metric(
                "IDR coords",
                f"{parsed.get('start', '?')}–{parsed.get('end', '?')}",
            )
            c4.metric(
                f"Cluster ({coll.name})",
                cid_display if cid_display else "—",
            )

            z = idr_zscores(idrome, choice)
            if z is not None:
                st.plotly_chart(
                    _feature_bar(z, f"Feature Z-score profile of {choice}"),
                    use_container_width=True,
                )

            seq = load_idr_sequence(zenodo_root, choice)
            if seq:
                with st.expander("IDR amino-acid sequence (HUMAN)", expanded=False):
                    st.code("\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)))
            else:
                st.caption(
                    "IDR sequence not found in IDR_ALN folder — check that "
                    "the alignment files are present."
                )

            if cid is not None:
                st.button(
                    "Jump to this cluster →",
                    on_click=_jump_to_cluster,
                    args=(cid,),
                    key=f"jump_{choice}",
                )
            if cid is None and coll.kind == "selected":
                st.caption(
                    "This IDR is not part of any of the manually selected "
                    "clusters in Dataset S2."
                )


# ------------------------------ Tab 2 -------------------------------------

elif page == "Browse clusters":
    st.subheader(f"Browsing clusters — {coll.name}")

    # Build a quick summary table (cluster id, size, top GO term if any)
    summary_rows: List[dict] = []
    for cid, members in coll.clusters.items():
        gdf = coll.go.get(cid)
        top_go = ""
        if gdf is not None and len(gdf):
            top_go = str(gdf.iloc[0]["GO_description"])
        summary_rows.append({
            "cluster_id": coll.display_id(cid),
            "_canonical": cid,
            "size": len(members),
            "unique_proteins": len({
                (parse_idrid(m) or {}).get("uniprot", m) for m in members
            }),
            "top_GO_term": top_go,
        })
    summary = pd.DataFrame(summary_rows).sort_values(
        ["size", "cluster_id"], ascending=[False, True]
    ).reset_index(drop=True)

    left, right = st.columns([1, 2])
    with left:
        min_size = st.number_input(
            "Minimum size",
            min_value=1,
            max_value=max(summary["size"].max() if len(summary) else 1, 1),
            value=1,
        )
        only_with_go = st.checkbox("Only clusters with overrepresented GO term",
                                   value=False)
        view = summary[summary["size"] >= min_size]
        if only_with_go:
            view = view[view["top_GO_term"].astype(bool)]
        st.caption(f"{len(view):,} clusters match")
        st.dataframe(
            view.drop(columns=["_canonical"]),
            use_container_width=True, hide_index=True, height=420,
        )

    preferred = st.session_state.pop("browse_cluster", None)
    canonical_ids_for_picker = view["_canonical"].tolist() if len(view) else list(
        coll.clusters.keys()
    )
    if preferred is not None:
        if preferred not in canonical_ids_for_picker:
            # The preferred cluster fell outside the current size/GO filters;
            # surface it anyway by prepending so the user actually lands on it.
            canonical_ids_for_picker = [preferred] + canonical_ids_for_picker
            st.info(
                f"Showing cluster {coll.display_id(preferred)} from your "
                "selection on the Find page (it was outside the current "
                "size/GO filters)."
            )
        # Force-override the selectbox's session state so the jump takes
        # effect even after the user previously interacted with the picker.
        st.session_state["cluster_picker"] = preferred

    default_index = 0
    current_pick = st.session_state.get("cluster_picker")
    if current_pick in canonical_ids_for_picker:
        default_index = canonical_ids_for_picker.index(current_pick)

    with right:
        if not canonical_ids_for_picker:
            st.info("No clusters match the current filters.")
            st.stop()
        chosen_cid = st.selectbox(
            "Cluster",
            options=canonical_ids_for_picker,
            index=default_index,
            format_func=lambda c: (
                f"cluster {coll.display_id(c)}  ·  n={len(coll.clusters.get(c, []))}"
            ),
            key="cluster_picker",
        )

        # Members + Z-score matrix. For automatic clusters we keep the
        # exact CDT row order (matches JavaTreeView); for selected
        # clusters we pull rows from the main IDR-ome map.
        if coll.kind == "auto" and coll.threshold is not None:
            # display_id strips the threshold suffix, parse it back to int
            try:
                cid_int = int(coll.display_id(chosen_cid))
                members_df = load_cluster_cdt(zenodo_root, coll.threshold, cid_int)
            except (ValueError, FileNotFoundError):
                members_df = member_zscores_df(idrome, coll.clusters[chosen_cid])
        else:
            members_df = member_zscores_df(idrome, coll.clusters[chosen_cid])
        cname = chosen_cid

        c1, c2, c3 = st.columns(3)
        c1.metric("Members", f"{len(members_df):,}")
        c2.metric("Unique proteins", f"{len({(parse_idrid(m) or {}).get('uniprot', m) for m in members_df.index}):,}")
        c3.metric("Threshold", f"{threshold}")

        # Member table with gene names and coords
        member_rows = []
        for iid in members_df.index:
            p = parse_idrid(iid) or {}
            member_rows.append({
                "IDRID": iid,
                "UniProt": p.get("uniprot", ""),
                "Gene": seq_idx.uniprot_to_gene.get(p.get("uniprot", ""), ""),
                "start": p.get("start", np.nan),
                "end": p.get("end", np.nan),
                "length": (
                    p.get("end", 0) - p.get("start", 0) + 1
                    if p.get("start") is not None and p.get("end") is not None else np.nan
                ),
            })
        mem_table = pd.DataFrame(member_rows)

        st.markdown("#### Mean feature profile")
        mean_z = members_df.mean(axis=0, skipna=True)
        st.plotly_chart(
            _feature_bar(mean_z, f"Mean Z-score across {len(members_df)} members"),
            use_container_width=True,
        )

        st.markdown("#### Feature heatmap")
        st.caption(
            "Scroll horizontally to see all 144 features."
        )
        _show_scrollable_plotly(
            _cluster_heatmap(members_df, seq_idx.uniprot_to_gene)
        )

        st.markdown("#### Members")
        st.dataframe(mem_table, use_container_width=True, hide_index=True)

        # Optional: inspect a member's sequence
        with st.expander("Inspect member sequence", expanded=False):
            m_sel = st.selectbox(
                "IDR",
                options=mem_table["IDRID"].tolist(),
                key=f"member_seq_{chosen_cid}",
            )
            if m_sel:
                seq = load_idr_sequence(zenodo_root, m_sel)
                if seq:
                    st.code("\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)))
                else:
                    st.caption("Sequence not found in IDR_ALN folder.")

        # GO + feature enrichments
        g_col, f_col = st.columns(2)
        with g_col:
            st.markdown("#### Overrepresented GO terms")
            gdf = coll.go.get(cname)
            if gdf is None or gdf.empty:
                st.caption("No significant GO terms recorded for this cluster.")
            else:
                st.dataframe(
                    gdf.sort_values("p_value"),
                    use_container_width=True,
                    hide_index=True,
                )
        with f_col:
            st.markdown("#### Significantly enriched features")
            fdf = coll.feats.get(cname)
            if fdf is None or fdf.empty:
                st.caption("No significant feature enrichments recorded.")
            else:
                st.dataframe(
                    fdf.sort_values(["direction", "p_value"],
                                    ascending=[False, True]),
                    use_container_width=True,
                    hide_index=True,
                )


# ------------------------------ Tab 3 -------------------------------------

elif page == "Supplementary datasets":
    st.subheader("Supplementary datasets")
    st.caption(
        "Datasets S1 – S8 from `<ZENODO>/DATASETS/`. Pick a dataset and a "
        "tab to inspect the table; use the download button for a CSV copy."
    )

    if not os.path.isdir(supp_datasets.datasets_dir(zenodo_root)):
        st.error(
            "No DATASETS folder found at "
            f"`{supp_datasets.datasets_dir(zenodo_root)}`. Add it to the "
            "Zenodo data folder to populate this tab."
        )
    else:
        ds_ids = supp_datasets.list_dataset_ids()
        ds_titles = {
            d: supp_datasets.get_spec(d).title for d in ds_ids
        }
        ds_id = st.selectbox(
            "Dataset",
            options=ds_ids,
            format_func=lambda d: ds_titles[d],
            key="dataset_picker",
        )
        spec = supp_datasets.get_spec(ds_id)

        with st.expander("About this dataset", expanded=False):
            desc = supp_datasets.load_description(zenodo_root, ds_id)
            st.markdown(desc or "_No description._")

        sheet_labels = [s.label for s in spec.sheets]
        sheet_label = st.selectbox(
            "Tab", options=sheet_labels, key=f"sheet_picker_{ds_id}"
        )
        sh = next(s for s in spec.sheets if s.label == sheet_label)
        if sh.note:
            st.caption(sh.note)

        if sh.layout == "matrix":
            _render_assignment_matrix(zenodo_root, ds_id, sh)
        else:
            try:
                df_sheet, layout = supp_datasets.load_sheet(
                    zenodo_root, ds_id, sheet_label
                )
            except Exception as e:
                st.error(f"Could not load sheet: {e}")
            else:
                ds_top_a, ds_top_b = st.columns([3, 1])
                ds_top_a.metric("Rows", f"{len(df_sheet):,}")
                ds_top_b.metric("Columns", f"{df_sheet.shape[1]:,}")

                search = st.text_input(
                    "Search (substring across all columns)",
                    key=f"dataset_search_{ds_id}_{sheet_label}",
                )
                if search:
                    mask = df_sheet.astype(str).apply(
                        lambda c: c.str.contains(search, case=False, na=False)
                    ).any(axis=1)
                    df_view = df_sheet[mask].reset_index(drop=True)
                    st.caption(f"{len(df_view):,} of {len(df_sheet):,} rows match.")
                else:
                    df_view = df_sheet

                st.dataframe(df_view, use_container_width=True, hide_index=True,
                             height=520)
                st.download_button(
                    "Download CSV",
                    data=df_view.to_csv(index=False).encode("utf-8"),
                    file_name=f"Dataset{ds_id}_{sh.sheet}.csv",
                    mime="text/csv",
                    key=f"dataset_dl_{ds_id}_{sheet_label}",
                )
