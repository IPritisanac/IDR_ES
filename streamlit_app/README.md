# IDR-ome Streamlit explorer

An interactive web app for exploring the human IDR-ome evolutionary
signature map and associated data (Pritišanac et al. PNAS 2026). The app reads the data directly from the published Zenodo data repository.

## 1. Get the data

The full dataset (IDR-ome map, Cluster archives, per-IDR alignments, Supplementary datasets) is hosted on Zenodo at [doi:10.5281/zenodo.10812874](https://doi.org/10.5281/zenodo.10812874). From the repository root:

```bash
python download_zenodo.py --target ./ZENODO
```

This will download ~350 MB across ten archives and unpack the `.zip` files into `./ZENODO/`. See the top-level [README](../README.md) for details and per-archive contents.

The app expects this layout:

```
<ZENODO>/
    MAP/HUMAN_ES.txt               # full IDR-ome Z-score map (CDT format)
    MAP/HUMAN_ES.gtr               # dendrogram from Cluster3.0
    CLUSTERS_AUTO/CLUSTERS_0p4.zip # lists of different clusters and associated IDRs, across clustering thresholds
    CLUSTERS_AUTO/AUTO_GO_FEATS.xlsx
    IDROME_SEQUENCES/UP000005640_9606_SPOTD_MIN_30AA.fasta # human IDR sequences in fasta forma
    IDR_ALN/<UniProt>_…_ALN_IDR_<start>_<end>.fa # alignments of orthologous IDR sequences extracted from the full-length protein alignment files
    DATASETS/DatasetS1.xlsx … S8.xlsx # Supplementary Datasets accompanying the paper (Pritišanac, et al. PNAS 2026)
```

## 2. Install dependencies

```bash
cd streamlit_app
pip install -r requirements.txt
```

Python 3.9+ is recommended.

## 3. Run

```bash
streamlit run app.py
```

Override the path to ZENODO directory either by typing it into the sidebar text field, or via an environment variable before launching:

```bash
export IDR_ES_ZENODO=/abs/path/to/ZENODO
streamlit run app.py
```

## What the app supports

- **Find a protein.** Search by gene name, UniProt ID, or `UNIPROT_IDR_<start>_<end>`. For a selected IDR the app shows the 144-feature Z-score profile, the cluster it belongs to at the current threshold, and its human IDR amino-acid sequence.
- **Browse clusters.** Pick a correlation threshold (0.4 – 0.8) and browse the clusters. For each cluster the app shows the mean feature profile, a Z-score heatmap across members × features, a member table with gene names and IDR coordinates, the overrepresented GO terms, and the significantly enriched (positively / negatively) sequence features.
- **Supplementary datasets.** Browse Datasets S1 – S8 from `<ZENODO>/DATASETS/`: feature lists (S1), exploratory and automatic cluster GO/feature blocks plus IDR↔cluster annotations (S2), IDR lists per-clusters featured in paper figures (S3), the IDP and unknown-function subsets of the IDR-ome map (S4–S5), FAIDR GO-classification performance and per-IDR/per-protein GO assignments (S6), disease-risk and biomolecular-condensate clusters with FAIDR top-10% predictions (S7), and SLiM-based clusters + interactors analyses (S8). Each tab supports substring search and CSV download.
