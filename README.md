# IDR_ES

A Python package for calculation of **Z-score signatures** of intrinsically disordered regions (IDRs) of proteins. Two complementary protocols are supported by the same codebase:

- **ES — Evolutionary Signatures.** For each molecular feature, a mean Z-score is computed. Z-score evaluates how far the mean feature value measured across a set of orthologous IDR sequences deviates from a null hypothesis built from pairwise simulations. Requires an MSA per IDR.
- **FS — Feature Signatures.** Sequence-only Z-scores: each feature value of an individual IDR is normalized by the global mean and standard deviation of that feature across a whole IDRome (one FASTA file with all IDRs of interest). Does *not* require orthologous sequences.

Both protocols are implemented in `src/core/es_pw_sim.py` (class `EVSign`) and are dispatched from the top-level script `run_es.py`.

## Main features

### Shared across ES and FS

- Parsing of a plain-text input file (`src/utils/input_file.txt`) that configures parameter files, input paths, and protocol variables.
- A fixed set of knowledge-based molecular features computed from the IDR sequence (amino-acid composition, short linear motif counts, homorepeat content, physicochemical features). The feature list is defined in `src/core/sequence_features.py` and can be extended by the user.
- Per-feature Z-score output as a tab-separated file, written to an `output/` directory that is auto-created when `run_es.py` is run.

### ES protocol (evolutionary signatures)

- Reads one alignment file (FASTA, `*.fa`/`*.fasta`) per IDR from a user-specified directory (`align_dir`).
- Selects a reference sequence inside each alignment either by index (`REF_NUM`) or by name/substring match (`REF_NAME`, e.g. `HUMAN`).
- Applies sequence quality-control heuristics (`L_MIN`, `L_FACTOR`, `D_RATIO`, `D_TOTAL`) to pick which orthologous sequences to keep.
- Computes a proxy pairwise evolutionary distance (F81 model) between the reference IDR (human) sequence and each orthologuous sequence — **no phylogenetic tree is used**.
- Builds a null distribution via `n_simulations` pairwise simulations per orthologue under JC69-style amino-acid substitution, with a simple indel model (optional, set with params - `use_indels`, `max_indel_size`, `indel_prob`).
- Returns per-feature **mean Z-scores** (how far the observed mean across orthologues deviates from the simulated null). Per-feature variance Z-scores are currently disabled in this version of the code.
- Requires each kept alignment to contain a reference IDR of length ≥ 30 aa and ≥ 10 orthologous sequences; alignments that do not meet these thresholds are skipped with a message.

### FS protocol (feature signatures, sequence-only)

- Reads **a single FASTA file** (the "IDRome") from a user-specified directory (`fasta_dir`). If multiple FASTA files are present in the directory, the first one in alphabetical order is used and a warning is printed.
- De-gaps sequences, drops empty sequences and sequences containing `X`.
- Computes the raw value of every feature for every IDR (`_compute_features_for_seq`), then normalizes each feature by its global mean and standard deviation across all IDRs in the FASTA file (standard deviation is floored at `MIN_SD`).
- Returns per-sequence **mean Z-scores** for all features.

## Companion data (Zenodo)

This repository contains the code for ES/FS computation. 
The full human IDR-ome dataset that the paper analyses — the global (ES) Z-score map, the cluster(s) archives, per-IDR alignments, and supplementary datasets are available on Zenodo:

> **Pritisanac, I.** *Data repository associated with 'A Functional Map of the Human Intrinsically Disordered Proteome'.* Zenodo. [https://doi.org/10.5281/zenodo.10812874](https://doi.org/10.5281/zenodo.10812874)

The DOI above is the *concept* DOI, which always resolves to the latest published version of the deposit. The current version (v4) bundles ten archives totalling ~350 MB:

| Archive | Size | What's in it |
|---|---|---|
| `ES_MAP.zip` | 23 MB | clustered IDR-ome map (`HUMAN_ES.txt`, `HUMAN_ES.gtr`, `HUMAN_ES.cdt`) |
| `CLUSTERS_AUTO.zip` | 76 MB | per-threshold (0.4–0.8) cluster archives + `AUTO_GO_FEATS.xlsx` |
| `CLUSTERS_EXPLORE.zip` | 10 MB | 93 manually exported / exploratory clusters |
| `DATASETS.zip` | 3 MB | supplementary datasets S1–S8 |
| `IDR_ALN.zip` | 200 MB | 19,459 ortholog alignments per IDR |
| `IDROME_SEQUENCES.zip` | 10 MB | proteome + IDRome FASTA, SPOT-Disorder boundaries |
| `FAIDR_TSTATS.zip` | <1 MB | FAIDR t-stat hierarchical clustering |
| `FAIDR_HIGH_AUC_PPV_GO.zip` | 22 MB | 148 high-quality FAIDR target files |
| `PROTEIN_GROUPS_FAIDR_TARGETS.zip` | 2 MB | FAIDR target groups |
| `TUTORIAL.zip` | <1 MB | Cluster3.0 / JavaTreeView tutorial PDF |

To fetch everything into a local folder you can point the Streamlit app at, run from the repo root:

    python download_zenodo.py --target ./ZENODO

The script calls Zenodo's public API, downloads each archive in the deposit, and unpacks the `.zip` files into `./ZENODO/`. Optional flags:

- `--no-extract` — keep the archives without unpacking.
- `--only ES_MAP DATASETS CLUSTERS_AUTO IDROME` — only download files whose names contain those substrings (if you want to skip the 200 MB `IDR_ALN.zip`).
- `--record <id>` — pin to a specific version, e.g. `--record 20076448` for v3, instead of the default concept ID `10812874`.

After the download you should have the layout the Streamlit app and the analysis scripts expect:

    ZENODO/
        MAP/HUMAN_ES.txt              # full IDR-ome Z-score map (CDT format)
        MAP/HUMAN_ES.gtr              # dendrogram from Cluster3.0
        CLUSTERS_AUTO/CLUSTERS_0p4.zip … 0p8.zip
        CLUSTERS_AUTO/AUTO_GO_FEATS.xlsx
        IDROME_SEQUENCES/UP000005640_9606_SPOTD_MIN_30AA.fasta
        IDR_ALN/<UniProt>_…_ALN_IDR_<start>_<end>.fa
        DATASETS/DatasetS1.xlsx … S8.xlsx

Set `IDR_ES_ZENODO=/abs/path/to/ZENODO` (or paste the path into the sidebar of the Streamlit app) and you're ready to go.

## Streamlit explorer

An interactive web app, [`streamlit_app/`](streamlit_app/README.md), lets you:

- look up an IDR by gene name / UniProt / IDR ID and see its 144-feature Z-score profile,
- browse the automatic clusters at any of the 0.4–0.8 correlation thresholds, _or_ the manually selected clusters from Dataset S2 (PNAS 2026, Supplementary Datasets),
- inspect the Supplementary datasets (S1–S8), including a directional-query viewer for the FAIDR per-IDR×GO-term assignment matrix.

After downloading the data:

    cd streamlit_app
    pip install -r requirements.txt
    streamlit run app.py

See [`streamlit_app/README.md`](streamlit_app/README.md) for details.

## Installation / dependencies

The code can be run on Linux / macOS / Windows with Python 3 and the `numpy` and `scipy` packages installed.

1. Install Python 3 from https://www.python.org/downloads/ (if not already installed).
2. Make sure `pip` is available:

         pip --version

3. Install NumPy and SciPy:

         pip install numpy scipy

Alternatively, with Anaconda:

1. Install Anaconda / Miniconda: https://www.anaconda.com/products/distribution
2. Create and activate an environment (optional but recommended):

         conda create --name idr_es python=3
         conda activate idr_es

3. Install dependencies:

         conda install numpy scipy

## Running

`run_es.py` takes **two arguments**: the input file, and the protocol to execute (`ES` or `FS`). It validates the protocol name, instantiates `EVSign`, checks that the corresponding path variable is set in the input file, and then calls either `compute_es_dir()` or `compute_fs_seq()`.

    python run_es.py <input_file> <ES|FS>

Examples:

    # Evolutionary Signatures from a directory of IDR alignments
    python run_es.py src/utils/input_file.txt ES

    # Feature Signatures from a FASTA file of IDRs
    python run_es.py src/utils/input_file.txt FS

If the wrong protocol is requested for the configured input, `run_es.py` will print an explanatory error (e.g. `protocol 'FS' requires 'fasta_dir' to be set in …`) and exit.

## Input file

A single plain-text input file configures both protocols. An example is provided in `src/utils/input_file.txt`. The most relevant entries are:

| Key | Meaning |
|---|---|
| `motifs_file` | Path to the motif definition table (`MOTIFS.txt`). |
| `exp_motifs_n_file` | Path to the pre-computed expected-motif-count table. |
| `repeats_file` | Path to the homorepeat definition table. |
| `aa_freq_file` | Path to the background amino-acid composition table. |
| `align_dir` | **ES protocol only.** Directory with one alignment file per IDR (`*.fa` / `*.fasta`). |
| `fasta_dir` | **FS protocol only.** Directory containing a single FASTA file with all IDR sequences. |
| `use_indels` | `on`/`off` — toggle the simple indel model in the ES simulations. |
| `n_simulations` | Number of pairwise simulations per orthologous sequence (ES only). |
| `REF_NUM` | Default reference-sequence index in each alignment (ES only). |
| `REF_NAME` | Substring used to auto-pick the reference sequence by name (ES only). |
| `MIN_SD` | Floor on the per-feature standard deviation used in Z-score denominators. |
| `L_MIN`, `L_FACTOR`, `D_RATIO`, `D_TOTAL` | Sequence-quality-control heuristics (ES only). |

At least one of `align_dir` / `fasta_dir` must be set. The parser validates that the configured path exists and contains at least one `*.fa` / `*.fasta` file before the protocol runs.

## Output

All output files are written to an `output/` directory created in the current working directory.

- **ES:** `output/ES_<basename-of-align_dir>.out.txt`. Tab-separated; first column is the IDR alignment file name, remaining columns are `<feature>_meanZ` values in a fixed order.
- **FS:** `output/FS_<fasta-filename>.out.txt`. Tab-separated; first column is the sequence ID (FASTA header), remaining columns are `<feature>_meanZ` values in the same fixed order.

## Authors

- Iva Pritisanac (iva.pritisanac[at]helmholtz-munich.de), Helmholtz Munich
- Alan Moses (alan.moses[at]utoronto.ca), University of Toronto
- Julie Forman-Kay (forman[at]sickkids.ca), The Hospital for Sick Children

## Known bugs

- `random_int` in `src/core/es_pw_sim.py`: if amino-acid probabilities are supplied at limited precision (e.g. 2 or 4 decimal places) the function can return `None` instead of an integer. **Current workaround:** supply amino-acid frequencies at higher precision in `AA_COMPOSITION.txt`. **Planned:** raise an exception upstream.
