# IDR_ES

Code and Streamlit explorer for **"A Functional Map of the Human Intrinsically Disordered Proteome"** (Pritišanac et al., PNAS 2026).

This repository contains:

1. **A Streamlit web app** (`streamlit_app/`) for interactive exploration of the human IDR-ome map and IDR functional predictions.

Use the app to:

    - find a protein/IDR and visualize the 144-feature Z-score profile
    - browse IDR clusters extracted from the IDR-ome map at any of the 0.4–0.8 correlation distance thresholds (or the manually selected clusters from Dataset S2)
    - query the Supplementary Datasets, including the FAIDR-derived per-IDR GO-term assignment matrix (i.e., S6, Tab D - GO term predictions at the IDR level)
      
4. **The Python analysis package** (`src/`, `run_es.py`) 

Compute the IDR Z-score signatures, using one of the two protocols:
   
   - **ES — Evolutionary Signatures** 
   For each molecular feature, a Z-score is computed by taking the difference of the mean of the molecular feature across orthologous IDRs and the mean expected under a null model with no evolutionary restraint, divided by the null-model standard deviation. Requires an input MSA per IDR (see provided examples in `idr_aln/`).
   
   - **FS — Feature Signatures**
   For each molecular feature, a Z-score is computed by taking the difference of the IDR's feature value minus the mean of that feature across all human IDRs, divided by its standard deviation across the human IDRome. No orthologous sequences are required.

The full dataset of pre-computed Z-scores and the accompanying analyses (~350 MB compressed) is hosted on Zenodo at [doi:10.5281/zenodo.10812874](https://doi.org/10.5281/zenodo.10812874).

---
## Quick start 

To use the Streamlit explorer

# 1. Clone the repository
```bash
git clone https://github.com/IPritisanac/IDR_ES.git
cd IDR_ES

# 2. Get the data (~350 MB → ~2 GB unpacked)
python download_zenodo.py --target ./ZENODO

# 3. Set up a clean Python env and install dependencies
conda create -n idr_es python=3.11 -y
conda activate idr_es
pip install -r streamlit_app/requirements.txt

# 4. Launch
cd streamlit_app
streamlit run app.py
```

The browser tab will open and load the explorer at `http://localhost:8501`.

If this does not work as outlined, jump to **[Troubleshooting](#troubleshooting)** below

---

## Accompanying data (Zenodo)

> **Pritišanac, I.** *Data repository associated with 'A Functional Map of the Human Intrinsically Disordered Proteome'.* Zenodo. [https://doi.org/10.5281/zenodo.10812874](https://doi.org/10.5281/zenodo.10812874)

The DOI above is the *concept* DOI, which resolves to the latest published version. The deposit contains ten archives:

| Archive | Size | What's in it |
|---|---|---|
| `ES_MAP.zip` | 12 MB | clustered IDR-ome map (`HUMAN_ES.txt`, `HUMAN_ES.gtr`, `HUMAN_ES.cdt` ) |
| `CLUSTERS_AUTO.zip` | 76 MB | per-threshold (0.4–0.8) cluster archives + `AUTO_GO_FEATS.xlsx` |
| `CLUSTERS_EXPLORE.zip` | 10 MB | 93 manually exported / exploratory clusters |
| `DATASETS.zip` | 3 MB | Supplementary Datasets S1–S8 |
| `IDR_ALN.zip` | 200 MB | 19,459 alignments of human IDRs to orthologous sequences|
| `IDROME_SEQUENCES.zip` | 10 MB | proteome + IDRome FASTA, based on SPOT-Disorder v1.0-derived IDR boundaries|
| `FAIDR_TSTATS.zip` | <1 MB | FAIDR t-statistics’ hierarchical clustering |
| `FAIDR_HIGH_AUC_PPV_GO.zip` | 22 MB | 148 high-quality FAIDR target files |
| `PROTEIN_GROUPS_FAIDR_TARGETS.zip` | 2 MB | FAIDR target groups |
| `TUTORIAL.zip` | <1 MB | Cluster3.0 / JavaTreeView tutorial PDF |

`download_zenodo.py` calls Zenodo's public API, downloads each archive and unpacks the `.zip` files into a target directory. 

Optional flags:

- `--no-extract` — keep archives without unpacking.
- `--only ES_MAP DATASETS CLUSTERS_AUTO IDROME` — pull a subset of files (e.g., skip the 200 MB `IDR_ALN.zip` archive).
- `--record <id>` — pin to a specific version

After download, the layout that the Streamlit app expects is set up:

    ZENODO/
        MAP/HUMAN_ES.txt              # IDR-ome Z-score map (CDT format)
        MAP/HUMAN_ES.gtr              # dendrogram from Cluster3.0
        CLUSTERS_AUTO/CLUSTERS_0p4.zip … 0p8.zip
        CLUSTERS_AUTO/AUTO_GO_FEATS.xlsx
        IDROME_SEQUENCES/UP000005640_9606_SPOTD_MIN_30AA.fasta
        IDR_ALN/<UniProt>_…_ALN_IDR_<start>_<end>.fa
        DATASETS/DatasetS1.xlsx … S8.xlsx
---

## Streamlit explorer

The app has three switchable pages:

- **Find a protein** — autocomplete by gene name / UniProt / IDR ID; visualize the IDR's 144-feature Z-score profile, the cluster it belongs to, and the IDR's amino-acid sequence.  Select **"Jump to this cluster →"** button to get straight to the cluster with the IDR of interest on the “Browse clusters” page.

- **Browse clusters** — use a sidebar radio to choose the cluster set: either *Automatic* (Cluster3.0 output at one of the chosen 0.4–0.8 correlation threshold, 0.7 by default) or *Selected (Dataset S2)* (the manually curated clusters featured in the Supplementary Figure 7 of the paper). 
Visualize the per cluster mean Z-score profile of features, a Z-score heatmap (IDR members × features, scrollable), a member table with gene/protein names, the overrepresented GO terms, and the significantly enriched (positively / negatively) features.

- **Supplementary datasets** — browse Datasets S1–S8 directly from the Zenodo `DATASETS.zip`. Includes a dedicated directional-query UI for **Dataset S6 Tab D** (the FAIDR per-IDR GO-term assignment matrix): pick a GO term and a consistency threshold to list all IDRs that meet the criteria, or pick an IDR to list all its predicted GO terms. 

### Install & run

```bash
conda create -n idr_es python=3.11 -y
conda activate idr_es
cd streamlit_app
pip install -r requirements.txt
streamlit run app.py
```

The app should work with Python 3.9+;  (tested with Python 3.11).

The app caches everything via `st.cache_data`. The first request for loading a given dataset can be slow (e.g. the FAIDR matrix from Dataset S6 takes ~25–30 s to parse). The subsequent requests are instant.

By default the app will look for the data folder at `<repo>/ZENODO`, falling back to `<repo>/../ZENODO` and `~/ZENODO`. 

To override this path set

```bash
export IDR_ES_ZENODO=/abs/path/to/ZENODO
streamlit run app.py
```

or paste the path into the left sidebar of the app as it appears in the browser

under “Data” -> “Zenodo data folder”

### Troubleshooting

Most common issues arise from broken or stale Python environments (e.g. long-lived conda `(base)` envs). To mitigate, create a fresh, project-specific env:

`conda create -n idr_es python=3.11`

Follow the below for more info on potential errors and how to mitigate them

#### macOS: `rosetta error: … _multiarray_umath.cpython-38-darwin.so`

`streamlit run` will open a browser tab that never loads, the message printed in the terminal will read:

```
rosetta error: Attachment of code signature supplement failed: 1
 .../_multiarray_umath.cpython-38-darwin.so.aot
zsh: abort   streamlit run app.py
```

Cause: the NumPy installed in your Python env was built for a different CPU architecture than the interpreter (e.g. x86_64 NumPy under Apple Silicon arm64 Python). This is a common state for `(base)` conda envs that have been in use for years.

To fix this, first install into a new env:

```bash
conda deactivate
conda create -n idr_es python=3.11 -y
conda activate idr_es
cd ~/path/to/IDR_ES/streamlit_app
pip install -r requirements.txt
streamlit run app.py
```

Verify version before launching:

```bash
python -c "import platform, numpy; print('arch:', platform.machine(), '| numpy:', numpy.__version__)"
```

The expected terminal output is e.g. `arch: arm64 | numpy: 2.x.x` (or `x86_64` on Intel Macs). If you do not see a `rosetta error`, `streamlit run app.py` will function as intended.

For Mac OS

If your terminal/iTerm itself is running under Rosetta:

Uncheck "Open using Rosetta" in Finder → Get Info on the terminal app. 

Restart the terminal. Recreate the env (see above).

#### `ModuleNotFoundError: No module named 'streamlit.cli'`

Reported as e.g.:

```
File "/opt/anaconda3/bin/streamlit", line 7, in <module>
    from streamlit.cli import main
ModuleNotFoundError: No module named 'streamlit.cli'
```

Cause: the `streamlit` script in `<env>/bin/` was generated by an old install (Streamlit ≤ 1.3); the package itself has since been upgraded and the entry point moved to `streamlit.web.cli`. The wrapper points at a module that no longer exists. This is common in long-lived conda `(base)` envs.

To fix, use either:

# 1. Bypass the stale wrapper (this works whenever streamlit is importable):

```bash
python -m streamlit run app.py
```

Or

# 2. Reinstall

```bash
pip install --force-reinstall --no-deps streamlit
streamlit run app.py
```

Alternatively, set up the new env `conda create -n idr_es python=3.11` as introduced above to sidestep these fixes.
 
#### `zsh: command not found: streamlit`

`pip install` finished, but the `streamlit` binary is not on `PATH`. 

Use either:

# 1. Always works as long as streamlit is importable:

```bash
python -m streamlit run app.py
```
Or 

# 2. Add the env's bin/ to PATH manually (or activate the env):

```bash
which python                # see where Python lives
ls "$(dirname "$(which python)")"/streamlit*    # confirm that streamlit exists alongside it
```

#### `pip install` crashes with `InvalidVersion: '4.0.0-unsupported'` (or similar) at the end

This is usually not an issue. Newer pip is strict about PEP 440 version numbers and can crash during its **post-install summary scan**, *after* the packages are already in place.

To check:

```bash
python -c "import streamlit; print(streamlit.__version__)"
```

If a version is printed to the terminal, the app should run as intended.

Try running:

`streamlit run app.py`. 

Using a clean env removes the issues with stale dist-info dirs and pip post-install errors.

---

## Computing new FS/ES using **The Python analysis package** 

To compute Z-score signatures from your own input files (independent of the precomputed data and the Streamlit explorer):

    python run_es.py <input_file> <ES|FS>

`run_es.py` takes two arguments — the input file and the protocol type (“evolutionary” or “feature” signatures)

Examples:

    # Evolutionary Signatures from a directory with IDR alignments
    python run_es.py src/utils/input_file.txt ES

    # Feature Signatures from a FASTA file with IDR sequences
    python run_es.py src/utils/input_file.txt FS

If the wrong protocol is requested for the configured input, the script will print an explanatory error and exit.

### Dependencies

NumPy, SciPy. Python 3.9+ version:

```bash
pip install numpy scipy
```
or 

`conda install numpy scipy`

### Input file

A single plain-text input file configures both protocols. An example is provided in `src/utils/input_file.txt`. The most relevant entries:

| Key | Meaning |
|---|---|
| `motifs_file` | Path to the motif definitions (`MOTIFS.txt`) |
| `exp_motifs_n_file` | Path to the pre-computed expected-motif-counts |
| `repeats_file` | Path to the repeats’ definitions |
| `aa_freq_file` | Path to the background amino-acid composition |
| `align_dir` | **ES protocol only.** Directory with one alignment file for each IDR sequence (`*.fa` / `*.fasta` format) |
| `fasta_dir` | **FS protocol only.** Directory containing a single FASTA file with all IDR sequences (`*.fa` / `*.fasta` format) |
| `use_indels` | `on`/`off` — toggle the simple indel model in the null model simulations (ES only) |
| `n_simulations` | Number of pairwise simulations per orthologous IDR sequence (ES only) |
| `REF_NUM` | Default reference-sequence in each alignment (ES only) |
| `REF_NAME` | Substring used to auto-pick the reference sequence by name (ES only) |
| `MIN_SD` | Floor on the per-feature standard deviation used in denominators for Z-score computations |
| `L_MIN`, `L_FACTOR`, `D_RATIO`, `D_TOTAL` | Sequence-quality-control heuristics (ES only) |

At least one of `align_dir` / `fasta_dir` must be set. The parser validates that the configured path exists and contains at least one `*.fa` / `*.fasta` file before the protocol runs.

### Output

All outputs are written to an `output/` directory in the current working directory:

- **ES:** `output/ES_<basename-of-align_dir>.out.txt`. Tab-separated; first column is the IDR alignment file name, remaining columns are `<feature>_meanZ` values in a fixed order.
- **FS:** `output/FS_<fasta-filename>.out.txt`. Tab-separated; first column is the sequence ID (FASTA header), remaining columns are `<feature>_meanZ` values in the same fixed order.

### Method details (ES)

- Reads one alignment file (FASTA, `*.fa` / `*.fasta`) per IDR from `align_dir`
- Selects a reference sequence from each alignment file either by index (`REF_NUM`) or by name/substring match (`REF_NAME`, e.g. `HUMAN`)
- Applies sequence quality-control heuristics (`L_MIN`, `L_FACTOR`, `D_RATIO`, `D_TOTAL`) to filter orthologous sequences
- Computes a proxy pairwise evolutionary distance (F81 model) between the reference IDR and each orthologous sequence (**no phylogenetic tree is used**)
- Builds a null distribution via `n_simulations` pairwise simulations under JC69-style amino-acid substitution, with an optional simple indel model
- Returns per-feature **mean Z-scores**. Per-feature variance Z-scores are currently disabled
- Skips alignments if the reference IDR is shorter than 30 aa or if fewer than 10 orthologous sequences are found

### Method details (FS)

- Reads a single FASTA file from `fasta_dir` (warns if multiple are present and uses the first alphabetically)
- De-gaps sequences, drops empty sequences and sequences containing unknown characters e.g. `X`
- Computes raw feature values and normalizes each feature by the global mean and standard deviation across all IDRs in the FASTA. Standard deviation is floored at `MIN_SD`
- Returns per-sequence **mean Z-scores** for all features

---

## Authors

- Iva Pritišanac (iva.pritisanac[at]helmholtz-munich.de), Helmholtz Munich
- Alan Moses (alan.moses[at]utoronto.ca), University of Toronto
- Julie Forman-Kay (forman[at]sickkids.ca), The Hospital for Sick Children

## Known bugs

- `random_int` in `src/core/es_pw_sim.py`: if amino-acid probabilities are supplied at limited precision (e.g. 2 or 4 decimal places) the function can return `None` instead of an integer. **Current workaround:** supply amino-acid frequencies at higher precision in `AA_COMPOSITION.txt`. **Planned:** raise an exception upstream.
