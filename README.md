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

`run_es.py` now takes **two arguments**: the input file, and the protocol to execute (`ES` or `FS`). It validates the protocol name, instantiates `EVSign`, checks that the corresponding path variable is set in the input file, and then calls either `compute_es_dir()` or `compute_fs_seq()`.

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
