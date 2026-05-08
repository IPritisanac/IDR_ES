#!/usr/bin/env python3
"""Download the IDR-ome data deposit from Zenodo.

The repository on GitHub ships only the analysis code. The companion data
(map, clusters, alignments, supplementary datasets) lives on Zenodo
because it is too large to keep in a git repository. This script fetches
the deposit's archive(s) and unpacks them into a local folder you can
point the Streamlit app at via the ``IDR_ES_ZENODO`` environment
variable or the sidebar text field.

Usage
-----

    # Using the record ID baked into ZENODO_RECORD_ID below
    python download_zenodo.py --target /path/to/where/you/want/ZENODO

    # Or override the record explicitly (e.g. to pin to a specific version)
    python download_zenodo.py --record 20076448 --target ./ZENODO

    # Skip extraction (just download the archives)
    python download_zenodo.py --target ./ZENODO --no-extract

    # Pick specific files (substring match against filename)
    python download_zenodo.py --target ./ZENODO --only MAP DATASETS

The script uses only the Python standard library — no extra dependencies.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import urllib.error
import urllib.request
import zipfile
import tarfile
from typing import Iterable, List, Optional


# ---------------------------------------------------------------------------
# Published Zenodo deposit accompanying the paper.
#
# Defaults to the *concept* DOI 10.5281/zenodo.10812874 — Zenodo's stable,
# always-latest pointer to whichever version of the deposit is current.
# The Zenodo API resolves this ID to the latest version automatically, so
# the script will keep working if the deposit is ever updated. If you want
# to pin to a specific version, pass --record <version_id> instead, e.g.
#   --record 20076448   (v3, DOI 10.5281/zenodo.20076448)
# ---------------------------------------------------------------------------
ZENODO_RECORD_ID: str = "10812874"

ZENODO_API = "https://zenodo.org/api/records/{record_id}"


def fetch_record(record_id: str) -> dict:
    url = ZENODO_API.format(record_id=record_id)
    print(f"Querying {url} …")
    try:
        with urllib.request.urlopen(url) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        sys.exit(
            f"ERROR: Zenodo API returned HTTP {e.code} for record "
            f"{record_id}. Check that the record ID is correct and that "
            "the deposit is public.\n"
            f"  {url}"
        )
    except urllib.error.URLError as e:
        sys.exit(f"ERROR: could not reach Zenodo: {e.reason}")


def filter_files(files: List[dict], only: Optional[Iterable[str]]) -> List[dict]:
    if not only:
        return files
    needles = [n.lower() for n in only]
    kept = [
        f for f in files
        if any(n in f["key"].lower() for n in needles)
    ]
    if not kept:
        names = ", ".join(f["key"] for f in files)
        sys.exit(
            "ERROR: --only matched no files in the deposit.\n"
            f"  Available files: {names}"
        )
    return kept


def humanize(n_bytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n_bytes < 1024:
            return f"{n_bytes:.1f} {unit}"
        n_bytes /= 1024
    return f"{n_bytes:.1f} TB"


def download_file(url: str, dest: str, expected_size: Optional[int] = None) -> None:
    if os.path.exists(dest) and expected_size is not None and os.path.getsize(dest) == expected_size:
        print(f"  ✓ {os.path.basename(dest)} already present, skipping download.")
        return
    tmp = dest + ".part"
    print(f"  ↓ {os.path.basename(dest)}", end="", flush=True)
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as out:
        total = expected_size or int(resp.headers.get("Content-Length") or 0)
        seen = 0
        chunk = 1024 * 256
        next_tick = 0
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            out.write(buf)
            seen += len(buf)
            if total and seen >= next_tick:
                pct = 100 * seen / total
                print(f"  {pct:5.1f}%", end="\r", flush=True)
                next_tick += total // 50
    os.replace(tmp, dest)
    print(f"  ↓ {os.path.basename(dest)}  done ({humanize(os.path.getsize(dest))})")


def extract_archive(path: str, target: str) -> bool:
    """Unpack ``path`` into ``target``. Returns True on success."""
    if path.endswith(".zip"):
        with zipfile.ZipFile(path) as zf:
            zf.extractall(target)
        return True
    if path.endswith((".tar", ".tar.gz", ".tgz", ".tar.bz2")):
        with tarfile.open(path) as tf:
            tf.extractall(target)
        return True
    return False


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--record", default=ZENODO_RECORD_ID,
                   help=f"Zenodo record ID (default: {ZENODO_RECORD_ID})")
    p.add_argument("--target", required=True,
                   help="Local folder to download into (will be created).")
    p.add_argument("--only", nargs="*",
                   help="Only download files whose name contains one of these substrings.")
    p.add_argument("--no-extract", action="store_true",
                   help="Don't unpack downloaded archives.")
    p.add_argument("--keep-archives", action="store_true",
                   help="Keep .zip / .tar files after extraction (default: delete).")
    args = p.parse_args()


    target = os.path.abspath(args.target)
    os.makedirs(target, exist_ok=True)

    rec = fetch_record(args.record)
    title = rec.get("metadata", {}).get("title", "(no title)")
    files = rec.get("files", [])
    if not files:
        sys.exit("ERROR: deposit has no files.")
    files = filter_files(files, args.only)

    total_bytes = sum(f.get("size", 0) for f in files)
    print(f"Record: {title}")
    print(f"Files to fetch: {len(files)}  ({humanize(total_bytes)} total)")
    for f in files:
        print(f"  · {f['key']}  ({humanize(f.get('size', 0))})")
    print()

    for f in files:
        key = f["key"]
        url = f.get("links", {}).get("self") or f.get("links", {}).get("download")
        if not url:
            print(f"  ! skipping {key}: no download URL in API response")
            continue
        dest = os.path.join(target, key)
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
        download_file(url, dest, expected_size=f.get("size"))

        if not args.no_extract and dest.endswith((".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2")):
            print(f"  ⇪ extracting {os.path.basename(dest)} into {target}")
            try:
                if extract_archive(dest, target):
                    if not args.keep_archives:
                        os.remove(dest)
                        print(f"     removed {os.path.basename(dest)}")
            except Exception as e:
                print(f"     extraction failed: {e}")

    print()
    print(f"Done. Data folder: {target}")
    print("Point the Streamlit app at it with one of:")
    print(f"  export IDR_ES_ZENODO={target}")
    print( "  # or paste the path into the sidebar text field")


if __name__ == "__main__":
    main()
