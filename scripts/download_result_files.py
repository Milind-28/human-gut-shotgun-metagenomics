#!/usr/bin/env python3
"""
Taxonomic Profile Downloader
-----------------------------
Downloads taxonomic profile TSV files from the HMP2/IBDMDB data portal
for each run listed in SraRunTable.csv.

Steps:
  1. Reads SraRunTable.csv
  2. Builds download URL using the "Library Name" column
  3. Downloads the file and renames it to the "Run" accession (SRR...)
  4. Moves each file into a disease-category folder:
       CD      -> Crohn's disease
       UC      -> ulcerative colitis
       NonIBD  -> empty / not specified
"""

import os
import csv
import time
import logging
import requests
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────

CSV_FILE       = "SraRunTable.csv"          # path to the metadata table
OUTPUT_DIR     = "functional_profiles"       # root output directory
BASE_URL       = (
    "https://g-227ca.190ebd.75bc.data.globus.org"
    "/ibdmdb/products/HMP2/MGX/func_profile3_WGS"
)

RETRY_LIMIT    = 3       # number of download retries on failure
RETRY_DELAY    = 5       # seconds to wait between retries
REQUEST_TIMEOUT = 120    # seconds before a single request times out

# ── Logging ───────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Helpers ───────────────────────────────────────────────────────────────────

def disease_folder(host_disease: str) -> str:
    """Map a Host_disease value to its output folder name."""
    val = host_disease.strip()
    if val.lower() == "crohn's disease":
        return "CD"
    elif val.lower() == "ulcerative colitis":
        return "UC"
    else:                        # empty or anything else
        return "NonIBD"


def download_file(url: str, dest: Path) -> bool:
    """
    Download *url* to *dest*.  Retries up to RETRY_LIMIT times.
    Returns True on success, False on permanent failure.
    """
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            log.info("  Attempt %d/%d: GET %s", attempt, RETRY_LIMIT, url)
            response = requests.get(url, timeout=REQUEST_TIMEOUT, stream=True)

            if response.status_code == 404:
                log.warning("  404 Not Found – skipping this file.")
                return False

            response.raise_for_status()

            dest.parent.mkdir(parents=True, exist_ok=True)
            with open(dest, "wb") as fh:
                for chunk in response.iter_content(chunk_size=1024 * 256):
                    fh.write(chunk)

            log.info("  Saved → %s", dest)
            return True

        except requests.RequestException as exc:
            log.warning("  Download error: %s", exc)
            if attempt < RETRY_LIMIT:
                log.info("  Retrying in %d s …", RETRY_DELAY)
                time.sleep(RETRY_DELAY)

    log.error("  All %d attempts failed for %s", RETRY_LIMIT, url)
    return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    csv_path = Path(CSV_FILE)
    if not csv_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {csv_path}")

    # Create category folders up-front
    root = Path(OUTPUT_DIR)
    for folder in ("CD", "UC", "NonIBD"):
        (root / folder).mkdir(parents=True, exist_ok=True)

    # Summary counters
    total = success = skipped = failed = 0

    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        # Validate required columns
        required = {"Run", "Library Name", "Host_disease"}
        missing  = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"CSV is missing columns: {missing}")

        for row in reader:
            total += 1
            run_id       = row["Run"].strip()
            library_name = row["Library Name"].strip()
            host_disease = row["Host_disease"]          # may be empty

            if not run_id or not library_name:
                log.warning("Row %d: empty Run or Library Name – skipping.", total)
                skipped += 1
                continue

            url      = f"{BASE_URL}/{library_name}_humann3.tar.bz2"
            folder   = disease_folder(host_disease)
            dest     = root / folder / f"{run_id}_humann3.tar.bz2"

            log.info(
                "[%s]  Library=%s  Disease=%s  → %s/",
                run_id, library_name, host_disease or "<empty>", folder,
            )

            if dest.exists():
                log.info("  Already exists – skipping download.")
                skipped += 1
                continue

            ok = download_file(url, dest)
            if ok:
                success += 1
            else:
                failed += 1

    # ── Summary ──
    log.info("=" * 60)
    log.info("Download complete.")
    log.info("  Total rows  : %d", total)
    log.info("  Downloaded  : %d", success)
    log.info("  Skipped     : %d", skipped)
    log.info("  Failed      : %d", failed)
    log.info("Files saved under: %s/", OUTPUT_DIR)
    log.info("  CD/     → Crohn's disease")
    log.info("  UC/     → Ulcerative colitis")
    log.info("  NonIBD/ → No disease / not specified")


if __name__ == "__main__":
    main()
