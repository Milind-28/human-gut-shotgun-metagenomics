#!/usr/bin/env python3
import os
import shutil
import subprocess
import time
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    raise SystemExit(
        "Error: 'tqdm' library is missing. Please run 'pip install tqdm' before executing this script."
    )


def check_dependencies():
    """Ensure the required SRA Toolkit binaries are accessible."""
    for tool in ["prefetch", "fasterq-dump"]:
        if not shutil.which(tool):
            raise SystemExit(
                f"Error: '{tool}' is not installed or not in your PATH. "
                "Please install SRA Toolkit before running this script."
            )


def compress_file(file_path):
    """Compresses a file using pigz if available, otherwise falls back to gzip."""
    if shutil.which("pigz"):
        cmd = ["pigz", str(file_path)]
    else:
        cmd = ["gzip", str(file_path)]

    subprocess.run(
        cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


def process_cohort(acc_file_path, cohort_type, base_raw_dir):
    """Downloads SRA samples with an integrated checkpoint skip system."""
    target_dir = base_raw_dir / cohort_type
    target_dir.mkdir(parents=True, exist_ok=True)

    if not acc_file_path.exists():
        print(f"\n❌ Error: Accession file not found at {acc_file_path}")
        return

    with open(acc_file_path, "r") as f:
        sra_ids = [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

    total_samples = len(sra_ids)
    if total_samples == 0:
        print(f"\n⚠️ No accessions found in {acc_file_path.name}")
        return

    print("\n" + "=" * 65)
    print(
        f"🚀 Initializing Pipeline for {cohort_type} Cohort ({total_samples} Samples)"
    )
    print("=" * 65)

    with tqdm(
        total=total_samples,
        desc=f"📊 Analyzing {cohort_type}",
        unit="sample",
        bar_format="{l_bar}{bar:30}{r_bar}{bar:-10b}",
    ) as pbar:

        for sra_id in sra_ids:
            pbar.set_postfix_str(f"Checking: {sra_id}")

            # -----------------------------------------------------------------
            # 🔍 CHECKPOINT METHOD
            # -----------------------------------------------------------------
            # Define the expected final outputs of a successful run
            expected_r1 = target_dir / f"{sra_id}_1.fastq.gz"
            expected_r2 = target_dir / f"{sra_id}_2.fastq.gz"

            # If both compressed fastq files exist, skip this sample completely
            if expected_r1.exists() and expected_r2.exists():
                tqdm.write(
                    f"⏭️  [CHECKPOINT] {sra_id} fastq.gz files found. Skipping download..."
                )
                pbar.update(1)
                continue
            # -----------------------------------------------------------------

            pbar.set_postfix_str(f"Processing: {sra_id}")
            try:
                # Step 1: Prefetch SRA file
                prefetch_cmd = [
                    "prefetch",
                    sra_id,
                    "--output-directory",
                    str(target_dir),
                ]
                subprocess.run(
                    prefetch_cmd,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                sra_file_path = target_dir / sra_id / f"{sra_id}.sra"
                if not sra_file_path.exists():
                    sra_file_path = target_dir / f"{sra_id}.sra"

                # Step 2: Extract via fasterq-dump
                fasterq_cmd = [
                    "fasterq-dump",
                    "--split-files",
                    str(sra_file_path),
                    "--outdir",
                    str(target_dir),
                    "--temp",
                    str(target_dir),
                ]
                subprocess.run(
                    fasterq_cmd,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                # Step 3: Compress FASTQ files
                r1_file = target_dir / f"{sra_id}_1.fastq"
                r2_file = target_dir / f"{sra_id}_2.fastq"

                if r1_file.exists() and r2_file.exists():
                    compress_file(r1_file)
                    compress_file(r2_file)

                    # Clean up temporary SRA files to preserve disk space
                    sra_folder = target_dir / sra_id
                    if sra_folder.exists() and sra_folder.is_dir():
                        shutil.rmtree(sra_folder)

            except subprocess.CalledProcessError:
                tqdm.write(f"❌ Pipeline failed unexpectedly for sample {sra_id}")

            pbar.update(1)


def main():
    check_dependencies()

    script_dir = Path(__file__).resolve().parent
    base_data_dir = script_dir.parent / "data"
    metadata_dir = base_data_dir / "metadata"
    base_raw_dir = base_data_dir / "raw"

    cd_accessions = metadata_dir / "cd_acc.txt"
    uc_accessions = metadata_dir / "uc_acc.txt"
    nonibd_accessions = metadata_dir / "nonibd_acc.txt"

    start_time = time.time()

    process_cohort(cd_accessions, "CD", base_raw_dir)
    process_cohort(uc_accessions, "UC", base_raw_dir)
    process_cohort(nonibd_accessions, "NonIBD", base_raw_dir)

    total_duration = (time.time() - start_time) / 60
    print("\n" + "=" * 65)
    print(f"🎉 Pipeline finished successfully in {total_duration:.2f} minutes!")
    print("=" * 65)


if __name__ == "__main__":
    main()