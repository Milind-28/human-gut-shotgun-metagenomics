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

try:
    import multiqc
except ImportError:
    raise SystemExit(
        "Error: 'multiqc' library is missing. Please run 'pip install multiqc' before executing this script."
    )


def check_dependencies():
    """Ensure both fastp and multiqc are installed and accessible."""
    for tool in ["fastp", "multiqc"]:
        if not shutil.which(tool):
            raise SystemExit(
                f"Error: '{tool}' is not installed or accessible in your PATH.\n"
                "Please verify your Conda environment is active."
            )


def process_cohort_fastp(cohort_type, base_data_dir):
    """Runs fastp for simultaneous quality profiling, adapter trimming, and filtering."""
    raw_cohort_dir = base_data_dir / "raw" / cohort_type
    trimmed_cohort_dir = base_data_dir / "trimmed" / cohort_type
    qc_cohort_dir = base_data_dir / "results" / "qc" / cohort_type

    # Ensure output structures exist
    trimmed_cohort_dir.mkdir(parents=True, exist_ok=True)
    qc_cohort_dir.mkdir(parents=True, exist_ok=True)

    # Gather all forward read files (*_1.fastq.gz)
    forward_reads = sorted(list(raw_cohort_dir.glob("*_1.fastq.gz")))
    total_samples = len(forward_reads)

    if total_samples == 0:
        print(f"⚠️ No raw reads found for cohort: {cohort_type}")
        return

    print("\n" + "=" * 70)
    print(
        f"⚡ Running fastp QC & Trimming for {cohort_type} ({total_samples} Samples)"
    )
    print("=" * 70)

    with tqdm(
        total=total_samples,
        desc=f"⚡ Processing {cohort_type}",
        unit="sample",
        bar_format="{l_bar}{bar:30}{r_bar}",
    ) as pbar:

        for r1_path in forward_reads:
            sra_id = r1_path.name.split("_1.fastq.gz")[0]
            r2_path = raw_cohort_dir / f"{sra_id}_2.fastq.gz"

            if not r2_path.exists():
                tqdm.write(
                    f"⚠️ Warning: Missing reverse pair for {sra_id}. Skipping..."
                )
                pbar.update(1)
                continue

            # Define final paired outputs
            final_out1 = trimmed_cohort_dir / f"{sra_id}_1_paired.fastq.gz"
            final_out2 = trimmed_cohort_dir / f"{sra_id}_2_paired.fastq.gz"

            # Checkpoint Validation: Skip if outputs already exist
            if final_out1.exists() and final_out2.exists():
                tqdm.write(
                    f"skip  [CHECKPOINT] {sra_id} is already processed. Skipping..."
                )
                pbar.update(1)
                continue

            pbar.set_postfix_str(f"Processing: {sra_id}")

            # Define paths for HTML and JSON quality report generation
            html_report = qc_cohort_dir / f"{sra_id}_fastp.html"
            json_report = qc_cohort_dir / f"{sra_id}_fastp.json"

            try:
                # Constructing the fastp execution command
                fastp_cmd = [
                    "fastp",
                    "-i",
                    str(r1_path),
                    "-I",
                    str(r2_path),
                    "-o",
                    str(final_out1),
                    "-O",
                    str(final_out2),
                    "-h",
                    str(html_report),
                    "-j",
                    str(json_report),
                    "--detect_adapter_for_pe",
                    "--qualified_quality_phred",
                    "15",
                    "--unqualified_percent_limit",
                    "40",
                    "--length_required",
                    "36",
                    "--thread",
                    "4",
                ]

                # Run fastp and suppress stdout/stderr to prevent progress bar distortion
                subprocess.run(
                    fastp_cmd,
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

            except subprocess.CalledProcessError as e:
                tqdm.write(
                    f"❌ fastp processing failed execution for sample {sra_id}. Error: {e}"
                )

            pbar.update(1)


def run_multiqc_aggregation(base_data_dir):
    """Compiles all individual fastp reports into a single cohesive MultiQC dashboard."""
    qc_results_dir = base_data_dir / "results" / "qc"
    report_output_dir = base_data_dir / "results"

    print("\n" + "=" * 70)
    print("📊 Compiling Integrated Quality Control Report via MultiQC")
    print("=" * 70)

    # Construct the MultiQC command
    multiqc_cmd = [
        "multiqc",
        str(qc_results_dir),
        "-o",
        str(report_output_dir),
        "-n",
        "multiqc_report.html",
        "--force",
    ]

    try:
        subprocess.run(multiqc_cmd, check=True)
        print("\n" + "=" * 70)
        print("🎉 MultiQC Dashboard compiled successfully!")
        print(f"📂 Report Location: {report_output_dir / 'multiqc_report.html'}")
        print("=" * 70)
    except subprocess.CalledProcessError as e:
        print(f"\n❌ MultiQC compilation failed. Error: {e}")


def main():
    check_dependencies()

    script_dir = Path(__file__).resolve().parent
    base_data_dir = script_dir.parent / "data"

    pipeline_start = time.time()

    # Step 1: Execute fastp trimming sequentially across cohorts
    process_cohort_fastp("CD", base_data_dir)
    process_cohort_fastp("UC", base_data_dir)
    process_cohort_fastp("NonIBD", base_data_dir)
    # Step 2: Automatically aggregate fastp outputs into a single dashboard
    run_multiqc_aggregation(base_data_dir)

    total_duration = (time.time() - pipeline_start) / 60
    print(
        f"\n🚀 Quality Control Pipeline finished in {total_duration:.2f} minutes!"
    )


if __name__ == "__main__":
    main()