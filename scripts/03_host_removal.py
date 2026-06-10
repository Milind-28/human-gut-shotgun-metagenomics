import os
import shutil
import subprocess
import time
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    raise SystemExit("Error: 'tqdm' library is missing. Run 'pip install tqdm'.")


def check_dependencies():
    """Verify Bowtie2 installation visibility."""
    if not shutil.which("bowtie2"):
        raise SystemExit(
            "Error: 'bowtie2' is not accessible. Activate your Conda environment."
        )


def deplete_host_dna(cohort_type, base_data_dir, index_prefix):
    """Aligns reads to GRCh38 index and isolates unmapped microbial pairs."""
    trimmed_dir = base_data_dir / "trimmed" / cohort_type
    host_removed_dir = base_data_dir / "host_removed" / cohort_type
    host_removed_dir.mkdir(parents=True, exist_ok=True)

    forward_reads = sorted(list(trimmed_dir.glob("*_1_paired.fastq.gz")))
    total_samples = len(forward_reads)

    if total_samples == 0:
        print(f"⚠️ No trimmed input reads discovered for cohort: {cohort_type}")
        return

    print("\n" + "=" * 75)
    print(
        f"🧬 Subtracting Host Contamination: {cohort_type} Cohort ({total_samples} Samples)"
    )
    print("=" * 75)

    with tqdm(
        total=total_samples,
        desc=f"🧬 Aligning {cohort_type}",
        unit="sample",
        bar_format="{l_bar}{bar:30}{r_bar}",
    ) as pbar:
        for r1_path in forward_reads:
            sra_id = r1_path.name.split("_1_paired.fastq.gz")[0]
            r2_path = trimmed_dir / f"{sra_id}_2_paired.fastq.gz"

            # Checkpoint target verification
            final_r1 = host_removed_dir / f"{sra_id}_1_microbial.fastq.gz"
            final_r2 = host_removed_dir / f"{sra_id}_2_microbial.fastq.gz"

            if final_r1.exists() and final_r2.exists():
                tqdm.write(
                    f"skip  [CHECKPOINT] {sra_id} is already depleted. Skipping..."
                )
                pbar.update(1)
                continue

            pbar.set_postfix_str(f"Mapping: {sra_id}")
            bt2_log = host_removed_dir / f"{sra_id}_bowtie2.log"

            try:
                # Bowtie2 command structure optimized for space limitations
                bowtie2_cmd = [
                    "bowtie2",
                    "-x",
                    str(index_prefix),
                    "-1",
                    str(r1_path),
                    "-2",
                    str(r2_path),
                    # Output unmapped paired reads using % wildcard format
                    "--un-conc-gz",
                    str(host_removed_dir / f"{sra_id}_%_microbial.fastq.gz"),
                    "--very-sensitive",  # Maximize alignment accuracy
                    "--threads",
                    "4",  # Parallel process execution
                    "-S",
                    "/dev/null",  # Avoid massive intermediary .sam storage allocation
                ]

                with open(bt2_log, "w") as log:
                    subprocess.run(
                        bowtie2_cmd,
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=log,
                    )

            except subprocess.CalledProcessError as e:
                tqdm.write(
                    f"❌ Alignment engine failed for sample {sra_id}: {e}"
                )

            pbar.update(1)


def main():
    check_dependencies()
    script_dir = Path(__file__).resolve().parent
    base_data_dir = script_dir.parent / "data"
    index_prefix = base_data_dir / "reference" / "GRCh38_noalt_as"/ "GRCh38_noalt_as"

    if not index_prefix.with_suffix(".1.bt2").exists():
        raise SystemExit(
            f"Error: Bowtie2 index files missing at {index_prefix.parent}."
        )

    start_time = time.time()
    #deplete_host_dna("CD", base_data_dir, index_prefix)
    #deplete_host_dna("UC", base_data_dir, index_prefix)
    deplete_host_dna("NonIBD", base_data_dir, index_prefix)

    print("\n" + "=" * 75)
    print(
        f"🎉 Host Depletion Complete in {(time.time() - start_time)/60:.2f} minutes!"
    )
    print("=" * 75)


if __name__ == "__main__":
    main()