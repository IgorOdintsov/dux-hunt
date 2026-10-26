"""Command-line interface for DUX-Hunt."""
from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path
from typing import Optional

import click

from .detector import (
    FusionConfig, FusionResult, screen_fusion,
    DUX4_REF, DUX4_CONTIG, DUX4_OFFSET,
    detect_chr_prefix,
)

log = logging.getLogger(__name__)

TSV_HEADER = [
    "sample", "fusion", "discordant_pairs", "split_reads", "supporting_reads",
    "total_reads", "fusion_index", "breakpoint_A", "breakpoint_B", "status",
]

# hg19 / GRCh37 CIC locus scanned for discordant and soft-clipped reads
CIC_WINDOW = "19:42772404-42800948"


def make_config(chr_prefix: str) -> FusionConfig:
    """CIC::DUX4 screening configuration for a BAM with the given contig prefix."""
    p = chr_prefix
    dux4_contig = f"{p}{DUX4_CONTIG}"
    chr4 = f"{p}4"
    chr10 = f"{p}10"

    def disc_dux4(read) -> bool:
        # mate on a DUX4-containing locus: 4q35, 10q26 or the chrUn_gl000228 scaffold
        m = read.next_reference_name
        mp = read.next_reference_start or 0
        return (
            m == dux4_contig
            or (m == chr4 and mp >= 190_900_000)
            or (m == chr10 and mp >= 135_400_000)
        )

    return FusionConfig(
        name="CIC::DUX4",
        gene_a_chr=f"{p}19",
        query_window=f"{p}{CIC_WINDOW}",
        disc_filter=disc_dux4,
        miniref=DUX4_REF,
        miniref_offset=DUX4_OFFSET,
        miniref_chr=dux4_contig,
        disc_b_chr=dux4_contig,
    )


def _print_result(result: FusionResult) -> None:
    sep = "═" * 46
    click.echo(f"\n{sep}")
    click.echo(f"  Sample:            {result.sample}")
    click.echo(f"  Fusion:            {result.fusion}")
    click.echo(f"  Discordant pairs:  {result.discordant_pairs}")
    click.echo(f"  Split reads:       {result.split_reads}")
    click.echo(f"  Total reads:       {result.total_reads}")
    click.echo(f"  Fusion index:      {result.fusion_index:.2f} per 1000 reads")
    click.echo(f"  Breakpoint A:      {result.breakpoint_a}")
    click.echo(f"  Breakpoint B:      {result.breakpoint_b}")
    click.echo(f"  Status:            {result.status}")
    click.echo(sep)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("bams", nargs=-1, required=True, metavar="BAM [BAM ...]")
@click.option(
    "-o", "--output", default="dux_hunt_results.tsv", show_default=True,
    help="Output TSV. Appends to an existing file (resume-safe).",
)
@click.option(
    "--bai", default=None, envvar="SAMTOOLS_BAI",
    help="BAI index path for remote HTTPS BAMs. "
         "Also reads the SAMTOOLS_BAI environment variable.",
)
@click.option(
    "-t", "--min-total", default=10, show_default=True,
    help="Minimum supporting reads (discordant pairs + split reads) to call POSITIVE.",
)
@click.option(
    "-s", "--min-split", default=2, show_default=True,
    help="Minimum split reads to call POSITIVE.",
)
@click.option(
    "-d", "--min-disc", default=0, show_default=True,
    help="Minimum discordant pairs to call POSITIVE.",
)
@click.option("-v", "--verbose", is_flag=True, help="Verbose logging (DEBUG level).")
@click.version_option(package_name="dux-hunt")
def main(
    bams: tuple[str, ...],
    output: str,
    bai: Optional[str],
    min_total: int,
    min_split: int,
    min_disc: int,
    verbose: bool,
) -> None:
    """Detect CIC::DUX4 fusions from targeted DNA panel BAMs.

    BAMs may be local file paths or HTTPS URLs. For remote BAMs pass the local
    BAI path via --bai or the SAMTOOLS_BAI environment variable; htslib streams
    only the CIC region rather than downloading the whole file.

    Coordinates: hg19 / GRCh37 (with or without the "chr" prefix).
    Requires `bwa` on PATH.

    \b
    Examples:
      dux-hunt sample.bam
      dux-hunt *.bam -o cohort.tsv
    """
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
        stream=sys.stderr,
    )

    out_path = Path(output)
    write_header = not out_path.exists()

    first_bam = bams[0]
    try:
        chr_prefix = detect_chr_prefix(first_bam)
    except Exception as exc:
        click.echo(f"ERROR: cannot open {first_bam}: {exc}", err=True)
        sys.exit(1)

    config = make_config(chr_prefix)
    n_ok = n_failed = n_skipped = 0

    with out_path.open("a", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if write_header:
            writer.writerow(TSV_HEADER)

        for bam_path in bams:
            if not bam_path.startswith("https://") and not Path(bam_path).exists():
                click.echo(f"[WARN] BAM not found: {bam_path} — skipping.", err=True)
                n_skipped += 1
                continue

            try:
                result = screen_fusion(
                    bam_path, config,
                    min_total=min_total,
                    min_split=min_split,
                    min_disc=min_disc,
                    bai_path=bai,
                )
                writer.writerow([
                    result.sample,
                    result.fusion,
                    result.discordant_pairs,
                    result.split_reads,
                    result.discordant_pairs + result.split_reads,
                    result.total_reads,
                    f"{result.fusion_index:.2f}",
                    result.breakpoint_a,
                    result.breakpoint_b,
                    result.status,
                ])
                fh.flush()
                _print_result(result)
                n_ok += 1

            except Exception as exc:
                click.echo(f"[ERROR] {bam_path}: {exc}", err=True)
                log.debug("Traceback:", exc_info=True)
                n_failed += 1

    click.echo(
        f"\nDone — OK: {n_ok}  Skipped: {n_skipped}  Failed: {n_failed}\n"
        f"Results: {output}"
    )
