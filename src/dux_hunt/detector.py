"""Core CIC::DUX4 detection logic: discordant pairs, split reads, concordance filter."""
from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional

import pysam

log = logging.getLogger(__name__)

# ── Bundled mini-reference paths and genomic offsets ──────────────────────────

_DATA = Path(__file__).parent / "data"

DUX4_REF     = _DATA / "dux4_miniref.fa"
DUX4_CONTIG  = "Un_gl000228"   # without chr prefix — callers prepend as needed
DUX4_START   = 111_467
DUX4_END     = 125_289
DUX4_OFFSET  = DUX4_START - 1  # BWA 1-based pos + offset → chrUn_gl000228 coordinate


def detect_chr_prefix(bam_path: str) -> str:
    """Return 'chr' if the BAM uses chr-prefixed contigs, otherwise ''."""
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for sq in bam.header.get("SQ", []):
            if sq["SN"].startswith("chr"):
                return "chr"
    return ""

# Tuning constants
_MIN_CLIP_LONG       = 15     # bp; shorter clips go to rescue
_MIN_CLIP_SHORT      = 5      # bp; minimum for rescue matching
_MAX_CLIP_SHORT      = 14     # bp; maximum for rescue matching
_DISC_CLUSTER_WIN    = 500    # bp; discard discordant reads outside this window of the median
_CONCORDANCE_WIN     = 10     # bp; A-side and B-side modal-breakpoint tolerance
_RESCUE_GATE         = 3      # minimum concordant long clips required to enable short-clip rescue
_BWA_MIN_SCORE       = 30     # -T flag passed to bwa mem

_REF_OPS = frozenset("MDN=X")


# ── Helpers ────────────────────────────────────────────────────────────────────

def _ref_len(cigar: str) -> int:
    """Number of reference-consuming bases in a CIGAR string."""
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
               if op in _REF_OPS)


def _rc(seq: str) -> str:
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _sample_name(bam_path: str) -> str:
    """Derive sample name from a BAM path or HTTPS URL (basename without .bam)."""
    from urllib.parse import unquote
    base = os.path.basename(bam_path.rstrip("/"))
    base = unquote(base)  # decode %2B → + etc. from URL-encoded paths
    return base[:-4] if base.endswith(".bam") else base


def _open_bam(bam_path: str, bai_path: Optional[str] = None) -> pysam.AlignmentFile:
    kwargs: dict = {"mode": "rb"}
    if bai_path and Path(bai_path).exists():
        kwargs["index_filename"] = bai_path
    return pysam.AlignmentFile(bam_path, **kwargs)


def _fetch(bam_path: str, region: str,
           bai_path: Optional[str] = None) -> List[pysam.AlignedSegment]:
    """Return all reads in *region* as a list (safe for repeated re-use)."""
    with _open_bam(bam_path, bai_path) as bam:
        return list(bam.fetch(region=region))


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass
class FusionConfig:
    name: str                        # display name, e.g. "CIC::DUX4"
    gene_a_chr: str                  # A-side chromosome, e.g. "chr19"
    query_window: str                # samtools-style region for A-side, e.g. "chr19:42791000-42840000"
    disc_filter: Callable            # (AlignedSegment) -> bool — True if mate is on B-side
    miniref: Path                    # BWA mini-reference FASTA for B-side soft-clip alignment
    miniref_offset: int              # add to BWA 1-based pos to get genomic B-side coordinate
    miniref_chr: str                 # chromosome label for B-side output, e.g. "chrUn_gl000228"
    disc_b_chr: str                  # primary B-side chromosome for fallback mate-position lookup


@dataclass
class FusionResult:
    sample: str
    fusion: str
    discordant_pairs: int
    split_reads: int
    total_reads: int
    fusion_index: float
    breakpoint_a: str
    breakpoint_b: str
    status: str


# ── Discordant pairs ───────────────────────────────────────────────────────────

def _count_discordant(reads: list, disc_filter: Callable) -> int:
    """
    Count discordant read pairs that cluster within ±500 bp of their median POS.
    Genuine fusions produce a tight cluster; panel ligation artifacts are scattered.
    Reads with identical sequences are counted once (PCR-duplicate guard).
    """
    seen: set[str] = set()
    positions = []
    for r in reads:
        if r.is_unmapped or not disc_filter(r):
            continue
        seq = r.query_sequence or ""
        if seq in seen:
            continue
        seen.add(seq)
        positions.append(r.reference_start)
    if not positions:
        return 0
    positions.sort()
    median = positions[len(positions) // 2]
    return sum(1 for p in positions if abs(p - median) <= _DISC_CLUSTER_WIN)


# ── Long soft-clip extraction ──────────────────────────────────────────────────

def _extract_long_clips(reads: list, out_fa: Path) -> int:
    """
    Write ≥15 bp soft clips to a FASTA.
    Read name encodes position and CIGAR for downstream breakpoint derivation:
      {readname}##{1-based-SAM-POS}##{cigar}_L  (left clip)
      {readname}##{1-based-SAM-POS}##{cigar}_R  (right clip)
    Returns the number of clips written.
    """
    count = 0
    seen_clips: set[str] = set()
    with out_fa.open("w") as fh:
        for read in reads:
            if not read.cigarstring or "S" not in read.cigarstring:
                continue
            seq = read.query_sequence
            if seq is None:
                continue
            cigar = read.cigarstring
            pos = read.reference_start + 1  # convert to 1-based SAM POS

            m = re.match(r"^(\d+)S", cigar)
            if m:
                n = int(m.group(1))
                if n >= _MIN_CLIP_LONG:
                    clip = seq[:n]
                    if clip not in seen_clips:
                        seen_clips.add(clip)
                        fh.write(f">{read.query_name}##{pos}##{cigar}_L\n{clip}\n")
                        count += 1

            m = re.search(r"(\d+)S$", cigar)
            if m:
                n = int(m.group(1))
                if n >= _MIN_CLIP_LONG:
                    clip = seq[-n:]
                    if clip not in seen_clips:
                        seen_clips.add(clip)
                        fh.write(f">{read.query_name}##{pos}##{cigar}_R\n{clip}\n")
                        count += 1
    return count


# ── BWA alignment ──────────────────────────────────────────────────────────────

def _bwa_align(clips_fa: Path, miniref: Path, out_sam: Path) -> int:
    """Align clips_fa to miniref with bwa mem. Returns number of mapped reads."""
    result = subprocess.run(
        ["bwa", "mem", "-T", str(_BWA_MIN_SCORE), str(miniref), str(clips_fa)],
        capture_output=True,
        text=True,
        check=True,
    )
    out_sam.write_text(result.stdout)
    mapped = 0
    for line in result.stdout.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) > 4 and not (int(fields[1]) & 4):
            mapped += 1
    return mapped


# ── Double concordance filter ──────────────────────────────────────────────────

def _double_concordance(
    clips_fa: Path, bwa_sam: Path
) -> tuple[str, int, int, int, int]:
    """
    Filter clips by double concordance:
      A-side: CIGAR-derived breakpoint within ±10 bp of the modal A-side position.
      B-side: BWA mapping position within ±10 bp of the modal B-side position.

    Chimeric ligation artifacts from capture panels cluster at the probe edge
    (A-side concordant) but scatter across the DUX4 locus (B-side discordant).
    Genuine fusions satisfy both constraints.

    Returns:
        (template_seq, rescue_bp_min, rescue_bp_max, modal_a_bp, concordant_count)
    """
    # Read BWA SAM → {read_name: bwa_1based_pos} for mapped reads only
    bwa_pos: dict[str, int] = {}
    for line in bwa_sam.read_text().splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 5 or int(f[1]) & 4:
            continue
        bwa_pos[f[0]] = int(f[3])

    if not bwa_pos:
        return "", 0, 0, 0, 0

    names_mapped = set(bwa_pos)
    confirmed: list[tuple[str, int, str]] = []  # (clip_seq, a_bp, name)
    raw_positions: list[int] = []
    cur_name: Optional[str] = None
    cur_bp: Optional[int] = None

    for line in clips_fa.read_text().splitlines():
        if line.startswith(">"):
            cur_name = line[1:].strip()
            cur_bp = None
            if cur_name in names_mapped:
                try:
                    parts = cur_name.split("##")
                    pos = int(parts[1])
                    cigar_side = parts[2]
                    side = cigar_side.rsplit("_", 1)[1]
                    cigar = cigar_side.rsplit("_", 1)[0]
                    raw_positions.append(pos)
                    if side == "R":
                        cur_bp = pos + _ref_len(re.sub(r"\d+S$", "", cigar))
                    else:
                        cur_bp = pos
                except (IndexError, ValueError):
                    pass
        elif cur_bp is not None and line.strip():
            confirmed.append((line.strip().upper(), cur_bp, cur_name))  # type: ignore[arg-type]
            cur_bp = None

    if not confirmed:
        return "", 0, 0, 0, 0

    # A-side concordance: keep clips near the modal A-side breakpoint
    modal_a = Counter(c[1] for c in confirmed).most_common(1)[0][0]
    a_conc = [(s, b, n) for s, b, n in confirmed if abs(b - modal_a) <= _CONCORDANCE_WIN]

    # B-side concordance: within A-concordant clips, keep those near the modal BWA position
    b_pos_list = [bwa_pos[n] for _, _, n in a_conc if n in bwa_pos]
    if b_pos_list:
        modal_b = Counter(b_pos_list).most_common(1)[0][0]
        doubly = [
            (s, b, n) for s, b, n in a_conc
            if abs(bwa_pos.get(n, -999_999) - modal_b) <= _CONCORDANCE_WIN
        ]
    else:
        doubly = a_conc

    if not doubly:
        return "", 0, 0, 0, 0

    # Longest doubly-concordant clip becomes the rescue template
    doubly.sort(key=lambda x: len(x[0]), reverse=True)
    template = doubly[0][0]
    median_pos = sorted(raw_positions)[len(raw_positions) // 2]

    return template, median_pos - 200, median_pos + 200, modal_a, len(doubly)


# ── Short-clip rescue ──────────────────────────────────────────────────────────

def _rescue_short_clips(
    reads: list, template: str, bp_min: int, bp_max: int
) -> int:
    """
    Count 5–14 bp clips within the rescue window that are exact substrings
    of the junction template (or its reverse complement).
    Gate: only called when ≥3 concordant long clips already confirmed.
    """
    template = template.upper()
    template_rc = _rc(template)
    matched: set[str] = set()  # unique clip sequences

    for read in reads:
        pos = read.reference_start + 1
        if pos < bp_min or pos > bp_max:
            continue
        if not read.cigarstring or "S" not in read.cigarstring:
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        cigar = read.cigarstring

        m = re.match(r"^(\d+)S", cigar)
        if m:
            n = int(m.group(1))
            if _MIN_CLIP_SHORT <= n < _MIN_CLIP_LONG:
                clip = seq[:n].upper()
                if clip in template or _rc(clip) in template or clip in template_rc:
                    matched.add(clip)

        m = re.search(r"(\d+)S$", cigar)
        if m:
            n = int(m.group(1))
            if _MIN_CLIP_SHORT <= n < _MIN_CLIP_LONG:
                clip = seq[-n:].upper()
                if clip in template or _rc(clip) in template or clip in template_rc:
                    matched.add(clip)

    return len(matched)


# ── Breakpoint helpers ─────────────────────────────────────────────────────────

def _modal_pos(values: list[int], one_based: bool = True) -> Optional[int]:
    if not values:
        return None
    modal = Counter(values).most_common(1)[0][0]
    return modal + (1 if one_based else 0)


def _fallback_a_bp(reads: list, disc_filter: Callable, gene_a_chr: str) -> str:
    positions = [r.reference_start for r in reads
                 if not r.is_unmapped and disc_filter(r)]
    p = _modal_pos(positions)
    return f"{gene_a_chr}:{p}" if p else "NA"


def _fallback_b_bp(reads: list, disc_b_chr: str) -> str:
    positions = [
        r.next_reference_start
        for r in reads
        if not r.is_unmapped
        and r.next_reference_name is not None
        and r.next_reference_name == disc_b_chr
        and r.next_reference_start is not None
    ]
    p = _modal_pos(positions)
    return f"{disc_b_chr}:{p}" if p else "NA"


# ── Main screening function ────────────────────────────────────────────────────

def screen_fusion(
    bam_path: str,
    config: FusionConfig,
    min_total: int = 10,
    min_split: int = 2,
    min_disc: int = 0,
    bai_path: Optional[str] = None,
) -> FusionResult:
    """
    Detect one fusion type in one BAM and return a FusionResult.

    Steps:
      1. Count discordant pairs (clustered within ±500 bp of median)
      2. Extract ≥15 bp soft clips and BWA-align to the B-side mini-reference
      3. Apply double concordance filter (A-side ±10 bp AND B-side ±10 bp)
      4. Short-clip rescue (5–14 bp exact substring match, gated on ≥3 concordant)
      5. Classify: POSITIVE (discordant + split ≥ min_total AND split ≥ min_split
         AND discordant ≥ min_disc) / NEGATIVE
    """
    sample = _sample_name(bam_path)
    log.info("[%s] Analysing %s", config.name, sample)
    log.info("[%s] Query window: %s", config.name, config.query_window)

    reads = _fetch(bam_path, config.query_window, bai_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        clips_fa = tmp / "clips_long.fa"
        bwa_sam  = tmp / "long_aln.sam"

        # ── Step 1: discordant pairs ──────────────────────────────────────────
        disc = _count_discordant(reads, config.disc_filter)

        # ── Step 2: long split reads via BWA ─────────────────────────────────
        n_clips = _extract_long_clips(reads, clips_fa)
        split_long = concordant_long = 0
        template = ""
        bp_rescue_min = bp_rescue_max = 0
        modal_a_bp = 0
        b_side_bwa_pos = 0

        if n_clips > 0:
            try:
                mapped = _bwa_align(clips_fa, config.miniref, bwa_sam)
            except subprocess.CalledProcessError as exc:
                log.warning("[%s] bwa mem failed: %s", config.name, exc)
                mapped = 0

            if mapped > 0:
                # B-side: modal BWA position → genomic coordinate
                bwa_positions = []
                for line in bwa_sam.read_text().splitlines():
                    if line.startswith("@"):
                        continue
                    f = line.split("\t")
                    if len(f) > 3 and not (int(f[1]) & 4):
                        bwa_positions.append(int(f[3]))
                if bwa_positions:
                    modal_bwa = Counter(bwa_positions).most_common(1)[0][0]
                    b_side_bwa_pos = modal_bwa + config.miniref_offset

                # ── Step 3: double concordance filter ─────────────────────────
                (template, bp_rescue_min, bp_rescue_max,
                 modal_a_bp, concordant_long) = _double_concordance(clips_fa, bwa_sam)

        split_long = concordant_long
        log.info(
            "[%s] Long clips: %d concordant at breakpoint (of %d BWA-mapped)",
            config.name, split_long, n_clips,
        )

        # ── Step 4: short-clip rescue ─────────────────────────────────────────
        split_short = 0
        if (template and len(template) >= 8
                and bp_rescue_min > 0
                and concordant_long >= _RESCUE_GATE):
            split_short = _rescue_short_clips(reads, template, bp_rescue_min, bp_rescue_max)

        split_total = split_long + split_short
        log.info(
            "[%s] Split reads: %d long + %d short = %d total",
            config.name, split_long, split_short, split_total,
        )

        # ── Total reads in window ─────────────────────────────────────────────
        total = max(len(reads), 1)
        fusion_index = round((disc + split_total) / total * 1000, 2)

        # ── Breakpoints ───────────────────────────────────────────────────────
        bp_a = (
            f"{config.gene_a_chr}:{modal_a_bp}"
            if modal_a_bp > 0
            else _fallback_a_bp(reads, config.disc_filter, config.gene_a_chr)
        )
        bp_b = (
            f"{config.miniref_chr}:{b_side_bwa_pos}"
            if b_side_bwa_pos > 0
            else _fallback_b_bp(reads, config.disc_b_chr)
        )

        # ── Step 5: classify ──────────────────────────────────────────────────
        evidence = disc + split_total
        if evidence >= min_total and split_total >= min_split and disc >= min_disc:
            status = "POSITIVE"
        else:
            status = "NEGATIVE"

    log.info(
        "[%s] disc=%d splits=%d total=%d index=%.2f status=%s",
        config.name, disc, split_total, total, fusion_index, status,
    )

    return FusionResult(
        sample=sample,
        fusion=config.name,
        discordant_pairs=disc,
        split_reads=split_total,
        total_reads=total,
        fusion_index=fusion_index,
        breakpoint_a=bp_a,
        breakpoint_b=bp_b,
        status=status,
    )
