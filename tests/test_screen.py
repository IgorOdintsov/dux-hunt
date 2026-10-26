"""End-to-end tests on small synthetic BAMs (no patient data)."""
import random
import shutil

import pysam
import pytest

from dux_hunt.cli import make_config
from dux_hunt.detector import DUX4_OFFSET, DUX4_REF, _rc, _ref_len, screen_fusion

pytestmark = pytest.mark.skipif(shutil.which("bwa") is None, reason="bwa not on PATH")

CONTIGS = [("chr4", 191_154_276), ("chr10", 135_534_747), ("chr19", 59_128_983),
           ("chrUn_gl000228", 129_120)]
TID = {name: i for i, (name, _) in enumerate(CONTIGS)}
JUNCTION = 42_799_250          # 0-based first DUX4 base after the CIC part of the read
DUX4_JUNCTION = 1_500          # 0-based offset of the fused DUX4 base in the mini-reference


def _dux4_seq():
    with open(DUX4_REF) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def _random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _read(name, tid, pos, seq, cigar, mate_tid, mate_pos, flag):
    r = pysam.AlignedSegment()
    r.query_name, r.query_sequence, r.flag = name, seq, flag
    r.reference_id, r.reference_start, r.cigarstring = tid, pos, cigar
    r.mapping_quality = 60
    r.next_reference_id, r.next_reference_start = mate_tid, mate_pos
    r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return r


def _write_bam(path, reads):
    unsorted = str(path) + ".unsorted.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": n, "LN": ln} for n, ln in CONTIGS]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    return str(path)


def _background(rng, n=60):
    """Properly paired reads across the CIC window."""
    reads = []
    for i in range(n):
        pos = rng.randint(42_773_000, 42_800_000)
        reads.append(_read(f"bg{i}", TID["chr19"], pos, _random_seq(rng, 100), "100M",
                           TID["chr19"], pos + 250, 0x1 | 0x2 | 0x20 | 0x40))
    return reads


def _fusion_reads(rng, n_disc=12, n_split=8):
    dux4 = _dux4_seq()
    reads = []
    for i in range(n_disc):  # mates land on the DUX4 scaffold
        pos = rng.randint(42_798_900, 42_799_100)
        reads.append(_read(f"disc{i}", TID["chr19"], pos, _random_seq(rng, 100), "100M",
                           TID["chrUn_gl000228"], 113_000 + i, 0x1 | 0x40))
    for i in range(n_split):  # CIC bases, then a soft clip of DUX4 sequence
        m = 50 + i
        clip = dux4[DUX4_JUNCTION:DUX4_JUNCTION + 100 - m]
        reads.append(_read(f"split{i}", TID["chr19"], JUNCTION - m,
                           _random_seq(rng, m) + clip, f"{m}M{100 - m}S",
                           TID["chr19"], JUNCTION + 300, 0x1 | 0x2 | 0x20 | 0x40))
    return reads


def test_helpers():
    assert _ref_len("50M3I20M2D10S") == 72
    assert _rc("ACGTN") == "NACGT"


def test_positive_bam(tmp_path):
    rng = random.Random(1)
    bam = _write_bam(tmp_path / "fusion.bam", _background(rng) + _fusion_reads(rng))
    result = screen_fusion(bam, make_config("chr"))
    assert result.status == "POSITIVE"
    assert result.discordant_pairs == 12
    assert result.split_reads >= 8
    assert result.breakpoint_a == f"chr19:{JUNCTION + 1}"
    assert result.breakpoint_b == f"chrUn_gl000228:{DUX4_JUNCTION + 1 + DUX4_OFFSET}"


def test_negative_bam(tmp_path):
    rng = random.Random(2)
    bam = _write_bam(tmp_path / "normal.bam", _background(rng))
    result = screen_fusion(bam, make_config("chr"))
    assert result.status == "NEGATIVE"
    assert result.discordant_pairs == result.split_reads == 0


def test_discordant_only_is_negative(tmp_path):
    """Discordant pairs without split-read support do not make a call (>=2 split reads)."""
    rng = random.Random(3)
    bam = _write_bam(tmp_path / "disc_only.bam",
                     _background(rng) + _fusion_reads(rng, n_disc=15, n_split=0))
    result = screen_fusion(bam, make_config("chr"))
    assert result.discordant_pairs == 15 and result.split_reads == 0
    assert result.status == "NEGATIVE"
