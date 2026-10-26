<p align="center">
  <img src="docs/dux-hunt-logo.svg" alt="DUX-Hunt" width="420">
</p>

# DUX-Hunt

**Rapid *CIC::DUX4* fusion detection from targeted DNA sequencing, using realignment to a canonical *DUX4* mini-reference.**

DUX-Hunt reads an aligned BAM from a targeted DNA capture panel (hg19 / GRCh37) and reports whether it carries a *CIC::DUX4* rearrangement, the defining event of CIC-rearranged sarcoma. It needs no realignment of the sample and no change to an existing pipeline, and it takes under a minute per sample on a standard workstation.

## Why

*DUX4* lies inside the D4Z4 macrosatellite repeat, near-identical 3.3-kb units on 4q35, 10q26 and unplaced hg19 scaffolds. Reads that span a *CIC::DUX4* junction therefore map to many equally good locations. They get mapping quality ≈ 0, their mates scatter across the repeat copies, and genome-wide structural-variant callers filter out or fragment the evidence.

DUX-Hunt avoids the problem by looking only at the *CIC* locus. Any mate that lands on a *DUX4* copy counts as *DUX4*, and soft-clipped junction sequence is realigned to **a single canonical *DUX4* unit**, where it has exactly one place to go.

## How it works

1. **Fetch** reads in the *CIC* window, `chr19:42,772,404-42,800,948` (hg19).
2. **Discordant pairs.** Count reads whose mate maps to 4q35 (chr4 ≥ 190.9 Mb), 10q26 (chr10 ≥ 135.4 Mb) or `chrUn_gl000228`. Only pairs clustered within ±500 bp are counted, and identical sequences are counted once.
3. **Split reads.** Soft clips of at least 15 bp are aligned with `bwa mem` to the bundled canonical *DUX4* mini-reference (`chrUn_gl000228:111,467-125,289`, 13.8 kb). A clip is kept only when **both** junction sides agree within ±10 bp of the modal breakpoint. This double-concordance filter removes capture-panel ligation artifacts.
4. **Short-clip rescue.** Clips of 5–14 bp are counted when they exactly match the reconstructed junction sequence, but only once at least 3 concordant long clips have confirmed that junction.
5. **Call.** A sample is **POSITIVE** when it has **≥ 10 supporting reads (discordant pairs + split reads), including ≥ 2 split reads**. The output also reports the *CIC* and *DUX4* breakpoints. These thresholds were established on UCSF500, the UCSF in-house targeted next-generation sequencing cancer panel. Other panels can differ in bait design, coverage depth and read length, so we recommend revalidating the thresholds on known positive and negative cases before using DUX-Hunt on other assays. All three call thresholds can be set on the command line: --min-total (supporting reads, default 10), --min-split (split reads, default 2) and --min-disc (discordant pairs, default 0).


The split-read requirement matters. In our validation, every confirmed positive had at least 2 split reads, while the two calls with ≤ 1 split read were negative on clinical review.

## Installation

DUX-Hunt needs Python ≥ 3.9 and [`bwa`](https://github.com/lh3/bwa) on your `PATH`.

```bash
# bwa, e.g. via conda / mamba
mamba install -c bioconda bwa

# DUX-Hunt
pip install git+https://github.com/IgorOdintsov/dux-hunt.git
```

The *DUX4* mini-reference and its BWA index are bundled with the package, so no reference genome download is needed.

## Usage

```bash
dux-hunt sample.bam                      # one sample
dux-hunt *.bam -o cohort.tsv             # many samples; appends to the TSV, so reruns resume
dux-hunt sample.bam --min-total 10 --min-split 2 --min-disc 0   # call thresholds (defaults shown)
```

The BAM must be coordinate-sorted and indexed (`.bai` next to it). Contigs can be named with or without the `chr` prefix. For a BAM served over HTTPS, pass its index with `--bai` (or `SAMTOOLS_BAI`); only the *CIC* region is streamed.

### Output

One row per BAM, tab-separated:

| column | meaning |
|---|---|
| `sample` | BAM file name without `.bam` |
| `fusion` | `CIC::DUX4` |
| `discordant_pairs` | clustered, de-duplicated pairs with a mate on a *DUX4* locus |
| `split_reads` | doubly concordant long clips + rescued short clips |
| `supporting_reads` | `discordant_pairs + split_reads` |
| `total_reads` | reads in the *CIC* window |
| `fusion_index` | supporting reads per 1,000 reads in the window |
| `breakpoint_A` | *CIC* breakpoint (`chr19:pos`) |
| `breakpoint_B` | *DUX4* breakpoint on `chrUn_gl000228` (canonical unit) |
| `status` | `POSITIVE` or `NEGATIVE` |

## Validation

We applied DUX-Hunt retrospectively to 27,185 BAMs from a 500-gene targeted DNA panel, covering 23,296 unique patients. It detected *CIC::DUX4* in 10 patients; positive calls had 12–799 supporting reads. On clinical review, both cases with at least 10 and 11 total supporting reads that had ≤ 1 split read were negative. Among all other patients, no one had more than 9 total supporting reads. Every *CIC* breakpoint fell within a 215-bp window and every *DUX4* breakpoint within a 699-bp segment of the canonical unit.

## Limitations

- hg19 / GRCh37 coordinates only. GRCh38 is not supported yet.
- Paired-end, targeted capture data with baits over *CIC*.
- Detects *CIC::DUX4* only. Other *CIC* partners (e.g. *NUTM1*, *FOXO4*) are not assessed.
- This is research software, not a validated in-vitro diagnostic. Confirm clinically actionable results with an orthogonal method.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

The tests build small synthetic BAMs (no patient data) and run the full pipeline, including `bwa`.

## Citation

Siegmund SE, Sirohi D, Mirchia K, Van Ziffle J, Devine WP, Odintsov I. *DUX-Hunt: Rapid CIC::DUX4 Fusion Detection from DNA Sequencing Using Canonical DUX4 Mini-Reference Realignment.* Association for Molecular Pathology (AMP) Annual Meeting, 2026.

## Acknowledgments

The code, tests and documentation were developed with the assistance of [Claude Code](https://claude.com/claude-code), Anthropic's AI coding tool. All algorithm design decisions, call thresholds and validation results were reviewed by the authors.

## License

[MIT](LICENSE)
