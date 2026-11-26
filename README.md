# RepBridge

RepBridge is a read-level, assembly-independent assay for the integrity of a
known antimicrobial resistance (AMR) module. It extracts a small target window
around two annotated anchor genes, maps short reads to that window, and counts
paired-end templates that physically connect the two anchors across their
junction.

RepBridge does **not** reconstruct plasmids and does **not** discover unknown
modules. It answers a narrower question: given a known module and its
annotation, is there physical template-level evidence that the two anchor
genes are linked in the original reads?

## Status

Manuscript under review. Patent application pending. This repository is
private; do not redistribute without permission.

## Requirements

- Linux or WSL2 (the workflow calls `bwa`, `samtools` and uses POSIX paths)
- Python 3.9 or newer
- `bwa` and `samtools` available on `PATH`
- Python packages: `numpy`, `pysam`, `scipy` (see `requirements.txt`)

## Installation

```bash
conda create -n repbridge python=3.11 -y
conda activate repbridge
conda install -c bioconda bwa samtools -y
pip install -r requirements.txt
```

## Quick start

Paired-end FASTQ:

```bash
python RepBridge_v1.py run \
  --sample SRR10094357 \
  --reads SRR10094357_1.fastq.gz SRR10094357_2.fastq.gz \
  --reference CP047975.1.fasta \
  --genbank CP047975.1.gb \
  --gene1 ble-MBL \
  --gene2 blaNDM-1 \
  --out results/SRR10094357 \
  --threads 8
```

The reference and annotation are user-supplied. For the Tn125 example, the
public records are NCBI accessions **CP047975.1** (pDETAB2) and its GenBank
annotation; for the mcr-1/pap2 example, **OM839890.1** (pEC1271-mcr).

Existing sorted and indexed BAM:

```bash
python RepBridge_v1.py run \
  --sample SAMPLE \
  --bam SAMPLE.sorted.bam \
  --reference reference.fasta \
  --genbank annotation.gb \
  --gene1 ble-MBL --gene2 blaNDM-1 \
  --out results/SAMPLE
```

Custom coordinates, without a GenBank/GFF file:

```bash
python RepBridge_v1.py run \
  --sample SAMPLE \
  --reads R1.fastq.gz R2.fastq.gz \
  --reference reference.fasta \
  --gene1-coords 61012-61377:- \
  --gene2-coords 61381-62193:- \
  --out results/SAMPLE
```

Coordinates are 1-based, inclusive, in the form `start-end:strand`.

## Panel mode

`panel` runs many samples from a TSV sheet:

```text
sample	r1	r2	bam
SRR10094357	/path/R1.fastq.gz	/path/R2.fastq.gz	
SAMPLE2		/path/SAMPLE2.sorted.bam
```

```bash
python RepBridge_v1.py panel \
  --samples samples.tsv \
  --reference reference.fasta \
  --genbank annotation.gb \
  --gene1 ble-MBL --gene2 blaNDM-1 \
  --out results/panel \
  --threads 8
```

Either `r1/r2` or `bam` must be provided for each row.

## Outputs

For each sample, `--out` receives:

- `summary.json` — all parameters, depths, linkage counts and verdict
- `summary.tsv` — one-row summary
- `report.md` — human-readable report
- `linkage_pairs.tsv` — supporting templates and junction coordinates
- `depth_regions.tsv` — depth statistics per region
- `run.log` — commands and timings
- `*.sorted.bam` / `.bai` — reads mapped to the target window

The key fields are `linkage.union_junction_support`, the two anchor median
depths, the junction median depth, the span spike ratio and `verdict`.

## Verdicts

- `positive`: both anchor median depths >= 5x, junction support >= 10
  templates, junction median depth >= 5x, and span spike ratio <= 5.
- `weak`: at least one supporting template, but not all positive criteria are
  met (for example, a coverage-spike flag or a depth shortfall).
- `absent`: both anchor median depths < 5x.
- `inconclusive`: none of the above.

For single-end reads, `present_single_end` is reported when either anchor
reaches the depth threshold; linkage cannot be assessed without pairs.

## Short-junction boundary correction

For short intergenic gaps, a read that merely touches the junction can be
counted as spanning it. Enable the optional overhang filter:

```bash
python RepBridge_v1.py run ... --min-overhang 10 --short-gap 10
```

The filter applies only when the intergenic gap is at most `--short-gap` bp;
a template is counted only if its alignment has at least `--min-overhang`
aligned bases on each side of the junction. The default is `0` (disabled).

## Reproducibility notes

- All mapping and filtering parameters are written to `summary.json`.
- Insert-size bounds are estimated from the same sample (median +/- 4 sigma,
  where sigma = 1.4826 x MAD).
- No reads, BAMs, reference genomes or patient data are included in this
  repository. Public accessions and the manuscript's data-availability
  statement should be used to obtain input data.
- Scripts assume `bwa` and `samtools` are installed; no binaries are bundled.

## License

Copyright (c) 2026 Wang Guanglin, Kang Xiuwen, Hu Rong. All rights reserved.
No permission is granted to use, copy, modify, distribute or commercialize
this software without written permission. Licensing will be finalized after
the pending patent application is filed.

## Contact

Hu Rong, 760020240051@xzhmu.edu.cn
Kang Xiuwen, 760020220009@xzhmu.edu.cn
