#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RepBridge v1.1.1
================
Assembly-independent, reference-guided detection of ARG module integrity
from short-read WGS data.

Core idea
---------
Do not trust the assembly. Map reads to a small, annotation-defined target
window and ask two questions:

  1. Is there read coverage over both anchor genes?
  2. Is there template-level physical linkage across their junction
     (single read spanning, or mate pair spanning)?

The answer is a verdict: positive / weak / absent / inconclusive.

Features
--------
- Annotation-driven coordinates (GenBank or GFF3); no hardcoded intervals.
- Automatic target-window extraction (gene1/gene2 +/- flank).
- Paired-end or single-end; FASTQ or existing BAM.
- Template-level counting (not read-level), MAPQ/proper-pair/insert filters,
  secondary/supplementary/duplicate/QC-fail exclusion.
- Depth, coverage spike ratio, Poisson CI.
- Optional negative-control gene pair.
- Single-sample and batch (panel) modes.
- Reports every parameter and source path for reproducibility.

Usage
-----
  # paired FASTQ, default CP047975.1 Tn125 module
  python RepBridge_v1.py run \
    --sample SRR21472814 \
    --reads R1.fastq.gz R2.fastq.gz \
    --out results/SRR21472814

  # existing BAM (sorted + indexed, against --reference)
  python RepBridge_v1.py run --sample X --bam X.sorted.bam --out results/X

  # batch
  python RepBridge_v1.py panel --samples samples.tsv --out results/panel

  # custom genes / coordinates
  python RepBridge_v1.py run ... --gene1 ble-MBL --gene2 blaNDM-1
  python RepBridge_v1.py run ... --gene1-coords 61012-61377:- \
                                 --gene2-coords 61381-62193:-

samples.tsv columns: sample  r1  r2  [bam]
  - r1/r2 may be empty if bam is provided.
"""

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np
import pysam

VERSION = "1.1.2"
DEFAULT_REF = "/mnt/d/RepBridge_Data/reference/CP047975.1.fasta"
DEFAULT_GBK = "/mnt/d/ENA/RepBridge_v3/CP047975.1.gb"
DEFAULT_GENE1 = "ble-MBL"
DEFAULT_GENE2 = "blaNDM-1"


# ============================================================================
# logging / shell
# ============================================================================
class Logger:
    def __init__(self, path=None):
        self.path = path
        self.fh = open(path, "w", encoding="utf-8") if path else None

    def __call__(self, msg=""):
        print(msg, flush=True)
        if self.fh:
            self.fh.write(str(msg) + "\n")
            self.fh.flush()

    def close(self):
        if self.fh:
            self.fh.close()


def run_cmd(cmd, log, check=True, shell=False):
    shown = cmd if isinstance(cmd, str) else " ".join(shlex.quote(x) for x in cmd)
    log(f"  $ {shown}")
    t0 = time.time()
    res = subprocess.run(cmd, shell=shell, capture_output=True, text=True,
                         executable="/bin/bash" if shell else None)
    if check and res.returncode != 0:
        log((res.stdout or "")[-1000:])
        log((res.stderr or "")[-2500:])
        raise RuntimeError(f"cmd failed ({res.returncode}): {shown}")
    log(f"    完成 ({time.time()-t0:.1f}s)")
    return res


def require_tools(mapping=True):
    need = ["samtools"] + (["bwa"] if mapping else [])
    missing = [t for t in need if not shutil.which(t)]
    if missing:
        raise SystemExit("PATH 缺少: " + ", ".join(missing))


# ============================================================================
# annotation
# ============================================================================
@dataclass
class Feature:
    seqid: str
    ftype: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    strand: str
    gene: str = ""
    product: str = ""
    locus_tag: str = ""
    joined: bool = False

    def label(self):
        bits = [f"{self.seqid}:{self.start}-{self.end}({self.strand})"]
        if self.gene:
            bits.append(f"gene={self.gene}")
        if self.product:
            bits.append(f"product={self.product}")
        return " ".join(bits)


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def parse_genbank(path):
    feats = []
    seqid = None
    in_feat = False
    cur = None

    def flush():
        nonlocal cur
        if not cur:
            return
        head, lines = cur
        m = re.match(r"^ {5}(\S+)\s+(.+)$", head)
        if not m:
            cur = None
            return
        ftype, loc = m.group(1), m.group(2).strip()
        for extra in lines:
            if extra.strip().startswith("/"):
                break
            loc += extra.strip()
        nums = re.findall(r"(\d+)\.\.(\d+)", loc)
        if nums:
            start, end = int(nums[0][0]), int(nums[-1][1])
            strand = "-" if "complement" in loc else "+"
            qual = "\n".join(lines)
            def q(key):
                mm = re.search(r'/' + key + r'="([^"]+)"', qual)
                return mm.group(1) if mm else ""
            feats.append(Feature(seqid or "unknown", ftype, start, end, strand,
                                 gene=q("gene"), product=q("product"),
                                 locus_tag=q("locus_tag"),
                                 joined=("join" in loc)))
        cur = None

    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("LOCUS"):
                p = line.split()
                seqid = p[1] if len(p) > 1 else seqid
                continue
            if line.startswith("VERSION"):
                p = line.split()
                seqid = p[1] if len(p) > 1 else seqid
                continue
            if line.startswith("FEATURES"):
                in_feat = True
                continue
            if not in_feat:
                continue
            if line.startswith("ORIGIN") or line.startswith("CONTIG"):
                break
            m = re.match(r"^ {5}(\S+)\s+(.*)$", line)
            if m:
                flush()
                cur = (line, [])
                continue
            if cur is not None:
                cur[1].append(line)
    flush()
    return [f for f in feats if f.ftype in ("CDS", "gene", "mRNA")]


def parse_gff3(path):
    feats = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("CDS", "gene", "mRNA"):
                continue
            attrs = {}
            for item in p[8].split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = v
            feats.append(Feature(
                p[0], p[2], int(p[3]), int(p[4]), p[6],
                gene=attrs.get("gene", ""), product=attrs.get("product", ""),
                locus_tag=attrs.get("locus_tag", attrs.get("gene_id", ""))))
    return feats


def select_feature(features, query, seqid=None):
    pool = [f for f in features if f.ftype in ("CDS", "gene", "mRNA")]
    if seqid:
        pool = [f for f in pool if f.seqid == seqid]
    qn = _norm(query)
    variants = {qn}
    if qn.startswith("bla") and len(qn) > 3:
        variants.add(qn[3:])
    scored = []
    for f in pool:
        text = _norm(" ".join([f.gene, f.product, f.locus_tag]))
        score = 0
        for q in variants:
            if q and q == text:
                score = max(score, 3)
            elif q and q in text:
                score = max(score, 2)
        if score:
            scored.append((score, f))
    if not scored:
        raise SystemExit(f"No annotation matched gene {query!r}. "
                         f"Use --gene*-coords or check the annotation.")
    top = max(s for s, _ in scored)
    best = [f for s, f in scored if s == top]
    uniq = {}
    for f in best:
        uniq.setdefault((f.seqid, f.start, f.end, f.strand), f)
    if len(uniq) > 1:
        msg = [f"Gene {query!r} is ambiguous:"]
        msg += ["  " + f.label() for f in uniq.values()]
        raise SystemExit("\n".join(msg))
    return next(iter(uniq.values()))


def parse_coords(text, seqid):
    t = text.replace(",", "").strip()
    m = re.match(r"^(\d+)[:-](\d+)(?::([+-]))?$", t)
    if not m:
        raise SystemExit(f"Bad coordinates: {text!r}; use START-END[:STRAND]")
    start, end = int(m.group(1)), int(m.group(2))
    if start < 1 or end < start:
        raise SystemExit(f"Bad interval: {text!r}")
    return Feature(seqid, "manual", start, end, m.group(3) or "+")


# ============================================================================
# reference / target
# ============================================================================
def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(path, name, seq, width=70):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def build_target(ref_seqs, ref_id, f1, f2, outdir, flank):
    seq = ref_seqs[ref_id]
    win_s = max(1, min(f1.start, f2.start) - flank)
    win_e = min(len(seq), max(f1.end, f2.end) + flank)
    name = "RepBridge_target"
    path = os.path.join(outdir, "target.fasta")
    write_fasta(path, name, seq[win_s - 1:win_e])

    def shift(f):
        return Feature(name, f.ftype, f.start - win_s + 1,
                       f.end - win_s + 1, f.strand, f.gene, f.product,
                       f.locus_tag, f.joined)

    return path, name, shift(f1), shift(f2), win_s, win_e


def index_ref(path, log):
    if not all(os.path.exists(path + e)
               for e in (".bwt", ".pac", ".ann", ".amb", ".sa")):
        run_cmd(["bwa", "index", path], log)
    if not os.path.exists(path + ".fai"):
        run_cmd(["samtools", "faidx", path], log)


# ============================================================================
# mapping
# ============================================================================
def map_reads(ref, r1, r2, bam, sample, threads, log, force=False):
    if os.path.exists(bam) and os.path.exists(bam + ".bai") and not force:
        log(f"  ℹ️  复用 BAM: {bam}")
        return
    os.makedirs(os.path.dirname(bam), exist_ok=True)
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    tmp = f"/tmp/repbridge_sort_{os.getpid()}_{re.sub(r'[^A-Za-z0-9._-]', '_', sample)}"
    sort = f"samtools sort -@ {threads} -T {shlex.quote(tmp)} -o {shlex.quote(bam)} -"
    view = "samtools view -b -F 2308 -"
    if r2:
        cmd = (f"bwa mem -t {threads} -R {shlex.quote(rg)} "
               f"{shlex.quote(ref)} {shlex.quote(r1)} {shlex.quote(r2)} "
               f"| {view} | {sort}")
    else:
        cmd = (f"bwa mem -t {threads} -R {shlex.quote(rg)} "
               f"{shlex.quote(ref)} {shlex.quote(r1)} | {view} | {sort}")
    run_cmd(cmd, log, shell=True)
    run_cmd(["samtools", "index", bam], log)


def flagstat(bam):
    res = subprocess.run(["samtools", "flagstat", bam],
                         capture_output=True, text=True)
    txt = res.stdout

    def grab(pat):
        m = re.search(pat, txt, re.M)
        return int(m.group(1)) if m else None

    total = grab(r"^(\d+) \+ \d+ in total")
    mapped = grab(r"^(\d+) \+ \d+ mapped")
    proper = grab(r"^(\d+) \+ \d+ properly paired")
    return {"total": total, "mapped": mapped, "properly_paired": proper,
            "mapping_rate": (mapped / total) if total else None}


# ============================================================================
# statistics
# ============================================================================
def insert_bounds(sam, ref_id, mapq, max_insert=None, sample_limit=200000,
                  mad_mult=None):
    vals = []
    for r in sam.fetch(ref_id):
        if len(vals) >= sample_limit:
            break
        if (r.is_paired and r.is_proper_pair and not r.is_unmapped
                and not r.is_secondary and not r.is_supplementary
                and not r.is_duplicate and not r.is_qcfail
                and r.mapping_quality >= mapq and not r.mate_is_unmapped
                and r.next_reference_id == r.reference_id):
            t = abs(r.template_length)
            if t:
                vals.append(t)
    if not vals:
        return {"n": 0, "median": None, "mad": None,
                "lower": None, "upper": max_insert}
    arr = np.asarray(vals, dtype=np.int64)
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    sigma = 1.4826 * mad if mad > 0 else float(np.std(arr))
    if mad_mult is not None:
        lower = max(1.0, med - float(mad_mult) * max(mad, 1.0))
        upper = (float(max_insert) if max_insert
                 else med + float(mad_mult) * max(mad, 1.0))
    else:
        upper = float(max_insert) if max_insert else med + 4 * sigma
        lower = max(1.0, med - 4 * sigma)
    return {"n": len(vals), "median": med, "mad": mad,
            "lower": lower, "upper": upper}


def depth_stats(sam, ref_id, start, end, mapq):
    def cb(r):
        return (not r.is_unmapped and not r.is_secondary
                and not r.is_supplementary and not r.is_duplicate
                and not r.is_qcfail and r.mapping_quality >= mapq)

    cov = sam.count_coverage(ref_id, start - 1, end,
                             quality_threshold=0, read_callback=cb)
    arr = np.asarray([sum(x) for x in zip(*cov)], dtype=np.int64)
    if arr.size == 0:
        return {"n": 0, "median": 0, "mean": 0, "max": 0,
                "p95": 0, "spike_ratio": None}
    med = float(np.median(arr))
    mx = int(arr.max())
    return {"n": int(arr.size), "median": med, "mean": float(arr.mean()),
            "max": mx, "p95": float(np.percentile(arr, 95)),
            "spike_ratio": (mx / med) if med else None}


def overlaps(s, e, f):
    return s < f.end and e > (f.start - 1)


def analyze_linkage(sam, ref_id, f1, f2, mapq, bounds, max_insert=None,
                    min_overhang=0, short_gap=10):
    if f1.end <= f2.start:
        jl, jr = f1.end, f2.start
    elif f2.end <= f1.start:
        jl, jr = f2.end, f1.start
    else:
        jl = jr = (max(f1.start, f2.start) + min(f1.end, f2.end)) // 2

    gap_bp = max(0, jr - jl - 1)
    apply_boundary = bool(min_overhang and min_overhang > 0
                          and gap_bp <= short_gap)

    lo = max(0, min(f1.start, f2.start) - 2000)
    hi = max(f1.end, f2.end) + 2000
    by_name = defaultdict(list)
    for r in sam.fetch(ref_id, lo, hi):
        by_name[r.query_name].append(r)

    rows = []
    pair_any, pair_junc, single_junc = set(), set(), set()
    reject = Counter()
    passed = 0

    for name, recs in by_name.items():
        primary = [r for r in recs
                   if not r.is_secondary and not r.is_supplementary]
        if not primary:
            reject["no_primary"] += 1
            continue
        r = next((x for x in primary if x.is_read1), primary[0])
        if not r.is_paired:
            reject["not_paired"] += 1
            continue
        if not r.is_proper_pair:
            reject["not_proper_pair"] += 1
            continue
        if r.is_unmapped or r.mate_is_unmapped:
            reject["unmapped"] += 1
            continue
        if r.mapping_quality < mapq:
            reject["mapq_read"] += 1
            continue
        if r.is_duplicate or r.is_qcfail:
            reject["duplicate_or_qcfail"] += 1
            continue
        if r.next_reference_id != r.reference_id:
            reject["mate_other_ref"] += 1
            continue
        try:
            m = sam.mate(r)
        except Exception:
            reject["mate_fetch_failed"] += 1
            continue
        if m.mapping_quality < mapq:
            reject["mapq_mate"] += 1
            continue
        tlen = abs(r.template_length)
        if tlen == 0:
            reject["zero_tlen"] += 1
            continue
        if max_insert and tlen > max_insert:
            reject["insert_too_large"] += 1
            continue
        if bounds.get("upper") is not None and tlen > bounds["upper"]:
            reject["insert_outlier"] += 1
            continue
        passed += 1

        rs1, re1 = r.reference_start, r.reference_end
        rs2, re2 = m.reference_start, m.reference_end
        in1a, in1b = overlaps(rs1, re1, f1), overlaps(rs1, re1, f2)
        in2a, in2b = overlaps(rs2, re2, f1), overlaps(rs2, re2, f2)

        any_pair = (in1a and in2b) or (in1b and in2a)
        span_j = (min(rs1, rs2) <= jl and max(re1, re2) >= jr)
        single_j = ((rs1 <= jl and re1 >= jr) or
                    (rs2 <= jl and re2 >= jr))

        # boundary-artefact suppression for short intergenic gaps
        pair_left_ov = max(0, jl - min(rs1, rs2))
        pair_right_ov = max(0, max(re1, re2) - jr)
        single_left_ov = 0
        single_right_ov = 0
        if apply_boundary:
            # Proper pairs that span the junction are retained; the overhang
            # requirement is applied to single-read spans, which are the main
            # source of boundary artefacts at very short intergenic gaps.
            single_ok = False
            for s, e in ((rs1, re1), (rs2, re2)):
                if s <= jl and e >= jr:
                    l_ov = max(0, jl - s)
                    r_ov = max(0, e - jr)
                    if l_ov >= min_overhang and r_ov >= min_overhang:
                        single_ok = True
                        single_left_ov = max(single_left_ov, l_ov)
                        single_right_ov = max(single_right_ov, r_ov)
            single_j = single_j and single_ok

        if not (any_pair or span_j or single_j):
            continue

        if any_pair:
            pair_any.add(name)
        if span_j:
            pair_junc.add(name)
        if single_j:
            single_junc.add(name)
        rows.append({
            "template": name, "read_start": rs1, "read_end": re1,
            "mate_start": rs2, "mate_end": re2, "insert_size": tlen,
            "mapq_read": r.mapping_quality, "mapq_mate": m.mapping_quality,
            "pair_in_genes": int(any_pair),
            "pair_spans_junction": int(span_j),
            "single_read_spans_junction": int(single_j),
            "pair_left_overhang": int(pair_left_ov),
            "pair_right_overhang": int(pair_right_ov),
            "single_left_overhang": int(single_left_ov),
            "single_right_overhang": int(single_right_ov),
            "boundary_filter_applied": int(apply_boundary),
            "junction_left": jl, "junction_right": jr,
        })

    union = pair_junc | single_junc
    return {
        "rows": rows,
        "counts": {
            "templates_fetched": len(by_name),
            "templates_passed_filters": passed,
            "pair_in_genes": len(pair_any),
            "pair_spans_junction": len(pair_junc),
            "single_read_spans_junction": len(single_junc),
            "union_junction_support": len(union),
        },
        "rejections": dict(reject),
        "junction": {"left": jl, "right": jr, "gap_bp": gap_bp},
        "boundary_filter": {
            "min_overhang": int(min_overhang),
            "short_gap": int(short_gap),
            "applied": apply_boundary,
            "gap_bp": gap_bp,
        },
    }


def poisson_ci(k):
    try:
        from scipy.stats import poisson
        lo, hi = poisson.interval(0.95, k)
        return [round(float(lo), 2), round(float(hi), 2)]
    except Exception:
        return [None, None]


def make_verdict(mode, depth, counts, min_depth, min_junction, max_spike):
    d1 = depth["gene1"]["median"]
    d2 = depth["gene2"]["median"]
    jd = depth["junction"]["median"]
    spike = depth["span"]["spike_ratio"]
    if mode == "SE":
        return "present_single_end" if max(d1, d2) >= min_depth else "absent"
    # An intact two-gene module requires both anchors. If either anchor is
    # below the depth threshold, the module cannot be called present, even
    # when the other anchor is well covered or a boundary artefact produces
    # a few junction-spanning templates.
    if d1 < min_depth or d2 < min_depth:
        return "absent"
    support = counts["union_junction_support"]
    if (d1 >= min_depth and d2 >= min_depth
            and support >= min_junction and jd >= min_depth
            and (spike is None or spike <= max_spike)):
        return "positive"
    if support >= 1:
        return "weak"
    return "inconclusive"


# ============================================================================
# single-sample run
# ============================================================================
def run_one(args, sample, r1, r2, bam_in, outdir, log):
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()

    # reference
    refs = read_fasta(args.reference)
    if args.ref_id:
        ref_id_full = args.ref_id
    elif len(refs) == 1:
        ref_id_full = next(iter(refs))
    else:
        raise SystemExit("Reference has multiple records; pass --ref-id.")
    if ref_id_full not in refs:
        raise SystemExit(f"{ref_id_full} not in reference")

    # annotation
    if args.genbank:
        feats = parse_genbank(args.genbank)
        ann = args.genbank
    elif args.gff:
        feats = parse_gff3(args.gff)
        ann = args.gff
    else:
        feats = []
        ann = "none"

    f1 = (parse_coords(args.gene1_coords, ref_id_full)
          if args.gene1_coords else select_feature(feats, args.gene1, ref_id_full))
    f2 = (parse_coords(args.gene2_coords, ref_id_full)
          if args.gene2_coords else select_feature(feats, args.gene2, ref_id_full))

    # target
    if bam_in:
        ref_for_bam = args.reference
        ref_id = ref_id_full
        t_f1, t_f2 = f1, f2
        target_fa = args.reference
        win_s, win_e = min(f1.start, f2.start), max(f1.end, f2.end)
        index_ref(ref_for_bam, log)
    else:
        target_fa, ref_id, t_f1, t_f2, win_s, win_e = build_target(
            refs, ref_id_full, f1, f2, outdir, args.flank)
        index_ref(target_fa, log)
        ref_for_bam = target_fa

    # mapping
    bam = bam_in or os.path.join(outdir, f"{sample}.sorted.bam")
    if not bam_in:
        if not r1:
            raise SystemExit("Need --reads or --bam")
        map_reads(ref_for_bam, r1, r2, bam, sample, args.threads, log,
                  force=args.force)
    if not os.path.exists(bam + ".bai"):
        run_cmd(["samtools", "index", bam], log)
    fs = flagstat(bam)
    if not bam_in:
        # BAM is produced with unmapped/secondary/supplementary reads removed,
        # so flagstat total == mapped and a "mapping rate" would be misleading.
        fs["mapping_rate"] = None
        fs["mapping_note"] = ("filtered BAM: unmapped/secondary/supplementary "
                              "removed; input read total not counted")

    sam = pysam.AlignmentFile(bam, "rb")
    if ref_id not in sam.references:
        raise SystemExit(f"{ref_id} not in BAM; rebuild the BAM against the "
                         f"same reference/target.")
    bounds = insert_bounds(sam, ref_id, args.mapq, args.max_insert,
                           mad_mult=args.insert_mad)
    mode = "PE" if r2 else "SE"
    if bam_in:
        mode = "PE" if (fs.get("properly_paired") or 0) > 0 else "SE"

    d1 = depth_stats(sam, ref_id, t_f1.start, t_f1.end, args.mapq)
    d2 = depth_stats(sam, ref_id, t_f2.start, t_f2.end, args.mapq)
    jl, jr = (t_f1.end, t_f2.start) if t_f1.end <= t_f2.start else (
        (t_f1.end + t_f2.start) // 2, (t_f1.end + t_f2.start) // 2)
    dj = depth_stats(sam, ref_id, max(1, jl - 100), jr + 100, args.mapq)
    span_s = min(t_f1.start, t_f2.start)
    span_e = max(t_f1.end, t_f2.end)
    ds = depth_stats(sam, ref_id, span_s, span_e, args.mapq)

    res = analyze_linkage(sam, ref_id, t_f1, t_f2, args.mapq, bounds,
                          args.max_insert, min_overhang=args.min_overhang,
                          short_gap=args.short_gap)
    if (args.min_overhang == 0
            and res["boundary_filter"]["gap_bp"] <= args.short_gap):
        log("WARNING: short intergenic gap (<=%d bp) with --min-overhang 0; "
            "boundary-proximal reads can create false junction support"
            % args.short_gap)
    counts = res["counts"]
    depth = {"gene1": d1, "gene2": d2, "junction": dj, "span": ds}
    verdict = make_verdict(mode, depth, counts, args.min_depth,
                           args.min_junction, args.max_spike)

    # optional negative control on the full reference
    neg = {"pair": None, "support": None}
    if args.negative_gene and ref_id_full in refs:
        try:
            fc = select_feature(feats, args.negative_gene, ref_id_full)
            if abs(fc.start - f1.start) > 10000:
                if bam_in:
                    c_f1, c_f2 = f1, fc
                    c_ref = ref_id_full
                else:
                    # target window does not contain the control; skip
                    raise ValueError("negative control outside target window")
                rr = analyze_linkage(sam, c_ref, c_f1, c_f2, args.mapq,
                                     bounds, args.max_insert,
                                     min_overhang=args.min_overhang,
                                     short_gap=args.short_gap)
                neg = {"pair": f"{args.gene1}/{args.negative_gene}",
                       "support": rr["counts"]["union_junction_support"]}
        except Exception:
            pass

    sam.close()

    ci = poisson_ci(counts["union_junction_support"])
    summary = {
        "version": VERSION,
        "sample": sample,
        "mode": mode,
        "reference": args.reference,
        "reference_id": ref_id_full,
        "annotation": ann,
        "gene1": {"name": args.gene1, "coords_1based": [f1.start, f1.end],
                  "strand": f1.strand, "product": f1.product},
        "gene2": {"name": args.gene2, "coords_1based": [f2.start, f2.end],
                  "strand": f2.strand, "product": f2.product},
        "target_window_1based": [win_s, win_e] if not bam_in else None,
        "bam": bam,
        "mapping": fs,
        "insert_size": bounds,
        "depth": depth,
        "linkage": counts,
        "junction": res["junction"],
        "rejections": res["rejections"],
        "junction_support_ci95": ci,
        "negative_control": neg,
        "params": {
            "mapq": args.mapq, "min_depth": args.min_depth,
            "min_junction": args.min_junction, "max_spike": args.max_spike,
            "flank": args.flank, "max_insert": args.max_insert,
            "insert_mad": args.insert_mad, "threads": args.threads,
            "min_overhang": args.min_overhang,
            "short_gap": args.short_gap,
        },
        "verdict": verdict,
        "runtime_s": round(time.time() - t0, 1),
    }

    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)

    cols = ["metric", "value"]
    with open(os.path.join(outdir, "summary.tsv"), "w", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        w.writerow(["sample", sample])
        w.writerow(["verdict", verdict])
        w.writerow(["junction_support", counts["union_junction_support"]])
        w.writerow(["junction_support_ci95", ci])
        w.writerow(["gene1_depth_median", d1["median"]])
        w.writerow(["gene2_depth_median", d2["median"]])
        w.writerow(["junction_depth_median", dj["median"]])
        w.writerow(["spike_ratio", ds["spike_ratio"]])
        w.writerow(["mapping_rate", fs.get("mapping_rate")])
        w.writerow(["negative_control_support", neg["support"]])

    with open(os.path.join(outdir, "linkage_pairs.tsv"), "w",
              encoding="utf-8") as fh:
        if res["rows"]:
            w = csv.DictWriter(fh, fieldnames=list(res["rows"][0].keys()),
                               delimiter="\t")
            w.writeheader()
            w.writerows(res["rows"])
        else:
            fh.write("template\tjunction_left\tjunction_right\n")

    with open(os.path.join(outdir, "depth_regions.tsv"), "w",
              encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region", "start", "end", "n", "median", "mean",
                    "max", "p95", "spike_ratio"])
        for name, d, s, e in (
                ("gene1", d1, t_f1.start, t_f1.end),
                ("gene2", d2, t_f2.start, t_f2.end),
                ("junction", dj, max(1, jl - 100), jr + 100),
                ("span", ds, span_s, span_e)):
            w.writerow([name, s, e, d["n"], d["median"], d["mean"],
                        d["max"], d["p95"], d["spike_ratio"]])

    # report
    lines = [
        f"# RepBridge v{VERSION} — {sample}\n",
        f"- Verdict: **{verdict}**",
        f"- Mode: {mode}",
        f"- Reference: `{args.reference}` ({ref_id_full})",
        f"- Annotation: `{ann}`",
        f"- Gene1: {f1.label()}",
        f"- Gene2: {f2.label()}",
        f"- BAM: `{bam}`",
        "",
        "## Mapping\n",
        f"- reads in BAM: {fs.get('mapped')}",
        (f"- mapping rate: {100*fs['mapping_rate']:.2f}%"
         if fs.get("mapping_rate") is not None
         else "- mapping rate: not computed (filtered BAM)"),
        f"- properly paired: {fs.get('properly_paired')}",
        f"- insert median/MAD: {bounds.get('median')} / {bounds.get('mad')}",
        "",
        "## Depth\n",
        "| region | median | mean | max | p95 | spike |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, d in (("gene1", d1), ("gene2", d2),
                    ("junction", dj), ("span", ds)):
        lines.append(f"| {name} | {d['median']} | {d['mean']} | {d['max']} | "
                     f"{d['p95']} | {d['spike_ratio']} |")
    lines += [
        "",
        "## Linkage\n",
        f"- templates fetched: {counts['templates_fetched']}",
        f"- passed filters: {counts['templates_passed_filters']}",
        f"- pair in genes: {counts['pair_in_genes']}",
        f"- pair spans junction: {counts['pair_spans_junction']}",
        f"- single read spans junction: {counts['single_read_spans_junction']}",
        f"- union junction support: **{counts['union_junction_support']}** "
        f"(95% CI {ci})",
        f"- negative control: {neg}",
        "",
        "## Rejections\n",
        f"`{res['rejections']}`",
        "",
        "## Interpretation\n",
        "- `positive`: both genes covered and ≥ threshold junction templates.",
        "- `weak`: 1–threshold-1 junction templates.",
        "- `absent`: both gene depths below min-depth; no evidence in reads.",
        "- `inconclusive`: coverage present but no junction support.",
        "- `present_single_end`: single-end data; linkage cannot be assessed.",
        "",
        "> This is a targeted, reference-guided read-level assay. It does not "
        "require an assembly.",
    ]
    with open(os.path.join(outdir, "report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    log(f"  >>> {sample}: {verdict} "
        f"(support={counts['union_junction_support']}, "
        f"depth1={d1['median']}, depth2={d2['median']})")
    return summary


# ============================================================================
# CLI
# ============================================================================
def add_common(p):
    p.add_argument("--reference", default=DEFAULT_REF)
    p.add_argument("--ref-id", default=None)
    p.add_argument("--genbank", default=DEFAULT_GBK)
    p.add_argument("--gff", default=None)
    p.add_argument("--gene1", default=DEFAULT_GENE1)
    p.add_argument("--gene2", default=DEFAULT_GENE2)
    p.add_argument("--gene1-coords", default=None)
    p.add_argument("--gene2-coords", default=None)
    p.add_argument("--negative-gene", default=None)
    p.add_argument("--flank", type=int, default=100)
    p.add_argument("--mapq", type=int, default=20)
    p.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 4))
    p.add_argument("--max-insert", type=int, default=None)
    p.add_argument("--insert-mad", type=float, default=None,
                   help="optional insert-size half-window in MAD units; "
                        "default keeps the historical 4*sigma (about 6 MAD) rule")
    p.add_argument("--min-depth", type=float, default=5.0)
    p.add_argument("--min-junction", type=int, default=10)
    p.add_argument("--max-spike", type=float, default=5.0)
    p.add_argument("--min-overhang", type=int, default=10,
                   help="minimum aligned reference bases on each side of a "
                        "short intergenic junction for single-read spans; "
                        "proper pairs spanning the junction are retained; "
                        "0 disables the filter")
    p.add_argument("--short-gap", type=int, default=10,
                   help="apply --min-overhang only when the intergenic gap "
                        "is at most this many bp")
    p.add_argument("--force", action="store_true")


def main():
    ap = argparse.ArgumentParser(
        description=f"RepBridge v{VERSION}: assembly-independent ARG module integrity")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="single sample")
    p_run.add_argument("--sample", required=True)
    p_run.add_argument("--reads", nargs="+", default=None,
                       help="R1 [R2] FASTQ(.gz)")
    p_run.add_argument("--bam", default=None)
    p_run.add_argument("--out", required=True)
    add_common(p_run)

    p_panel = sub.add_parser("panel", help="batch sample sheet")
    p_panel.add_argument("--samples", required=True,
                         help="TSV: sample, r1, r2, [bam]")
    p_panel.add_argument("--out", required=True)
    add_common(p_panel)

    args = ap.parse_args()
    if args.cmd == "run":
        if not args.reads and not args.bam:
            ap.error("run needs --reads or --bam")
        os.makedirs(args.out, exist_ok=True)
        log = Logger(os.path.join(args.out, "run.log"))
        log(f"RepBridge v{VERSION}")
        require_tools(mapping=not args.bam)
        r1 = args.reads[0] if args.reads else None
        r2 = args.reads[1] if args.reads and len(args.reads) > 1 else None
        try:
            run_one(args, args.sample, r1, r2, args.bam, args.out, log)
        finally:
            log.close()
        return

    # panel
    os.makedirs(args.out, exist_ok=True)
    log = Logger(os.path.join(args.out, "panel.log"))
    log(f"RepBridge v{VERSION} panel")
    require_tools(mapping=True)
    with open(args.samples, encoding="utf-8", errors="replace") as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    if not samples:
        raise SystemExit("empty sample sheet")
    rows = []
    for s in samples:
        name = s["sample"]
        r1 = (s.get("r1") or "").strip() or None
        r2 = (s.get("r2") or "").strip() or None
        bam = (s.get("bam") or "").strip() or None
        subdir = os.path.join(args.out, name)
        try:
            summ = run_one(args, name, r1, r2, bam, subdir, log)
            rows.append(summ)
        except Exception as e:
            log(f"  ❌ {name}: {e}")
            rows.append({"sample": name, "verdict": "failed", "error": str(e)})
    with open(os.path.join(args.out, "panel_summary.tsv"), "w",
              encoding="utf-8") as fh:
        cols = ["sample", "verdict", "junction_support", "gene1_depth_median",
                "gene2_depth_median", "junction_depth_median", "spike_ratio",
                "mapping_rate", "negative_control"]
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            if r.get("verdict") == "failed":
                w.writerow([r.get("sample"), "failed", "", "", "", "", "", "", r.get("error")])
                continue
            d = r["depth"]
            w.writerow([r["sample"], r["verdict"],
                        r["linkage"]["union_junction_support"],
                        d["gene1"]["median"], d["gene2"]["median"],
                        d["junction"]["median"], d["span"]["spike_ratio"],
                        r["mapping"].get("mapping_rate"),
                        r["negative_control"].get("support")])
    counts = Counter(r.get("verdict", "") for r in rows)
    lines = [f"# RepBridge panel\n", f"- samples: {len(rows)}",
             f"- verdicts: {dict(counts)}", "",
             "| sample | verdict | junction support | gene1 depth | gene2 depth |",
             "|---|---|---:|---:|---:|"]
    for r in rows:
        if r.get("verdict") == "failed":
            lines.append(f"| {r['sample']} | failed | | | |")
        else:
            d = r["depth"]
            lines.append(f"| {r['sample']} | **{r['verdict']}** | "
                         f"{r['linkage']['union_junction_support']} | "
                         f"{d['gene1']['median']} | {d['gene2']['median']} |")
    with open(os.path.join(args.out, "panel_report.md"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"✅ panel done: {os.path.join(args.out, 'panel_report.md')}")
    log.close()


if __name__ == "__main__":
    main()
