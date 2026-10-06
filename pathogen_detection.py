#!/usr/bin/env python3
"""pathogen_detection: estimate which taxa of a microbe are present in a sequencing sample.

Pipeline
    1. Pick the best assembly (highest contig N50) per tax ID for the microbe, plus a host
       genome, from NCBI with the `datasets` CLI, and/or use assemblies pinned by accession.
    2. Download them, concatenate them into one gzipped reference and index it for the
       aligner. Indexes are reused on later runs for as long as the genomes stay the same.
    3. Map the sample reads with `bwa mem` (single-end or paired-end Illumina) or, with
       --ont, `minimap2 -x map-ont` (Oxford Nanopore), filtering the alignment stream as it
       is produced rather than storing the whole SAM.
    4. Count each read for the genome it was placed on: confident primary alignments (MAPQ
       >= 30 by default), and reads whose equally good hits all lie in one genome. Reads that
       fit several genomes equally well are shared between them by EM, in proportion to
       each genome's abundance.
    5. Write a TSV/JSON report with the relative abundance of every genome, plus a breakdown
       of what happened to each read.

The host genome is included as a decoy so host reads do not get mis-assigned to the microbe.
Only the standard library is required; `datasets` and `bwa` (or `minimap2` for --ont) must
be on PATH.
"""

import argparse
import contextlib
import gzip
import json
import logging
import shutil
import subprocess
import sys
import time
import zipfile
from collections import Counter
from pathlib import Path

log = logging.getLogger("pathogen_detection")

ASSEMBLY_LEVELS = ["complete", "chromosome", "scaffold", "contig"]

# NCBI infraspecific name fields, in order of preference, for the report's strain column
STRAIN_FIELDS = ["strain", "isolate", "cultivar", "ecotype", "breed"]

# files each aligner's index adds next to the reference; mapping needs only these
INDEX_SUFFIXES = {
    "bwa": [".amb", ".ann", ".bwt", ".pac", ".sa"],
    "minimap2": [".mmi"],
}

# SAM flag bits: 0x4 unmapped, 0x100 secondary, 0x800 supplementary
UNMAPPED = 0x4
NOT_PRIMARY = 0x100 | 0x800

# bwa lists a read's alternative hits (XA tag) only if it has fewer than this many; reads
# with more, such as host repeats, cannot be shared out and are dropped below --min-mapq
MAX_LISTED_HITS = 50

# reads seen while mapping, and the outcome for each (every read gets exactly one)
STAT_KEYS = ["reads", "unmapped", "low_mapq", "unknown_reference", "shared", "assigned"]

REPORT_COLUMNS = ["tax_id", "name", "strain", "role", "accession",
                  "alignments", "shared", "pct_of_assigned", "pct_of_microbe"]


# --------------------------------------------------------------------------- helpers
def run(cmd, **kwargs):
    """Run a command, raising a readable error on failure.

    Captured stderr is included in the error. Stderr that is not captured goes straight to
    the terminal, which is how long downloads and indexing show their progress.
    """
    log.debug("running: %s", " ".join(map(str, cmd)))
    try:
        return subprocess.run(cmd, check=True, text=True, **kwargs)
    except subprocess.CalledProcessError as exc:
        stderr = (exc.stderr or "").strip()
        raise RuntimeError(
            f"command failed ({exc.returncode}): {' '.join(map(str, cmd))}\n{stderr}") from exc


def check_tools(tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise SystemExit(f"Missing required tools on PATH: {', '.join(missing)}")


def elapsed(start):
    """Time since `start` (a time.monotonic() reading), e.g. '4m 07s'."""
    seconds = round(time.monotonic() - start)
    return f"{seconds // 60}m {seconds % 60:02d}s"


# --------------------------------------------------------------------------- NCBI download
def download_genome_summary(organism, outdir):
    """Return the summary JSON (dict) for the best assembly level that has data, else None."""
    for level in ASSEMBLY_LEVELS:
        log.info("Looking for %s-level genomes of '%s'", level, organism)
        try:
            result = run(
                ["datasets", "summary", "genome", "taxon", organism, "--assembly-level", level],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            data = json.loads(result.stdout)
        except (RuntimeError, json.JSONDecodeError) as exc:
            log.warning("Summary failed for %s at %s level: %s", organism, level, exc)
            continue
        if data.get("reports"):
            path = Path(outdir) / f"{organism.replace(' ', '_')}_{level}.json"
            path.write_text(json.dumps(data, indent=2))
            log.info("Saved summary to %s", path)
            return data
        log.info("No %s-level genomes found for '%s'", level, organism)
    return None


def assembly_info(report):
    """The fields used from one NCBI assembly report, or None without a tax ID or accession."""
    organism = report.get("organism", {})
    tax_id = organism.get("tax_id")
    accession = report.get("accession")
    if tax_id is None or not accession:
        return None
    names = organism.get("infraspecific_names") or {}
    stats = report.get("assembly_stats", {})
    return {
        "tax_id": tax_id,
        "name": organism.get("organism_name", str(tax_id)),
        "strain": next((names[f] for f in STRAIN_FIELDS if names.get(f)), ""),
        "accession": accession,
        "contig_n50": stats.get("contig_n50", 0) or 0,
        "length": int(stats.get("total_sequence_length", 0) or 0),  # a string in NCBI's JSON
    }


def top_strain_per_taxon(summary):
    """Pick the assembly with the highest contig N50 for each tax ID.

    Returns a list of assembly_info() dicts, best assembly first, so that capping the number
    of genomes keeps the best ones rather than arbitrary ones.
    """
    best = {}
    for report in summary.get("reports", []):
        info = assembly_info(report)
        if info is None:
            continue
        tax_id = info["tax_id"]
        if tax_id not in best or info["contig_n50"] > best[tax_id]["contig_n50"]:
            best[tax_id] = info
    return sorted(best.values(), key=lambda g: g["contig_n50"], reverse=True)


def summarize_accessions(accessions):
    """Return assembly_info() dicts for specific accessions, in the order given.

    An accession without a version (GCF_041053125) matches the version NCBI returns. NCBI
    answers an unknown accession with an empty result rather than an error, so missing
    accessions are checked for here.
    """
    result = run(["datasets", "summary", "genome", "accession", *accessions],
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    found = {}
    for report in json.loads(result.stdout).get("reports", []):
        info = assembly_info(report)
        if info is not None:
            found[info["accession"]] = info
            found.setdefault(info["accession"].split(".")[0], info)
    missing = [a for a in accessions if a not in found]
    if missing:
        raise SystemExit(f"Assembly accession(s) not found at NCBI: {', '.join(missing)}")
    return [found[a] for a in accessions]


def describe_assembly(info):
    strain = f" [{info['strain']}]" if info["strain"] else ""
    return f"'{info['name']}': {info['accession']}{strain} (contig N50 {info['contig_n50']})"


def select_genomes(organism, role, workdir, max_genomes=0, accessions=()):
    """Choose the genomes for one role ("microbe" or "host"), without downloading them.

    With an organism name, takes the best assembly per tax ID, capped at max_genomes
    (0 = all). Pinned accessions are added to those, or used alone when there is no name.
    Returns assembly_info() dicts with the role added.
    """
    candidates = []
    if organism:
        summary = download_genome_summary(organism, workdir)
        if summary is None:
            raise SystemExit(f"No genomes found for '{organism}'")
        candidates = top_strain_per_taxon(summary)
        if max_genomes > 0:
            candidates = candidates[:max_genomes]
        for info in candidates:
            log.info("Best assembly for %s", describe_assembly(info))
    if accessions:
        for info in summarize_accessions(accessions):
            if info["accession"] not in {c["accession"] for c in candidates}:
                log.info("Pinned assembly for %s", describe_assembly(info))
                candidates.append(info)
    return [{**info, "role": role} for info in candidates]


def fasta_member(zf):
    """Name of the genome FASTA inside an NCBI datasets zip."""
    names = sorted(n for n in zf.namelist() if n.endswith(".fna"))
    if not names:
        raise RuntimeError(f"No .fna file in {zf.filename}")
    return names[0]


def genome_zip_ok(path):
    """True if path is a complete datasets zip with a genome FASTA in it."""
    try:
        with zipfile.ZipFile(path) as zf:
            fasta_member(zf)
    except (OSError, zipfile.BadZipFile, RuntimeError):
        return False
    return True


def download_genome(accession, workdir):
    """Download a genome by accession and return the path of its zip.

    The zip is kept compressed and the FASTA is read straight out of it, so a genome takes
    its compressed size on disk once, rather than its full size twice (extracted, then
    copied into the reference). A zip left by an earlier run is reused if it is complete.
    """
    zip_path = Path(workdir) / f"{accession}.zip"
    if genome_zip_ok(zip_path):
        return zip_path
    if zip_path.exists():
        log.warning("Download of %s from an earlier run is incomplete; fetching it again",
                    accession)
        zip_path.unlink()
    log.info("Downloading %s", accession)
    start = time.monotonic()
    # datasets reports progress on stderr, which is left on the terminal
    run(["datasets", "download", "genome", "accession", accession,
         "--filename", str(zip_path)], stdout=subprocess.PIPE)
    if not genome_zip_ok(zip_path):
        raise RuntimeError(f"Corrupt download for {accession}: {zip_path}")
    log.info("Downloaded %s (%.0f MB) in %s",
             accession, zip_path.stat().st_size / 1e6, elapsed(start))
    return zip_path


# --------------------------------------------------------------------------- reference
@contextlib.contextmanager
def open_fasta(zip_path):
    """Open the genome FASTA inside a datasets zip for binary reading."""
    with zipfile.ZipFile(zip_path) as zf, zf.open(fasta_member(zf)) as handle:
        yield handle


def build_reference(genomes, reference_path):
    """Concatenate the genomes' FASTAs into one gzipped FASTA; return {sequence_id: genome}.

    Lines are streamed verbatim in binary, which keeps large host genomes fast to
    concatenate while still recording which genome each sequence id came from.
    """
    seq_to_genome = {}
    # level 1: ~10x faster than gzip's default on DNA, for ~13% more space
    with gzip.open(reference_path, "wb", compresslevel=1) as out:
        for genome in genomes:
            last = b"\n"
            with open_fasta(genome["zip"]) as handle:
                for line in handle:
                    if line.startswith(b">"):
                        seq_id = line[1:].split()[0].decode("utf-8", "replace")
                        if seq_id in seq_to_genome:
                            raise RuntimeError(
                                f"Duplicate sequence id {seq_id} across genomes")
                        seq_to_genome[seq_id] = genome
                    out.write(line)
                    last = line
            if not last.endswith(b"\n"):  # a missing trailing newline would fuse two records
                out.write(b"\n")
    return seq_to_genome


def build_index(reference, aligner, length=0):
    """Index the reference for bwa (short reads) or minimap2 (Nanopore reads).

    bwa chooses its own algorithm: `is` for small genomes and `bwtsw` for large ones.
    Forcing `is` on a plant genome fails or runs out of memory, since it needs about 5x the
    genome size in RAM.
    """
    size = f"{length / 1e6:,.1f} Mbp" if length else "size unknown"
    log.info("Indexing the reference (%s) for %s; for a plant genome this is the slowest step",
             size, aligner)
    start = time.monotonic()
    if aligner == "minimap2":
        cmd = ["minimap2", "-x", "map-ont", "-d", f"{reference}.mmi", str(reference)]
    else:
        cmd = ["bwa", "index", str(reference)]
    # progress goes to stderr, which is left on the terminal
    run(cmd, stdout=subprocess.PIPE)
    log.info("Indexed in %s", elapsed(start))


def load_manifest(reference, genomes):
    """Return the reference's manifest if it was built from exactly these genomes, else None.

    The manifest maps each sequence id to its genome's accession and lists the aligners
    whose index is complete. Mapping needs only the manifest and the index: aligners never
    read the FASTA, so the reference and the downloads may have been deleted to save space.
    """
    try:
        manifest = json.loads(Path(f"{reference}.json").read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if (sorted(manifest.get("accessions", [])) != sorted(g["accession"] for g in genomes)
            or "sequences" not in manifest):
        return None
    return manifest


def index_ready(reference, aligner, manifest):
    """True if the aligner's index of this reference finished building and is still there."""
    return (aligner in manifest.get("indexes", [])
            and all(Path(f"{reference}{s}").is_file() for s in INDEX_SUFFIXES[aligner]))


def write_reference(genomes, reference, workdir):
    """Download the genomes and concatenate them into the reference; return its manifest."""
    genomes = [{**g, "zip": download_genome(g["accession"], workdir)} for g in genomes]
    start = time.monotonic()
    seq_to_genome = build_reference(genomes, reference)
    log.info("Wrote %s (%d sequences) in %s", reference, len(seq_to_genome), elapsed(start))
    return {
        "accessions": [g["accession"] for g in genomes],
        "sequences": {seq_id: g["accession"] for seq_id, g in seq_to_genome.items()},
        "indexes": [],
    }


def prepare_reference(genomes, reference, workdir, aligner="bwa"):
    """Make sure the aligner has an index of exactly these genomes; return {seq_id: genome}.

    The reference is rebuilt only when the genomes change, and each aligner's index is built
    once per reference, so switching between bwa and minimap2 reuses what is there. The
    manifest is updated only after each step finishes, so an interrupted step is redone on
    the next run instead of being reused half-written.
    """
    manifest_path = Path(f"{reference}.json")
    manifest = load_manifest(reference, genomes)
    if manifest is None:
        manifest_path.unlink(missing_ok=True)
        for suffixes in INDEX_SUFFIXES.values():  # indexes of other genomes are stale
            for suffix in suffixes:
                Path(f"{reference}{suffix}").unlink(missing_ok=True)
        manifest = write_reference(genomes, reference, workdir)
        manifest_path.write_text(json.dumps(manifest))
    elif index_ready(reference, aligner, manifest):
        log.info("Reusing the %s index of %s, which already holds these %d genomes",
                 aligner, reference, len(genomes))
    elif not Path(reference).exists():
        write_reference(genomes, reference, workdir)  # deleted to save space; same genomes
    if not index_ready(reference, aligner, manifest):
        build_index(reference, aligner, sum(g["length"] for g in genomes))
        manifest["indexes"] = sorted({*manifest.get("indexes", []), aligner})
        manifest_path.write_text(json.dumps(manifest))
    by_accession = {g["accession"]: g for g in genomes}
    return {seq_id: by_accession[acc] for seq_id, acc in manifest["sequences"].items()}


# --------------------------------------------------------------------------- SAM parsing
def classify_sam_record(line, min_mapq=30):
    """Return (outcome, reference names) for one SAM line.

    Headers, malformed lines and secondary or supplementary records (extra alignments of a
    read, not reads) give (None, None), and an unmapped read ("unmapped", None). A
    confident read (MAPQ >= min_mapq) gives ("pass", [its reference]); any other gives
    ("low_mapq", every reference it fits as well as its primary one), or ("low_mapq", None)
    when the aligner did not list its other hits. A confident read counts even if it also
    has alternative hits, since MAPQ already accounts for them.
    """
    if not line or line.startswith("@"):
        return None, None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return None, None
    try:
        flag, mapq = int(fields[1]), int(fields[4])
    except ValueError:
        return None, None
    if flag & NOT_PRIMARY:
        return None, None
    ref = fields[2]
    if flag & UNMAPPED or ref == "*":
        return "unmapped", None
    if mapq >= min_mapq:
        return "pass", [ref]
    return "low_mapq", equally_good_hits(ref, fields[11:])


def equally_good_hits(ref, tags):
    """The primary reference plus each alternative hit (bwa's XA tag) with no more edits.

    Returns None without an XA and NM tag, as the read's other hits are then unknown.
    """
    tags = {tag[:2]: tag[5:] for tag in tags}
    if "XA" not in tags or "NM" not in tags:
        return None
    edits = int(tags["NM"])
    refs = [ref]
    for hit in tags["XA"].rstrip(";").split(";"):  # each hit is "ref,pos,CIGAR,NM"
        hit_ref, *_, hit_edits = hit.split(",")
        if int(hit_edits) <= edits:
            refs.append(hit_ref)
    return refs


def count_lines(lines, seq_to_genome, min_mapq=30, sam_out=None):
    """Count the reads placed on each genome over an iterable of SAM lines.

    Returns (counts, ambiguous, stats). counts has the reads placed on a single genome,
    keyed by accession: confident reads, plus low-MAPQ reads whose equally good hits are all
    in one genome (a repeat within it). ambiguous has the reads that fit several genomes
    equally well, keyed by the sorted tuple of their accessions, for share_ambiguous().
    stats counts the reads (primary records, so each mate of a pair is one read) and the
    outcome for each, as in STAT_KEYS.
    """
    counts, ambiguous = Counter(), Counter()
    stats = dict.fromkeys(STAT_KEYS, 0)
    for line in lines:
        if sam_out is not None:
            sam_out.write(line)
        outcome, refs = classify_sam_record(line, min_mapq)
        if outcome is None:
            continue
        stats["reads"] += 1
        if refs is None:  # unmapped, or low MAPQ with no list of its other hits
            stats[outcome] += 1
            continue
        genomes = [seq_to_genome.get(ref) for ref in refs]
        if None in genomes:
            stats["unknown_reference"] += 1
            continue
        accessions = tuple(sorted({g["accession"] for g in genomes}))
        if len(accessions) == 1:
            counts[accessions[0]] += 1
            stats["assigned"] += 1
        else:
            ambiguous[accessions] += 1
            stats["shared"] += 1
    return counts, ambiguous, stats


def mapping_command(aligner, reference, reads, threads):
    """The command that maps the reads to the indexed reference, writing SAM to stdout."""
    if aligner == "minimap2":
        # secondary alignments are never counted, so minimap2 need not report them
        return ["minimap2", "-a", "-x", "map-ont", "--secondary=no", "-t", str(threads),
                f"{reference}.mmi", *map(str, reads)]
    # -h makes bwa list up to MAX_LISTED_HITS - 1 alternative hits, so reads can be shared
    return ["bwa", "mem", "-h", str(MAX_LISTED_HITS), "-t", str(threads), str(reference),
            *map(str, reads)]


def map_and_count(reference, reads, threads, seq_to_genome, min_mapq=30, sam_path=None, aligner="bwa"):
    """Run the aligner and count its output as it streams; return count_lines()'s result.

    The full SAM is only written to disk when sam_path is given, because a complete SAM
    for a real read set is often tens of gigabytes.
    """
    cmd = mapping_command(aligner, reference, reads, threads)
    log.info("Mapping %d read file(s) with %s, %s thread(s)", len(reads), aligner, threads)
    log.debug("running: %s", " ".join(cmd))
    # progress goes to stderr and is left on the terminal so long runs stay visible.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
    sam_out = open(sam_path, "w") if sam_path else None
    try:
        result = count_lines(proc.stdout, seq_to_genome, min_mapq, sam_out)
    finally:
        if sam_out:
            sam_out.close()
        proc.stdout.close()
        returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(f"{aligner} failed with exit code {returncode}")
    return result


# --------------------------------------------------------------------------- abundance
def share_ambiguous(counts, ambiguous, max_rounds=10000):
    """Split the reads that fit several genomes equally well between them (EM).

    Each such read is shared in proportion to the genomes' abundance, starting from equal
    abundances; the abundances are then re-estimated from the result and the reads shared
    again, until nothing changes. A genome with more reads of its own thus gets more of the
    reads it shares, and one with none gets almost none. Returns {accession: shared reads}.
    """
    if not ambiguous:
        return {}
    accessions = set(counts) | {a for group in ambiguous for a in group}
    total = sum(counts.values()) + sum(ambiguous.values())
    abundance = dict.fromkeys(accessions, 1 / len(accessions))
    for _ in range(max_rounds):
        shared = dict.fromkeys(accessions, 0.0)
        for group, n in ambiguous.items():
            weight = sum(abundance[a] for a in group)
            for a in group:
                shared[a] += n * abundance[a] / weight
        previous = abundance
        abundance = {a: (counts.get(a, 0) + shared[a]) / total for a in accessions}
        if max(abs(abundance[a] - previous[a]) for a in accessions) < 1e-12:
            break
    return shared


# --------------------------------------------------------------------------- reporting
def build_report(counts, genomes, shared):
    """One row per genome, most reads first, with relative abundance.

    counts are the reads placed on a single genome, and shared each genome's part of the
    reads that fit several genomes (from share_ambiguous). Genomes without reads are listed
    too, so an absent taxon shows as 0 rather than as a missing row. The host is kept out
    of the microbe percentages.
    """
    reads = {g["accession"]: counts.get(g["accession"], 0) + shared.get(g["accession"], 0)
             for g in genomes}
    total_assigned = sum(reads.values())
    microbe_total = sum(reads[g["accession"]] for g in genomes if g["role"] == "microbe")
    rows = []
    for genome in sorted(genomes, key=lambda g: reads[g["accession"]], reverse=True):
        n = reads[genome["accession"]]
        is_microbe = genome["role"] == "microbe"
        rows.append({
            "tax_id": genome["tax_id"],
            "name": genome["name"],
            "strain": genome["strain"],
            "role": genome["role"],
            "accession": genome["accession"],
            "alignments": round(n),
            "shared": round(shared.get(genome["accession"], 0)),
            "pct_of_assigned": round(100 * n / total_assigned, 3) if total_assigned else 0.0,
            "pct_of_microbe": (round(100 * n / microbe_total, 3) if is_microbe and microbe_total else None),
        })
    return rows


def write_report(rows, outdir):
    tsv = Path(outdir) / "abundance.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join("" if row[c] is None else str(row[c])
                               for c in REPORT_COLUMNS) + "\n")
    (Path(outdir) / "abundance.json").write_text(json.dumps(rows, indent=2))
    return tsv


def filter_summary(stats, min_mapq, paired=False):
    """Lines saying what happened to the reads: placed, shared out, or dropped and why."""
    labels = {
        "unmapped": "unmapped",
        "low_mapq": f"MAPQ below {min_mapq} and could not be shared out, dropped",
        "unknown_reference": "on a sequence missing from the reference map",
        "shared": "fit several genomes equally well, shared between them",
        "assigned": "placed on one genome",
    }
    reads = stats["reads"]
    pairs = " (each mate of a pair counts as one)" if paired else ""
    lines = [f"{reads} reads{pairs}:"]
    for key, label in labels.items():
        if key == "unknown_reference" and not stats[key]:
            continue  # only possible when the index and its manifest disagree
        pct = 100 * stats[key] / reads if reads else 0.0
        lines.append(f"  {stats[key]:>10}  {pct:5.1f}%  {label}")
    return lines


# --------------------------------------------------------------------------- CLI
def parse_arguments(argv=None):
    p = argparse.ArgumentParser(
        description="Detect and quantify microbe taxa in a sample, using a host genome as decoy.")
    p.add_argument("-m", "--microbe", help="Microbe taxon name, e.g. 'Xylella fastidiosa'")
    p.add_argument("-H", "--host", help="Host taxon name, e.g. 'Olea europaea'")
    p.add_argument("--microbe-accession", action="append", default=[], metavar="ACC",
                   help="Microbe assembly to include, e.g. GCF_041053125.1 (repeatable). "
                        "Added to what -m finds, or used alone without -m")
    p.add_argument("--host-accession", action="append", default=[], metavar="ACC",
                   help="Host assembly to include, e.g. GCA_902713445.1 (repeatable). "
                        "Added to what -H finds, or used alone without -H")
    p.add_argument("-s", "--sample", required=True, nargs="+",
                   help="1 (single-end) or 2 (paired-end) read files; 1 with --ont")
    p.add_argument("--ont", action="store_true",
                   help="Reads are Oxford Nanopore: map them with minimap2 instead of bwa")
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("-o", "--outdir", default="pathogen_detection_out", help="Output directory")
    p.add_argument("--min-mapq", type=int, default=30,
                   help="Minimum mapping quality (default: 30)")
    p.add_argument("--max-host-genomes", type=int, default=1,
                   help="Host genomes to download, best first (0 = all)")
    p.add_argument("--max-microbe-genomes", type=int, default=0,
                   help="Microbe genomes to download, best first (0 = all)")
    p.add_argument("--keep-sam", action="store_true",
                   help="Also write the full SAM (large) for debugging")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    if not args.microbe and not args.microbe_accession:
        p.error("give a microbe taxon (-m), --microbe-accession, or both")
    if not args.host and not args.host_accession:
        p.error("give a host taxon (-H), --host-accession, or both")
    if len(args.sample) > 2:
        p.error("--sample takes 1 or 2 files")
    if args.ont and len(args.sample) != 1:
        p.error("--ont takes 1 read file; concatenate multiple FASTQ files into one first")
    for f in args.sample:
        if not Path(f).is_file():
            p.error(f"sample file not found: {f}")
    return args


def main(argv=None):
    args = parse_arguments(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    aligner = "minimap2" if args.ont else "bwa"
    check_tools(["datasets", aligner])

    outdir = Path(args.outdir)
    workdir = outdir / "genomes"
    workdir.mkdir(parents=True, exist_ok=True)

    genomes = select_genomes(args.microbe, "microbe", workdir, args.max_microbe_genomes,
                             args.microbe_accession)
    genomes += select_genomes(args.host, "host", workdir, args.max_host_genomes,
                              args.host_accession)
    log.info("Reference has %d genome(s)", len(genomes))

    reference = outdir / "reference.fasta.gz"
    seq_to_genome = prepare_reference(genomes, reference, workdir, aligner)

    start = time.monotonic()
    counts, ambiguous, stats = map_and_count(
        reference, args.sample, args.threads, seq_to_genome, args.min_mapq,
        sam_path=outdir / "mapping.sam" if args.keep_sam else None, aligner=aligner)
    log.info("Mapped in %s", elapsed(start))
    shared = share_ambiguous(counts, ambiguous)

    rows = build_report(counts, genomes, shared)
    tsv = write_report(rows, outdir)
    (outdir / "mapping_stats.json").write_text(
        json.dumps({**stats, "min_mapq": args.min_mapq}, indent=2))

    for line in filter_summary(stats, args.min_mapq, paired=len(args.sample) == 2):
        log.info("%s", line)
    if stats["unknown_reference"]:
        log.warning("The index does not match %s.json; delete that file to rebuild it",
                    reference)
    for row in rows:
        pct = row["pct_of_microbe"] if row["pct_of_microbe"] is not None else row["pct_of_assigned"]
        strain = f" [{row['strain']}]" if row["strain"] else ""
        of_which = f" ({row['shared']} shared)" if row["shared"] else ""
        print(f"{row['name']}{strain} ({row['role']}, taxid {row['tax_id']}, "
              f"{row['accession']}): {row['alignments']} alignments{of_which}, {pct}%")
    log.info("Report written to %s", tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
