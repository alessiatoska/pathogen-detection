# pathogen-detection

Estimate which taxa of a microbe are present in a sequencing sample

## Workflow

1. Pick the best assembly (highest contig N50) per tax ID for the microbe, plus one host genome, from NCBI, or uses assemblies you pin by accession.
2. Download them, builds one gzipped reference and indexes it for the aligner.
3. Maps the reads with `bwa mem` (single- or paired-end), or with `minimap2` for Nanopore reads (`--ont`), filtering alignments as they stream instead of writing a SAM.
4. Counts each read for the genome it fits best (MAPQ ≥ 30, adjustable). Reads that fit several genomes equally well are shared between them in proportion to each genome's own reads (EM).
5. Report alignments and relative abundance for every genome

## Install

```bash
conda env create -f environment.yml
conda activate pathogen-detection
```

## Usage

```bash
python pathogen_detection.py -m "Xylella fastidiosa" -H "Olea europaea" \
    -s sample_R1.fastq.gz sample_R2.fastq.gz -t 8 -o results

# pin the host to a specific assembly (cultivated olive 'Farga')
python pathogen_detection.py -m "Xylella fastidiosa" --host-accession GCA_902713445.1 \
    -s sample_R1.fastq.gz sample_R2.fastq.gz -t 8 -o results

# Nanopore reads
python pathogen_detection.py -m "Xylella fastidiosa" -H "Olea europaea" --ont \
    -s nanopore_reads.fastq.gz -t 8 -o results
```

| Flag | |
|---|---|
| `-m, --microbe` / `-H, --host` | Looks up the microbe / host on NCBI by species name |
| `--microbe-accession` / `--host-accession` | Adds a specific NCBI assembly (repeatable); without `-m` / `-H`, only these are used |
| `-s, --sample` | Takes your reads: one file (single-end) or two (paired-end); one with `--ont` |
| `--ont` | Maps Oxford Nanopore reads with minimap2 instead of bwa |
| `-t, --threads` | Sets the aligner's threads (default 1) |
| `-o, --outdir` | Sets the output directory (default `pathogen_detection_out`) |
| `--min-mapq` | Sets the minimum mapping quality for a read to count (default 30) |
| `--max-microbe-genomes` / `--max-host-genomes` | Caps how many genomes the name search takes, best first; 0 = all (defaults 0 / 1) |
| `--keep-sam` | Also saves the full alignment file (SAM) |
| `-v, --verbose` | Prints debug logs, including every command run |

You always need a microbe and a host; give each by name, by accession, or both.

## Output

| File | |
|---|---|
| `abundance.tsv`, `abundance.json` | One row per genome, most alignments first |
| `mapping_stats.json` | How many reads were placed on one genome, shared between genomes, unmapped, or dropped (also printed at the end of the run) |
| `reference.fasta.gz` and its index | The combined reference, indexed for bwa and/or minimap2 (`.mmi`); `reference.fasta.gz.json` lists the genomes in it |
| `genomes/` | Downloaded NCBI zips and summary JSONs |
| `mapping.sam` | Every alignment; only with `--keep-sam` |

Abundance columns: `tax_id`, `name`, `strain`, `role` (microbe or host), `accession`, `alignments` (each mate of a pair counts once), `shared` (how many of those came from reads shared with other genomes), `pct_of_assigned` (host included) and `pct_of_microbe` (blank for the host).

Re-running in the same output directory reuses the downloads, and the index too if the genomes are unchanged. Mapping needs only the index, so `genomes/` and `reference.fasta.gz` can be deleted to save space.

## Acknowledgments

This project originated as coursework for Bioinformatics in Plant Pathology (2024) at Sapienza University of Rome, and extended here.
