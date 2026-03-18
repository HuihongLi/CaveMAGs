# CaveMAGs

A comprehensive prokaryotic genomic catalog from 37 geographically diverse cave environments, built with an optimized metagenome-assembled genome (MAG) reconstruction pipeline.

---

## Overview

Cave microorganisms are unique extremophiles that have evolved in isolated, nutrient-limited environments and harbor exceptional metabolic capabilities. Yet knowledge of cave microbial diversity at the genomic level remains limited — previous studies have focused on individual caves and do not provide a global picture.

Here we present the first prokaryotic cave metagenomic catalog spanning **37 geographically diverse cave environments**. Using an optimized genome reconstruction pipeline applied to **241 public metagenomic datasets**, we recovered **3,837 medium-to-high quality metagenome-assembled genomes (MAGs)**. These were dereplicated into **1,979 species-level representative clusters** covering **67 phyla** of Bacteria (*n* = 1,858) and Archaea (*n* = 121).

### Key Findings

| Metric | Value |
|---|---|
| Cave environments | 37 |
| Metagenomic samples | 241 |
| MAGs recovered | 3,837 |
| Dereplicated species clusters | 1,979 |
| Phyla represented | 67 |
| Bacterial genomes | 1,858 |
| Archaeal genomes | 121 |
| Novel genomes (< 95% ANI to named species) | 98.7% |
| Genomes with biosynthetic gene clusters (BGCs) | 98.0% |
| Genomes with antibiotic resistance genes (ARGs) | 95.0% |

This catalog provides a foundational resource for exploring cave microbial diversity, secondary metabolism, and the evolutionary origins of antibiotic resistance in subterranean ecosystems.

---

## Repository Structure

```
CaveMAGs/
├── Script.sh          # Complete bioinformatics pipeline (14 steps)
├── Output/
│   ├── 1.sra_list.txt # NCBI SRA accession numbers for all 241 samples
│   └── 2.samples.tsv  # Sample metadata: SampleID, front-trim bases (F), tail-trim bases (T)
└── README.md
```

### Output Files

| File | Description |
|---|---|
| `Output/1.sra_list.txt` | 241 NCBI SRA accession numbers (prefetch input) |
| `Output/2.samples.tsv` | Three-column TSV: `SampleID`, `F`, `T`. `F` is passed to fastp as `-f`/`-F` (hard-trim *F* bases from the 5′ end of read 1 and read 2); `T` is passed as `-t`/`-T` (hard-trim *T* bases from the 3′ end of read 1 and read 2). Values were chosen per sample to remove low-quality cycles at the read termini. |

---

## Pipeline

The pipeline is documented as a series of shell commands in [`Script.sh`](Script.sh). It covers 14 sequential steps:

### Step 1 — Download & Convert SRA Data

Downloads raw sequencing data with `prefetch` using the accession list (`Output/1.sra_list.txt`) and converts each `.sra` file to paired-end FASTQ with `fasterq-dump`. Output files are compressed with `pigz` (or `gzip` as fallback).

### Step 2 — Quality Control

- **FastQC**: Per-sample quality metrics.
- **MultiQC**: Aggregate QC summary across all samples.
- **fastp**: Adapter trimming, quality filtering (Phred ≥ 20, < 25% uncalled bases), length filtering (50–300 bp), poly-G/X trimming, and optical deduplication. Per-sample hard-trim lengths are read from `Output/2.samples.tsv`: column `F` sets the number of bases removed from the 5′ end (`-f`/`-F`) and column `T` sets the number removed from the 3′ end (`-t`/`-T`) to eliminate low-quality cycles specific to each sequencing run.

### Step 3 — Metagenomic Assembly

Assembles cleaned paired-end reads with **MEGAHIT** (`--presets meta-sensitive`, minimum contig length 1,500 bp).

### Step 4 — Read Alignment

Maps sample reads back to their assembled contigs with **Bowtie2** and produces sorted, indexed BAM files using **SAMtools**. The BAM files provide coverage depth information required for binning.

### Step 5 — Genome Binning

Three independent binning strategies are applied to each assembly:

| Tool | Approach |
|---|---|
| **MetaBAT2** | Differential-abundance binning (minimum contig 2,500 bp, minimum cluster size 200 kbp) |
| **SemiBin2** | Machine-learning binning using a global pretrained model |
| **MetaDecoder** | Seed-based, coverage-guided clustering |

### Step 6 — Bin Refinement

**MAGScoT** (R script) integrates the three binning results and removes contamination using single-copy marker genes (maximum contamination 1%, score threshold 0.5).

### Step 7 — Quality Assessment

**CheckM2** estimates completeness and contamination of all refined bins. Medium-to-high quality MAGs (completeness ≥ 50%, contamination ≤ 10% per MIMAG standards) are retained for downstream analysis.

### Step 8 — Dereplication

**dRep** clusters MAGs at 90% ANI (primary clustering) and 95% ANI (secondary/species-level clustering) using the **fastANI** algorithm, yielding 1,979 representative genomes.

### Step 9 — Taxonomic Classification

**GTDB-Tk** (`classify_wf`, full tree mode) assigns taxonomy to all representative MAGs according to the Genome Taxonomy Database.

### Step 10 — Phylogenetic Analysis

**IQ-TREE** infers a maximum-likelihood phylogenetic tree from the GTDB-Tk multiple sequence alignment (MSA) with automatic model selection (`MFP`), 1,000 ultrafast bootstrap replicates (`-B 1000`), and 1,000 SH-aLRT replicates (`-alrt 1000`).

### Step 11 — Genome Coverage Analysis

**CoverM** (`genome` mode) calculates per-genome mean coverage, relative abundance, covered fraction, variance, read counts, reads per base, and genome length across all samples. Minimum covered fraction is set to 0 so all genomes are reported.

### Step 12 — Biosynthetic Gene Cluster (BGC) Detection

**antiSMASH** annotates each representative MAG for secondary metabolite biosynthetic gene clusters with:
- Gene prediction via **Prodigal**
- MIBiG database comparison (`--cc-mibig`)
- ClusterBlast comparisons against general, sub-cluster, and known-cluster databases
- Active site finder, cluster HMMer, TIGRFAM, Pfam2GO, RRE, and TFBS modules

### Step 13 — Antibiotic Resistance Gene (ARG) Annotation

**DRAMMA** identifies ARGs using a multi-feature approach:
- **Prodigal** predicts protein-coding genes (`.faa`, `.ffn`, `.gff`) from each MAG.
- `run_DRAMMA_pipeline.py` searches for ARGs via **HMMER**, **MMseqs2**, and **TMHMM** (transmembrane topology), with an E-value threshold of 1 × 10⁻⁵ and gene/nucleotide context windows of 20 genes / 20,000 bp.

### Step 14 — Functional Annotation (DRAM)

**DRAM** provides metabolic and functional annotation of all representative MAGs using KOfam, Pfam, dbCAN (CAZymes), MEROPS peptidases, and AMG/ETC module databases.

---

## Dependencies

| Tool | Purpose |
|---|---|
| [SRA Toolkit](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`) | Download and convert SRA data |
| [pigz](https://zlib.net/pigz/) | Parallel gzip compression |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-sample quality metrics |
| [MultiQC](https://multiqc.info/) | Aggregate QC reports |
| [fastp](https://github.com/OpenGene/fastp) | Adapter trimming and quality filtering |
| [MEGAHIT](https://github.com/voutcn/megahit) | Metagenomic assembly |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | Read alignment |
| [SAMtools](http://www.htslib.org/) | BAM/SAM manipulation |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/) | Genome binning |
| [SemiBin2](https://github.com/BigDataBiology/SemiBin) | Genome binning |
| [MetaDecoder](https://github.com/liu-congcong/metadecoder) | Genome binning |
| [MAGScoT](https://github.com/ikmb/MAGScoT) | Bin refinement |
| [CheckM2](https://github.com/chklovski/CheckM2) | Bin quality assessment |
| [dRep](https://github.com/MrOlm/drep) | Genome dereplication |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | Taxonomic classification |
| [IQ-TREE](http://www.iqtree.org/) | Phylogenetic inference |
| [CoverM](https://github.com/wwood/CoverM) | Genome coverage and abundance |
| [antiSMASH](https://antismash.secondarymetabolites.org/) | BGC annotation |
| [Prodigal](https://github.com/hyattpd/Prodigal) | Gene prediction |
| [DRAMMA](https://github.com/BorisKovalev/DRAMMA) | ARG annotation |
| [HMMER](http://hmmer.org/) | Sequence homology search |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Sequence similarity search |
| [TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/) | Transmembrane topology prediction |
| [DRAM](https://github.com/WrightonLabCSU/DRAM) | Metabolic/functional annotation |

> **Note:** The pipeline was developed for a SLURM HPC environment. Thread counts reference `${SLURM_CPUS_PER_TASK}` throughout `Script.sh`. Adjust these values to match your system.

---

## Usage

1. **Clone this repository**
   ```bash
   git clone https://github.com/HuihongLi/CaveMAGs.git
   cd CaveMAGs
   ```

2. **Install dependencies** — install each tool listed above into your environment (conda, modules, or manual installation).

3. **Download raw data**
   ```bash
   prefetch --option-file Output/1.sra_list.txt --max-size 100G
   ```

4. **Run the pipeline** — follow the numbered steps in [`Script.sh`](Script.sh). Each step contains self-contained shell commands. Substitute shell variables (e.g., `${SAMPLE}`, `${READ1}`, `${READ2}`, `${SCRATCH}`) with paths appropriate for your system. The `F` and `T` values for fastp hard-trimming (5′ and 3′ end bases to remove per sample) are provided in `Output/2.samples.tsv`.
