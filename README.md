# CaveMAGs

A comprehensive prokaryotic genomic catalog from 37 geographically diverse cave environments, built with an optimized metagenome-assembled genome (MAG) reconstruction pipeline.

---

## Overview

Cave microorganisms are unique extremophiles that have evolved in isolated, nutrient-limited environments and harbor exceptional metabolic capabilities. Yet knowledge of cave microbial diversity at the genomic level remains limited — previous studies have focused on individual caves and do not provide a global picture.

Here we present the first prokaryotic cave metagenomic catalog spanning 37 geographically diverse cave environments. Using an optimized genome reconstruction pipeline applied to 241 public metagenomic datasets, we recovered 3,837 medium-to-high quality metagenome-assembled genomes (MAGs). These were dereplicated into 1,979 species-level representative clusters covering 67 phyla of Bacteria (*n* = 1,858) and Archaea (*n* = 121).
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

The pipeline is documented as a series of shell commands in [`Script.sh`](Script.sh)
