# 1. Download SRA files and convert them to fastq format
prefetch --option-file 1.sra_list.txt --output-directory /gpfs/radev/scratch/mcdougal/hl993/Cave --max-size 100G 
mkdir -p raw
find . -type f -name "*.sra" | while read -r sra_file; do
    base=$(basename "$sra_file" .sra)
    fasterq-dump --split-files --threads 64 --outdir raw "$sra_file"
    if command -v pigz &> /dev/null; then
        pigz -p 64 raw/"$base"_1.fastq
        pigz -p 64 raw/"$base"_2.fastq
    else
        gzip raw/"$base"_1.fastq
        gzip raw/"$base"_2.fastq
    fi
done

# 2. Quality control using FastQC and Fastp
fastqc -t 32 -o fastqc_results *.fastq.gz
multiqc fastqc_results -o fastqc_results

fastp -i raw/${SAMPLE}_1.fastq.gz -I raw/${SAMPLE}_2.fastq.gz \
      -o clean/${SAMPLE}_1.clean.fq.gz -O clean/${SAMPLE}_2.clean.fq.gz \
      --detect_adapter_for_pe \
      -q 20 -u 25 -n 5 \
      --length_required 50 --length_limit 300 \
      -f ${F} -t ${T} -F ${F} -T ${T} \
      --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
      --trim_poly_g --trim_poly_x \
      --dedup --dup_calc_accuracy 3 \
      -w 12 \                    
      -h clean/${SAMPLE}_fastp.html -j clean/${SAMPLE}_fastp.json

# 3. Metagenomic assembly using MEGAHIT
megahit \
  -1 "${READ1}" -2 "${READ2}" \
  --presets meta-sensitive \
  --min-contig-len 1500 \
  --tmp-dir "${SCRATCH}" \
  -t ${SLURM_CPUS_PER_TASK} \
  -o "${OUTDIR}"

# 4. Alignment using Bowtie2
bowtie2 -x ${FA} -1 ${R1} -2 ${R2} -p 16 \
 | samtools view -bS -@ 4 - \
 | samtools sort -@ 8 -o ${OUT}

samtools index ${OUT}

# 5. Binning using MetaBAT2, SemiBin2, and MetaDecoder

metabat2 \
    --inFile "${FA}" \
    --outFile "${TMP_RAW}/${SAMPLE}.bin" \
    --abdFile "${ABD_FILE}" \
    --numThreads ${SLURM_CPUS_PER_TASK} \
    --minContig 2500 \
    --minClsSize 200000 \
    --unbinned

SemiBin2 single_easy_bin \
    --input-fasta "${FA}" \
    --input-bam   "${BAM}" \
    --output      "${TMP_RAW}" \
    --threads     ${SLURM_CPUS_PER_TASK} \
    --environment global \
    --compression none

metadecoder coverage -b "${BAM}" -o "${COV}"
metadecoder seed -f "${FA}" -o "${SEED}" --threads ${SLURM_CPUS_PER_TASK}
metadecoder cluster -f "${FA}" -c "${COV}" -s "${SEED}" -o "${TMP_RAW}/md"

# 6. Refinement using MAGScoT
Rscript "${MAGS_DIR}/MAGScoT.R" \
        -i "${WORK_DIR}/${SAMPLE}.contig2bin.tsv" \
        --hmm "${WORK_DIR}/${SAMPLE}.hmm" \
        --out "${WORK_DIR}/${SAMPLE}" \
        --max_cont 1 \
        --threshold 0.5

# 7. Quality control using CheckM
checkm2 predict \
    -i "${IN_DIR}" \
    -o "${OUT_DIR}" \
    -x fa \
    --threads ${SLURM_CPUS_PER_TASK}

# 8. Dereplication using dRep
dRep dereplicate "${DREP_OUT}" \
    -g "${TMP_LINKS}"/*.fa "${TMP_LINKS}"/*.fasta \
    -p "${SLURM_CPUS_PER_TASK}" \
    -pa 0.90 -sa 0.95 \
    --S_algorithm fastANI

# 9. Taxonomic classification using GTDB-Tk
gtdbtk classify_wf \
    --genome_dir "$MAG_DIR" \
    --out_dir "$OUT_DIR" \
    --extension fa \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --full_tree \
    --skip_ani_screen 

# 10. Phylogenetic analysis using IQ-TREE
iqtree -s "$MSA" \
       -m MFP \
       -B 1000             \
       -alrt 1000          \
       -T ${SLURM_CPUS_PER_TASK} \
       --mem 0 \
       --prefix "${OUT_DIR}/${PREFIX}"

# 11. Genome coverage analysis using CoverM
coverm genome \
  --coupled "${READ1}" "${READ2}" \
  --genome-fasta-files "${GENOME_FILES[@]}" \
  -t ${SLURM_CPUS_PER_TASK} \
  --min-covered-fraction 0 \
  -m mean relative_abundance covered_fraction variance count reads_per_base length \
  -o "${OUTPUT_FILE}" 

# 12. AntiSMASH analysis
antismash \
    --taxon bacteria \
    --cpus ${SLURM_CPUS_PER_TASK} \
    --output-dir "$OUTDIR" \
    --output-basename "$SAMPLE" \
    --html-title "antiSMASH Analysis: ${SAMPLE}" \
    --genefinding-tool prodigal \
    --cc-mibig \
    --cb-general \
    --cb-subclusters \
    --cb-knownclusters \
    --asf \
    --clusterhmmer \
    --tigrfam \
    --pfam2go \
    --rre \
    --tfbs \
    "${WORKDIR}/${BASENAME}" 

# 13. DRAMMA pipeline for ARGs annotation
prodigal -i "$FASTA" \
         -o "${OUTDIR}/${SAMPLE}.gff" \
         -a "${OUTDIR}/${SAMPLE}.faa" \
         -d "${OUTDIR}/${SAMPLE}.ffn" \
         -p "$MODE"

python run_DRAMMA_pipeline.py \
    --dif_format_paths \
        "${SAMPLE_DIR}/${SAMPLE}.faa" \
        "${SAMPLE_DIR}/${SAMPLE}.gff" \
        "${SAMPLE_DIR}/${SAMPLE}.ffn" \
        "${SAMPLE_DIR}/${SAMPLE}.fa" \
    --hmmer_path  "$HMMSEARCH" \
    --mmseqs_path "$MMSEQS" \
    --tmhmm_path  "$TMHMM_BIN/tmhmm" \
    --feature_dir "${OUTDIR}/features" \
    --threshold_list 1e-5 \
    --gene_window 20 \
    --nucleotide_window 20000 \
    --ncpus "$SLURM_CPUS_PER_TASK" \
    --output_file "${OUTDIR}/${SAMPLE}.dramma.pkl"

# 14. Annotation using DRAM
DRAM-setup.py prepare_databases \
    --output_dir "${OUT_DIR}" \
    --select_db kofam_hmm \
    --select_db kofam_ko_list \
    --skip_uniref \
    --select_db pfam \
    --select_db pfam_hmm \
    --select_db dbcan \
    --select_db dbcan_fam_activities \
    --select_db dbcan_subfam_ec \
    --select_db peptidase \
    --select_db genome_summary_form \
    --select_db module_step_form \
    --select_db function_heatmap_form \
    --select_db amg_database \
    --select_db etc_module_database \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose 

DRAM.py annotate \
    -i "${WORKDIR}/${BASENAME}" \
    -o "$OUTDIR" \
    --threads ${SLURM_CPUS_PER_TASK}
