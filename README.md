# aDNA-snakemake

A Snakemake pipeline for **ancient DNA (aDNA)** shotgun / capture data. It takes raw FASTQ files and produces damage-masked, filtered BAMs, QC reports, mitochondrial haplogroups, genetic sex estimates and (optionally) pseudo-haploid genotypes at a SNP panel such as 1240k.

Every step runs in its own small conda environment (`envs/*.yaml`), which Snakemake builds automatically, so the only thing you install by hand is Snakemake itself.

---

## Table of contents

1. [What the pipeline does](#1-what-the-pipeline-does)
2. [Project layout](#2-project-layout)
3. [Requirements](#3-requirements)
4. [Installation](#4-installation)
5. [Prepare your inputs](#5-prepare-your-inputs)
6. [Configure the pipeline](#6-configure-the-pipeline)
7. [Run the pipeline](#7-run-the-pipeline)
8. [Outputs](#8-outputs)
9. [Interpreting the summary table](#9-interpreting-the-summary-table)
10. [Running on a cluster (SLURM)](#10-running-on-a-cluster-slurm)
11. [Common recipes](#11-common-recipes)
12. [Troubleshooting](#12-troubleshooting)
13. [Method notes and caveats](#13-method-notes-and-caveats)

---

## 1. What the pipeline does

```
FASTQ ──► FastQC
   │
   └──► AdapterRemoval (trim + merge pairs)
            │
            └──► bwa aln + samse (aDNA settings: -l 1024 -n 0.01 -o 2)
                     │
                     └──► filter (mapped, MAPQ, min length)
                              │
                              └──► dedup (samtools fixmate + markdup -r)
                                       │
                                       ├──► mapDamage (damage profiles)
                                       │
                                       └──► trimBam (soft-clip terminal bases)  ──► FINAL BAM
                                                    │
                                                    ├──► flagstat / idxstats / summary.tsv / MultiQC
                                                    ├──► bcftools -> HaploGrep3   (mt haplogroup)
                                                    └──► samtools mpileup -> pileupCaller -> (merge ds + ss) -> EIGENSTRAT
```

| Step | Rule(s) | Tool |
|---|---|---|
| Reference indexing | `bwa_index`, `faidx` | bwa, samtools |
| Raw read QC | `fastqc` | FastQC |
| Adapter trimming / pair merging | `trim_merge` | AdapterRemoval |
| Mapping | `bwa_aln`, `bwa_samse` | bwa |
| Filtering | `filter_bam` | samtools |
| Duplicate removal | `dedup` | samtools fixmate/markdup |
| Damage profile | `mapdamage` | mapDamage2 |
| Damage masking | `trim_ends` | bamUtil trimBam |
| Statistics | `bam_stats`, `sample_stats`, `summary`, `multiqc` | samtools, MultiQC |
| mtDNA haplogroup | `mt_vcf`, `haplogrep` | bcftools, HaploGrep3 |
| Genotyping | `pileupCaller`, `merge_geno` | sequenceTools, EIGENSOFT |

---

## 2. Project layout

```
adna-snakemake/
├── Snakefile                 # the workflow
├── environment.yml           # driver env (Snakemake + pandas)
├── config/
│   ├── config.yaml           # all parameters  <-- edit this
│   └── samples.tsv           # sample sheet    <-- edit this
├── envs/                     # one conda env per tool group (auto-built)
│   ├── map.yaml  trim.yaml  qc.yaml  damage.yaml
│   ├── bamutil.yaml  mt.yaml  geno.yaml  eigensoft.yaml
└── scripts/
    └── sample_stats.sh       # per-sample summary line
```

> **Important:** run Snakemake from the **project root** (the folder containing `Snakefile`). The Snakefile loads `config/config.yaml`, and relative paths in `samples.tsv` and `config.yaml` are resolved from the directory you launch from.

---

## 3. Requirements

- Linux or macOS (Apple Silicon: some bioconda packages need Rosetta / `CONDA_SUBDIR=osx-64`; Linux is strongly recommended)
- [Conda / Miniforge / Mamba](https://github.com/conda-forge/miniforge) installed
- Internet access on the first run (conda packages, plus the HaploGrep3 tree download)
- A reference genome FASTA (e.g. **hs37d5**), with write access to its folder (indexes are created next to it)
- Disk space: roughly 3–5x the size of your FASTQ files for intermediate BAMs, plus ~8 GB for the reference and its indexes
- Memory: 16 GB is typically enough for human aDNA samples; bwa indexing of a human reference needs ~5–6 GB

---

## 4. Installation

```bash
# 1. Get the project and enter it
cd adna-snakemake

# 2. Create the driver environment (Snakemake only)
mamba env create -f environment.yml        # or: conda env create -f environment.yml

# 3. Activate it
conda activate adna-snakemake

# 4. Check
snakemake --version                        # should be >= 8.4
```

The tool environments in `envs/` are **not** created now. Snakemake builds them the first time you run with conda enabled (this can take 10–20 minutes once; they are cached afterwards in `.snakemake/conda/`).

---

## 5. Prepare your inputs

### 5.1 Sample sheet: `config/samples.tsv`

Tab-separated with a header and exactly these four columns:

| Column | Meaning |
|---|---|
| `sample` | Unique sample ID (no `/`, no whitespace). Used in all output filenames. |
| `fq1` | Path to FASTQ R1 (`.fastq.gz`) |
| `fq2` | Path to FASTQ R2, or `NA` (or empty) for single-end data |
| `library` | `ds` = double-stranded library prep, `ss` = single-stranded library prep |

Example:

```
sample	fq1	fq2	library
IND001	data/IND001_R1.fastq.gz	data/IND001_R2.fastq.gz	ds
IND002	data/IND002_R1.fastq.gz	NA	ss
```

Rules enforced by the Snakefile (it will stop with a clear error otherwise):

- All four columns must be present.
- `library` must be `ds` or `ss`.
- Sample IDs must be unique.
- Lines starting with `#` are ignored (handy to comment out samples).

Why `library` matters: single-stranded libraries have a different damage pattern, so they are genotyped with `pileupCaller --singleStrandMode` in a separate run, then merged with the double-stranded set.

**One sample, multiple lanes/runs?** Concatenate the FASTQs per sample first (`cat lane1_R1.fastq.gz lane2_R1.fastq.gz > IND001_R1.fastq.gz`, and the same for R2), since the pipeline expects one FASTQ (pair) per sample.

### 5.2 Reference genome

Provide a FASTA, e.g. `hs37d5.fa`. The bwa index (`.amb .ann .bwt .pac .sa`) and `.fai` are generated automatically next to the FASTA if missing. To save time, place pre-built indexes beside it. The reference folder must be writable otherwise.

Contig names must match the settings in `config.yaml` (`mt_name`, `x_name`, `y_name`):

| Reference | MT | X | Y |
|---|---|---|---|
| hs37d5 / GRCh37 (default) | `MT` | `X` | `Y` |
| hg19 (UCSC) / GRCh38 with `chr` prefix | `chrM` | `chrX` | `chrY` |

Note: the haplogroup step uses the **rCRS** mitochondrial sequence (phylotree-rcrs tree). hs37d5 uses rCRS; the UCSC hg19 `chrM` is NOT rCRS, so use hs37d5 or GRCh38 for haplogrouping.

### 5.3 SNP panel (optional, for genotyping)

To get pseudo-haploid genotypes you need two files describing the same panel positions:

- `panel_bed`: BED file (chrom, 0-based start, end) of the SNP positions
- `panel_snp`: EIGENSTRAT `.snp` file of the panel

Contig names in both must match the reference (for 1240k on hs37d5: `1..22, X, Y`). If you leave either blank, genotyping is skipped.

---

## 6. Configure the pipeline

Edit `config/config.yaml`. The only setting you **must** change is `ref`.

| Key | Default | Description |
|---|---|---|
| `samples` | `config/samples.tsv` | Sample sheet path |
| `outdir` | `results` | Output directory |
| `ref` | `/path/to/hs37d5.fa` | **Reference FASTA (change this)** |
| `panel_bed` | `""` | Panel BED; blank = skip genotyping |
| `panel_snp` | `""` | Panel EIGENSTRAT `.snp`; blank = skip genotyping |
| `mt_name`, `x_name`, `y_name` | `MT`, `X`, `Y` | Contig names in the reference |
| `adapter1`, `adapter2` | Illumina TruSeq | Adapter sequences for AdapterRemoval |
| `min_len` | `30` | Minimum read length after trimming |
| `min_trim_quality` | `20` | Base quality for end trimming |
| `bwa_l`, `bwa_n`, `bwa_o` | `1024`, `0.01`, `2` | bwa aln settings for aDNA (seed disabled, relaxed edit distance) |
| `min_mapq` | `30` | Minimum mapping quality |
| `min_map_len` | `35` | Minimum alignment length on the reference |
| `min_bq` | `30` | Base quality threshold for pileupCaller |
| `trim_bp` | `2` | Terminal bases soft-clipped at both read ends (see below) |
| `run_haplogroup` | `true` | Set `false` to skip mtDNA haplogroup calling |
| `haplogrep_tree` | `phylotree-rcrs@17.1` | HaploGrep3 tree (downloaded on first use) |
| `mt_min_depth` | `3` | Minimum depth for mt variant calls |

### Choosing `trim_bp` (damage masking)

Set it according to how the libraries were treated:

| Library treatment | `trim_bp` |
|---|---|
| Non-UDG | 10 |
| UDG-half | 2 (default) |
| Full UDG | 0 (no trimming; BAM is just copied) |

`trim_bp` is global: if your samples have mixed treatments, run them in separate batches with separate sample sheets and `outdir`s.

You can also override any setting on the command line without editing the file:

```bash
snakemake --use-conda --cores 16 --config ref=/data/hs37d5.fa trim_bp=10 outdir=run2
```

---

## 7. Run the pipeline

Always do a **dry run** first to make sure the sample sheet and config are valid and to see what will be executed:

```bash
conda activate adna-snakemake
snakemake -n --use-conda            # Snakemake 8+: add --software-deployment-method conda if needed
```

Then run for real:

**Snakemake 8 or newer** (default from `environment.yml`):

```bash
snakemake --software-deployment-method conda --cores 16
```

**Snakemake 7:**

```bash
snakemake --use-conda --cores 16
```

(`--use-conda` also still works in Snakemake 8 as an alias.)

Useful options:

| Option | Purpose |
|---|---|
| `--cores N` | Total CPU cores to use (rules request 2–8 threads each; several samples run in parallel) |
| `--rerun-incomplete` | Redo jobs that were interrupted mid-way |
| `--keep-going` | Continue with other samples if one fails |
| `--conda-prefix ~/.conda-envs/adna` | Share built tool environments across projects |
| `--printshellcmds` | Print every shell command (good for debugging) |
| `--config key=value` | Override a config setting |
| `--forceall` | Redo everything (be careful) |
| `--dag \| dot -Tpdf > dag.pdf` | Draw the workflow graph (needs Graphviz) |

To build just the conda environments up front (useful before submitting to a cluster, or if compute nodes lack internet):

```bash
snakemake --software-deployment-method conda --conda-create-envs-only --cores 4
```

### Running a single sample or a single target

```bash
# Only one sample's final BAM
snakemake --use-conda --cores 8 results/final/IND001.final.bam

# Only the summary table
snakemake --use-conda --cores 8 results/stats/summary.tsv
```

### Test on a tiny dataset first

Before processing your full dataset, put 1–2 small samples (e.g. `zcat reads.fq.gz | head -n 400000 | gzip > small.fq.gz`) into `samples.tsv` and a small reference (e.g. just chr21 + MT with matching `x_name`/`y_name` adapted) to confirm everything works end to end.

---

## 8. Outputs

All results go under `outdir` (default `results/`):

```
results/
├── fastqc/<sample>/                    FastQC reports
├── trim/<sample>.merged.fq.gz          Trimmed (and merged) reads
├── trim/<sample>.settings              AdapterRemoval report
├── map/                                Temporary files (auto-deleted)
├── dedup/<sample>.dedup.bam            Deduplicated BAM (no damage masking)
├── dedup/<sample>.markdup.txt          Duplicate statistics
├── damage/<sample>/                    mapDamage output (misincorporation.txt, plots...)
├── final/<sample>.final.bam(.bai)      ★ Final damage-masked BAM
├── stats/<sample>.flagstat.txt         samtools flagstat
├── stats/<sample>.idxstats.txt         samtools idxstats
├── stats/summary.tsv                   ★ One-row-per-sample summary
├── stats/multiqc_report.html           ★ Aggregated QC report
├── mt/<sample>.vcf                     Mitochondrial variants
├── mt/<sample>.haplogrep.txt           ★ mtDNA haplogroup
├── geno/ds.* / ss.*                    Per-library EIGENSTRAT genotypes
├── geno/all_samples.{geno,snp,ind}     ★ Merged EIGENSTRAT dataset
└── logs/                               One log per rule and sample
```

(Genotype and haplogroup files only appear if those steps are enabled.)

---

## 9. Interpreting the summary table

`results/stats/summary.tsv` columns:

| Column | Meaning |
|---|---|
| `sample`, `lib` | Sample ID and library type |
| `raw_reads` | Reads in FASTQ R1 |
| `merged_reads` | Reads after trimming/merging |
| `mapped_q` | Reads passing mapping filters (mapped, MAPQ, length) |
| `endogenous_pct` | `mapped_q / merged_reads x 100` (proxy for endogenous DNA content) |
| `unique_reads` | Reads after duplicate removal |
| `dup_rate` | `1 - unique/mapped_q` |
| `mean_len` | Mean length of unique reads (aDNA is typically 40–70 bp) |
| `mean_depth` | Mean genome-wide depth in the final BAM |
| `cov_1x_pct` | Percentage of the genome covered at >= 1x |
| `Ry` | Fraction of X+Y reads that map to Y (Skoglund et al. 2013) |
| `sex` | `XY` if Ry > 0.075, `XX` if Ry < 0.016, `ambiguous` in between, `NA(low)` if fewer than 100 X+Y reads |
| `mt_depth` | Mean depth over the mitochondrial genome |

Check the damage plots in `damage/<sample>/` (and the 5' C→T frequency in `5pCtoT_freq.txt`): a clear C→T increase at 5' ends is evidence the DNA is authentically ancient. Values of roughly >= 10% at the first position are typical for non-UDG data.

---

## 10. Running on a cluster (SLURM)

Install the SLURM executor plugin once in the driver env:

```bash
pip install snakemake-executor-plugin-slurm
```

Then (Snakemake 8+):

```bash
snakemake --executor slurm --jobs 50 \
  --software-deployment-method conda \
  --default-resources slurm_account=<account> slurm_partition=<partition> mem_mb=8000 runtime=240
```

Notes:

- Build the conda environments on a login node first with `--conda-create-envs-only`.
- Run the `snakemake` command itself inside `screen`/`tmux` or as its own batch job.
- bwa and mapDamage can need more memory for deep samples: raise `mem_mb` with `--set-resources bwa_samse:mem_mb=16000`.

---

## 11. Common recipes

**Add samples later**: append rows to `samples.tsv` and rerun the same command. Only the new samples (and the summary / MultiQC / merged genotype files) are recomputed.

**Change a filter parameter** (e.g. `min_mapq`): rerun with `--rerun-triggers params` (Snakemake 8) or `--forcerun filter_bam` to redo from that step on.

**Skip haplogroup calling**: `--config run_haplogroup=false` (or set it in `config.yaml`).

**Skip genotyping**: leave `panel_bed` and `panel_snp` blank.

**Only ds or only ss libraries**: works as is; with a single library type the genotype set is simply copied to `all_samples.*` instead of being merged.

**Clean up large intermediates**: after you are happy with the results you can delete `results/dedup/*.bam` and `results/trim/*.fq.gz`. (Snakemake will want to regenerate them on the next run unless you use `snakemake --touch` or accept the rebuild.)

---

## 12. Troubleshooting

| Problem | Fix |
|---|---|
| `MissingInputException` for `config/config.yaml` | Run from the project root; the file must be at `config/config.yaml`. |
| `samples.tsv is missing columns` | Header must be tab-separated: `sample  fq1  fq2  library`. |
| `library column must be 'ds' or 'ss'` | Fix typos / trailing spaces in the `library` column. |
| `MissingInputException` for a FASTQ | Paths in `samples.tsv` are relative to where you launch Snakemake. Use absolute paths if unsure. |
| bwa index step fails with "Permission denied" | The reference directory is read-only. Copy the FASTA to a writable location and update `ref`. |
| Conda environment creation is very slow or fails | Use mamba (`snakemake --conda-frontend mamba`), check internet/proxy, and make sure the `bioconda` and `conda-forge` channels are reachable. |
| `haplogrep3` fails to download the tree | Compute node has no internet. Run once on a node with internet, or download the tree manually and point `haplogrep_tree` to it. |
| Empty `.vcf` / haplogroup fails | Too little mtDNA coverage. Check `mt_depth` in `summary.tsv`; lower `mt_min_depth` only if you accept noisier calls. |
| `sex` shows `NA(low)` | Fewer than 100 reads on X+Y; sequence deeper before relying on sex estimates. |
| `pileupCaller` / `mergeit` errors on contig names | Panel BED/SNP contig names do not match the reference (`chr1` vs `1`). Make them consistent. |
| Pipeline was killed mid-run | Rerun with `--rerun-incomplete`. If Snakemake reports a locked directory: `snakemake --unlock`. |
| Out of memory | Reduce `--cores`, or give heavy rules more memory (`--set-resources`). |

Each step writes its stdout/stderr to `results/logs/...`, which is the first place to look when a job fails.

---

## 13. Method notes and caveats

- **Merging and unmerged pairs.** For paired-end data, only collapsed reads (`.collapsed` + `.collapsed.truncated`) are kept; unmerged pairs are dropped, as is standard for short aDNA fragments. Single-end data are trimmed only.
- **Duplicate removal.** Because merged reads are treated as single-end, `samtools markdup` identifies duplicates by 5' position and strand. This is conservative and standard for aDNA, but it may also remove some independent molecules sharing a start site in very deep data.
- **Mapping.** `bwa aln` with `-l 1024 -n 0.01 -o 2` follows Schubert et al. (2012) for ancient DNA. Reads are mapped as single-end.
- **Damage masking.** Terminal bases are soft-clipped (not removed), so read lengths in `final.bam` are unchanged but clipped bases do not contribute to pileups.
- **Genotyping.** `pileupCaller --randomHaploid` samples one allele at random per site (pseudo-haploid). Single-stranded libraries use `--singleStrandMode`, so they are genotyped separately and merged with `mergeit`.
- **Sex determination** uses all reads on X and Y in the final BAM (Ry method) and is meant for human data aligned to a reference with standard X/Y contigs.
- **Reproducibility.** Tool versions are pinned loosely (`>=`) in `envs/*.yaml`. For fully reproducible runs, export the resolved environments (`conda env export`) or pin exact versions.

### References

- Schubert et al. 2012, *Improving ancient DNA read mapping against modern reference genomes*, BMC Genomics.
- Skoglund et al. 2013, *Accurate sex identification of ancient human remains using DNA shotgun sequencing*, J Archaeol Sci.
- Jónsson et al. 2013, *mapDamage2.0*, Bioinformatics.
- Mölder et al. 2021, *Sustainable data analysis with Snakemake*, F1000Research.

Gaurav Sablok \
gsablok@proton.me
