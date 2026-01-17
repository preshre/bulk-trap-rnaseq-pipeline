# RNA-seq Analysis Pipeline

A comprehensive R-based pipeline for differential expression analysis of RNA-seq data from Salmon quantification files.

---

## 0. Upstream Processing: FASTQ to Quantification

Before running this R pipeline, raw FASTQ files must be processed to generate Salmon quantification files. We used the **nf-core/rnaseq** pipeline for this step.

### Prerequisites for Upstream Processing

| Software | Purpose |
|----------|---------|
| Nextflow | Workflow management |
| Singularity/Docker | Container runtime |
| nf-core/rnaseq | RNA-seq processing pipeline |

### Reference Genome Files

Download mouse reference genome (mm39) from UCSC:
```bash
# Download reference FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz

# Download gene annotation GTF
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
```

### Sample Sheet Format

Create a CSV file (`sample.csv`) with the following format:
```csv
sample,fastq_1,fastq_2,strandedness
Sample1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz,auto
Sample2,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz,auto
```

### Running nf-core/rnaseq Pipeline

```bash
nextflow run nf-core/rnaseq \
    --input /path/to/sample.csv \
    --outdir /path/to/results/ \
    --fasta /path/to/mm39.fa.gz \
    --gtf /path/to/mm39.ncbiRefSeq.gtf.gz \
    --remove_ribo_rna \
    --save_merged_fastq \
    --trimmer fastp \
    --aligner star_salmon \
    -profile singularity \
    -bg
```

### Pipeline Parameters Explained

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--input` | sample.csv | Sample sheet with FASTQ paths |
| `--outdir` | results/ | Output directory |
| `--fasta` | mm39.fa.gz | Reference genome FASTA |
| `--gtf` | mm39.ncbiRefSeq.gtf.gz | Gene annotation GTF |
| `--remove_ribo_rna` | - | Remove ribosomal RNA reads |
| `--save_merged_fastq` | - | Save merged FASTQ files |
| `--trimmer` | fastp | Use fastp for adapter trimming |
| `--aligner` | star_salmon | Use STAR alignment + Salmon quantification |
| `-profile` | singularity | Use Singularity containers |
| `-bg` | - | Run in background |

### Optional: Contaminant Screening

For contaminant detection with Kraken2:
```bash
nextflow run nf-core/rnaseq \
    --input /path/to/sample.csv \
    --outdir /path/to/results/ \
    --fasta /path/to/mm39.fa.gz \
    --gtf /path/to/mm39.ncbiRefSeq.gtf.gz \
    --remove_ribo_rna \
    --save_merged_fastq \
    --trimmer fastp \
    --aligner star_salmon \
    --contaminant_screening kraken2_bracken \
    --kraken_db /path/to/kraken2/database \
    -profile singularity \
    -bg
```

### Output Structure from nf-core/rnaseq

After the pipeline completes, the relevant output for downstream analysis:
```
results/
└── star_salmon/
    ├── Sample1/
    │   └── quant.sf          # Salmon quantification (input for R pipeline)
    ├── Sample2/
    │   └── quant.sf
    ├── ...
    └── tx2gene.tsv           # Transcript-to-gene mapping
```

### Expected Runtime
- **6 samples:** 4-8 hours (depending on sequencing depth and cluster resources)
- **Storage:** ~50-100 GB per run

---

## 1. System Requirements

### Operating System
- **Tested on:** Rocky Linux 9.6 (Blue Onyx) - x86_64 architecture
- **Note:** This pipeline has been tested **only on Linux**. Compatibility with macOS or Windows is not guaranteed.

### Software Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| Conda/Miniconda | 25.3.1+ | Environment management |
| R | 4.2.3 | Statistical computing |
| Bioconductor | 3.16 | Bioinformatics packages |

### R Package Dependencies

#### Bioconductor Packages
| Package | Version | Purpose |
|---------|---------|---------|
| DESeq2 | 1.38.3 | Differential expression analysis |
| tximport | 1.26.1 | Import transcript-level quantifications |
| clusterProfiler | 4.6.2 | GO enrichment analysis |
| org.Mm.eg.db | 3.16.0 | Mouse annotation database |

#### CRAN Packages
| Package | Version | Purpose |
|---------|---------|---------|
| ggplot2 | 3.5.1 | Data visualization |
| dplyr | 1.1.4 | Data manipulation |
| readr | 2.1.5 | File reading |
| pheatmap | 1.0.12 | Heatmap visualization |
| ggrepel | 0.9.5 | Text label repulsion |
| scales | 1.3.0 | Scale functions for visualization |
| readxl | 1.4.5 | Excel file reading |
| svglite | 2.1.3 | SVG graphics device |
| ggbreak | 0.1.5 | Axis breaks in ggplot2 |
| matrixStats | 1.3.0 | Matrix statistics |
| BiocManager | 1.30.25 | Bioconductor package management |

#### Optional Packages
| Package | Version | Purpose |
|---------|---------|---------|
| RRHO2 | 1.0 | Rank-Rank Hypergeometric Overlap (optional) |

### Hardware Requirements
- **Minimum RAM:** 8 GB (16 GB recommended for large datasets)
- **Storage:** ~5 GB for environment + space for output files
- No specialized hardware required

---

## 2. Installation Guide

### Step 1: Install Miniconda (if not already installed)
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

### Step 2: Create the Conda Environment
```bash
# Create environment with R 4.2.3
conda create -n rnaseq-r-env -c conda-forge -c bioconda r-base=4.2.3

# Activate the environment
conda activate rnaseq-r-env
```

### Step 3: Install R Packages

#### Option A: Install via Conda (Recommended for core packages)
```bash
conda install -c conda-forge -c bioconda \
    bioconductor-deseq2=1.38.0 \
    bioconductor-annotationdbi \
    r-ggplot2=3.5.1 \
    r-pheatmap \
    r-ggrepel \
    r-scales \
    r-readxl \
    r-svglite \
    r-matrixstats
```

#### Option B: Install remaining packages from within R
```bash
# Start R
R
```

Then run in R:
```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Set Bioconductor version
BiocManager::install(version = "3.16")

# Install Bioconductor packages
BiocManager::install(c(
    "DESeq2",
    "tximport",
    "clusterProfiler",
    "org.Mm.eg.db"
))

# Install CRAN packages
install.packages(c(
    "ggplot2",
    "dplyr",
    "readr",
    "pheatmap",
    "ggrepel",
    "scales",
    "readxl",
    "svglite",
    "ggbreak",
    "matrixStats"
))

# Optional: Install RRHO2 for RRHO analysis
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("RRHO2/RRHO2")
```

### Step 4: Verify Installation
```bash
conda activate rnaseq-r-env
R --slave -e 'library(DESeq2); library(tximport); library(clusterProfiler); library(ggplot2); cat("All packages loaded successfully!\n")'
```

### Typical Installation Time
- **Conda environment creation:** 5-10 minutes
- **R package installation:** 15-30 minutes (depending on network speed)
- **Total:** ~30-45 minutes on a standard desktop computer

---

## 3. Demo

### Demo Data
The pipeline requires Salmon quantification files (`quant.sf`) organized as:
```
base_dir/
├── sample1/
│   └── quant.sf
├── sample2/
│   └── quant.sf
└── tx2gene.tsv
```

### Running the Demo
```bash
# Navigate to the scripts directory
cd /path/to/final_scripts

# Activate the environment
conda activate rnaseq-r-env

# Run with adjusted p-value (default)
Rscript run_pipeline.R

# Or run with both padj and raw p-value
Rscript run_pipeline.R --both
```

### Expected Output
```
results_20260116_131500/
├── vTRAP_BO_vTRAP_Trained/
│   ├── PCA_plot.pdf
│   ├── volcano_plot.pdf
│   ├── MA/
│   │   ├── MA_plot_*.pdf
│   │   └── DE_results_*.csv
│   ├── GO/
│   │   ├── GO_bubbles_BP_*.pdf
│   │   ├── GO_bubbles_MF_*.pdf
│   │   ├── GO_bubbles_CC_*.pdf
│   │   └── GO_*_results.csv
│   ├── DEG/
│   │   ├── DESeq2_full_results.csv
│   │   ├── Top2000_genes_by_padj.csv
│   │   ├── DEGs_all_padj.csv
│   │   ├── DEGs_up_in_*.csv
│   │   ├── Heatmap_TopDEGs.pdf
│   │   └── Heatmap_TopDEGs.svg
│   ├── RRHO/
│   │   └── *_RRHO_vector.rds
│   └── ReplicateCorrelations/
│       ├── correlation_heatmap_*.pdf
│       └── correlation_*.csv
├── vTRAP_BO_vTRAP_Trained_pvalue/
│   └── ... (same structure)
└── RRHO/
    └── vTRAP_vs_iPKR/
        ├── plots/
        │   ├── RRHO2_heatmap.pdf
        │   └── RRHO2_vennDiagram_*.pdf
        └── csv/
            └── genes_overlap_*.csv
```

### Expected Run Time
- **Single comparison:** 5-10 minutes
- **All 5 comparisons with --both:** ~29 minutes (actual: 28.81 minutes)
- Tested on: Rocky Linux 9.6 with 64 GB RAM on SLURM cluster

---

## 4. Instructions for Use

### Configuration
Edit `config.R` to customize your analysis:

#### Set Data Paths
```r
BASE_DIRS <- c(
    "/path/to/salmon/output1/star_salmon",
    "/path/to/salmon/output2/star_salmon"
)
OUTPUT_DIR <- "/path/to/results"
```

#### Adjust DEG Thresholds
```r
LOG2FC_THRESHOLD <- 0.58   # log2(1.5) fold change
PADJ_CUTOFF      <- 0.05   # Adjusted p-value cutoff
PVALUE_CUTOFF    <- 0.05   # Raw p-value cutoff
```

#### Define Sample Comparisons
```r
SAMPLES_my_comparison <- list(
    samples = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6"),
    conditions = c(rep("GroupA", 3), rep("GroupB", 3)),
    named_vector = c(
        "sample1" = "GroupA", "sample2" = "GroupA", "sample3" = "GroupA",
        "sample4" = "GroupB", "sample5" = "GroupB", "sample6" = "GroupB"
    ),
    label = "GroupA_vs_GroupB",
    ref_level = "GroupA"  # Reference/baseline group
)

# Add to comparisons list
ALL_COMPARISONS <- list(
    SAMPLES_my_comparison
)
```

### Running the Pipeline
```bash
# Activate environment
conda activate rnaseq-r-env

# Navigate to scripts
cd /path/to/final_scripts

# Option 1: Run with adjusted p-value (default)
Rscript run_pipeline.R

# Option 2: Run with raw p-value
Rscript run_pipeline.R --use-pvalue

# Option 3: Run both analyses
Rscript run_pipeline.R --both
```

### Pipeline Components

| Script | Function |
|--------|----------|
| `config.R` | Central configuration |
| `run_pipeline.R` | Master pipeline script |
| `correlation_analysis.R` | Replicate correlation heatmaps |
| `pca_volcano_plot.R` | PCA and volcano plots |
| `generate_DE_plots.R` | MA plots |
| `go_analysis.R` | GO enrichment analysis |
| `deg_analysis.R` | DEG lists and heatmaps |
| `rrho_analysis.R` | RRHO vector generation |
| `rrho_compute.R` | RRHO2 computation |
| `run_pipeline_job.sh` | SLURM job submission script |

### Analysis Parameters

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| log2FC threshold | 0.58 | Equivalent to 1.5-fold change |
| padj cutoff | 0.05 | Adjusted p-value significance |
| pvalue cutoff | 0.05 | Raw p-value significance |
| Top GO terms | 3 | Per direction (up/down) per ontology |
| baseMean minimum | 10 | Minimum expression for DEG |
| Top DEGs heatmap | 30 | Number of genes in heatmap |
| VST transform | blind=TRUE | For correlation analysis |

### DEG Analysis Output Files

| File | Description |
|------|-------------|
| `DESeq2_full_results.csv` | Complete results for all genes |
| `Top2000_genes_by_padj.csv` | Top 2000 genes ranked by significance |
| `DEGs_all_padj.csv` | All DEGs passing thresholds |
| `DEGs_up_in_{condition}.csv` | Upregulated genes per condition |
| `DEGs_top30.csv` | Top 30 DEGs used in heatmap |
| `Heatmap_TopDEGs.pdf` | Clustered heatmap of top DEGs |
| `Heatmap_TopDEGs.svg` | SVG version for publication |

---

## 5. Reproduction Instructions

To reproduce the results in the manuscript:

1. **Clone/copy the scripts** to your working directory
2. **Set up the conda environment** as described in Section 2
3. **Configure paths** in `config.R` to point to your Salmon output files
4. **Run the pipeline:**
   ```bash
   conda activate rnaseq-r-env
   cd final_scripts
   Rscript run_pipeline.R --both
   ```
5. **Outputs** will be generated in the specified `OUTPUT_DIR`

### Input Data Requirements
- Salmon quantification files (`quant.sf`) for each sample
- `tx2gene.tsv` file mapping transcript IDs to gene symbols
- Gene marker files in `../gene_data/` (for volcano plots):
  - `Glial marker genes.xlsx`
  - `Excitatory neuron marker genes.xlsx`
  - `Neuron marker genes.xlsx`
  - `PFC_Plasticity_Genes.xlsx`

---

## Troubleshooting

### Common Issues

1. **Package not found errors:**
   ```bash
   conda activate rnaseq-r-env
   R -e 'BiocManager::install("missing_package")'
   ```

2. **Cairo/PDF device errors:**
   ```bash
   conda install -c conda-forge r-cairo
   ```

3. **Memory errors:**
   - Reduce the number of comparisons in `ALL_COMPARISONS`
   - Run comparisons sequentially

### Contact
For issues related to this pipeline, please check the log files in the output directory.

---

## 6. Running on HPC Cluster (SLURM)

For large datasets or when running on a High-Performance Computing (HPC) cluster, use the provided SLURM job script.

### SLURM Job Script

The `run_pipeline_job.sh` script is configured for SLURM submission:

```bash
#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --partition=extended-40core
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
```

### Submitting a Job

```bash
# Navigate to scripts directory
cd /path/to/final_scripts

# Create logs directory (required)
mkdir -p logs

# Submit the job
sbatch run_pipeline_job.sh
```

### Monitoring Job Status

```bash
# Check if job is running
squeue -u $USER

# Check specific job
squeue -j JOBID

# Watch output log in real-time
tail -f logs/pipeline_JOBID.out

# Check detailed pipeline log
tail -f logs/pipeline_run_JOBID.log
```

### After Job Completes

```bash
# Check job exit status
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed,MaxRSS

# View full output
cat logs/pipeline_JOBID.out

# Check for errors
cat logs/pipeline_JOBID.err
```

### Canceling a Job

```bash
scancel JOBID
```

### SLURM Resource Recommendations

| Dataset Size | Memory | Time | CPUs |
|--------------|--------|------|------|
| Small (≤6 samples, 1-2 comparisons) | 32 GB | 4 hours | 2 |
| Medium (6-12 samples, 3-5 comparisons) | 64 GB | 8 hours | 4 |
| Large (>12 samples, >5 comparisons) | 128 GB | 12+ hours | 8 |

### Log Files

After job completion, logs are stored in `logs/`:

| File | Description |
|------|-------------|
| `pipeline_JOBID.out` | SLURM standard output |
| `pipeline_JOBID.err` | SLURM error output |
| `pipeline_run_JOBID.log` | Detailed R pipeline log |

### Customizing the Job Script

Edit `run_pipeline_job.sh` to:
- Change partition: `#SBATCH --partition=your_partition`
- Adjust memory: `#SBATCH --mem=128G`
- Modify time limit: `#SBATCH --time=24:00:00`
- Select specific nodes: `#SBATCH --nodelist=node01`

---

