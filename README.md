# BioVarFlow

Germline variant calling and reporting with standardized QC and functional annotation — designed for research, not diagnostic use.

---
## Overview

**BioVarFlow** is a Nextflow-based pipeline designed to produce high-quality germline variant calls and structured reports for research and data integration workflows.

Main Features:

- Alignment using **BWA-MEM2**
- Variant calling with:
  - DeepVariant
  - GATK HaplotypeCaller
- Consensus variant integration
- ACMG region filtering
- Coverage and QC analysis
- Functional annotation (VEP)
- Excel and HTML reporting outputs

---
## Inputs

The pipeline runs from FASTQ files using a standardized samplesheet.

### Required Parameters

- `--input` — Samplesheet CSV
- `--outdir` — Output directory

### Example samplesheet

```csv
patient,assay,sex,status,sample,fastq_1,fastq_2,lane
Patient001,WES,XX,0,Sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz,1
Patient002,WES,XY,0,Sample_B,data/sample_B_R1.fastq.gz,data/sample_B_R2.fastq.gz,1
Patient003,WES,XX,0,Sample_C,data/sample_C_R1.fastq.gz,data/sample_C_R2.fastq.gz,1
```

The CSV should include:
- **patient**: Patient identifier
- **assay**: Assay type (`WES` or `WGS`)
- **sex**: Sex (`XX` for female, `XY` for male, or `0` if unknown)
- **status**: Always set to `0`
- **sample**: Unique sample identifier
- **fastq_1**: Path to forward reads file
- **fastq_2**: Path to reverse reads file
- **lane**: Sequencing lane number

---

## Outputs

Results are organized by sample in the output directory:

```
outdir/
└── {SAMPLE}/
    ├── vcf/          # Variant files (filtered and annotated)
    ├── qc/           # Quality metrics and coverage statistics
    └── reports/      # Excel and HTML reports
```

**What you'll find:**

- **VCF files**:
  - Filtered and annotated VCF files
  - Consensus integration of DeepVariant and HaplotypeCaller calls
- **QC metrics**:
  - Coverage analysis (20x, 30x, 50x, 100x thresholds)
  - Coverage gaps identification
  - Alignment statistics (mapped reads, duplicates)
  - Read balance metrics (R1/R2 ratios, strand balance)
  - Sex determination from sequencing data
  - Variant statistics (transition/transversion ratios)
- **Reports**: Easy-to-read Excel summaries with variant classifications and clinical annotations

---

## How to Run

### Using the Platform UI

1. Navigate to the pipeline in your platform interface
2. Specify your **input samplesheet** path
3. Specify your **output directory** path
4. Click **Run Workflow**

### Command Line

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results
```

---
# Benchmarking

This workflow was benchmarked using Genome in a Bottle (GIAB) reference samples
to validate performance within clinically relevant ACMG Secondary Findings regions.

## Benchmarking Overview

Datasets:
- HG003 (GIAB)
- HG002 (GIAB)

Reference:
- GRCh38

Target Regions:
- ACMG SF v3.3 BED

Pipeline Configuration:
- Aligner: BWA-MEM2
- Variant Callers: DeepVariant + GATK HaplotypeCaller
- Strategy: Consensus variant integration (DV+HC)
- Filtering:
  - ACMG region restriction
  - Depth and genotype quality thresholds
  - PASS variants only

## Performance Summary (ACMG Regions)
### SNPs

| Metric | Result |
|---|---|
| Precision | 94.8% |
| Recall | 99.5% |
| F1 Score | 97.1% |

### INDELs

| Metric | Result |
|---|---|
| Precision | 71% |
| Recall | 100% |
| F1 Score | 83% |

Benchmarking was performed against GIAB high-confidence truth sets using hap.py.

## Reproducibility

To reproduce validated results, use the configuration shown in `docs/usage.md`.

## Support

For questions or assistance:
- Contact your bioinformatics team
- Email: khatib@bionl.ai

**Pipeline version**: 1.1.0
