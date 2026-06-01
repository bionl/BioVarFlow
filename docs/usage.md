# Usage Guide

---

## Pipeline Overview

**Bionl_Lean_Call** is a Nextflow workflow for germline variant calling,
quality control, consensus integration, and ACMG Secondary Findings reporting.

The pipeline runs from FASTQ files using a standardized samplesheet.

Main workflow steps:

1. Read preprocessing and alignment (BWA-MEM2)
2. Variant calling using DeepVariant and GATK HaplotypeCaller
3. Consensus variant integration
4. ACMG region filtering
5. Functional annotation (VEP)
6. QC and coverage analysis
7. Excel and HTML report generation

---

## Required Inputs

### Parameters

- `--input` — CSV samplesheet containing FASTQ paths
- `--outdir` — Output directory

---

## Samplesheet Format

Create a CSV file with the following columns:

```csv
patient,assay,sex,status,sample,fastq_1,fastq_2,lane
Patient001,WES,XX,0,Sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz,1
Patient002,WES,XY,0,Sample_B,data/sample_B_R1.fastq.gz,data/sample_B_R2.fastq.gz,1
```

**Column Descriptions:**
- `patient`: Patient identifier (can be the same for multiple samples)
- `assay`: Assay type - use `WES`(Whole Exome Sequencing) or `WGS` (Whole genome sequencing)
- `sex`: Biological sex - use `XX` (female), `XY` (male), or `0` (unknown)
- `status`: Disease status - `0` (unaffected) or `1` (affected) (Define always as 0)
- `sample`: Unique sample identifier
- `fastq_1`: Full path to R1 FASTQ file
- `fastq_2`: Full path to R2 FASTQ file
- `lane`: Sequencing lane number

#### Command Example

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results
```

---
## Output Structure

Results are organized by sample:

```
outdir/
├── SAMPLE1/
│   ├── vcf/
│   ├── qc/
│   └── reports/
└── pipeline_report.html
```

### Output Files Description

**VCF Directory:**
- Filtered consensus variant calls
- Annotated VCF files (VEP)
- ACMG Secondary Findings variants

**QC Directory:**
- Coverage statistics at 20x, 30x, 50x, 100x thresholds
- Per-exon coverage metrics
- Coverage gaps in target regions
- Sample quality summary (JSON)

**Reports Directory:**
- Excel report with variant classifications
- HTML report with interactive visualizations

---

## Complete Usage Examples

### Example: Run full pipeline with FASTQ files

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results \
```


---
# Validated Configuration (GIAB Benchmarking)

The workflow was benchmarked using Genome in a Bottle (GIAB) datasets
and achieved ≥95% precision within ACMG SF v3.3 regions using the
following validated configuration.

## Pipeline Components

- **Aligner:** BWA-MEM2
- **Variant Callers:** DeepVariant + GATK HaplotypeCaller
- **Variant Strategy:** Consensus integration (DV + HC)
- **Reference Genome:** GRCh38
- **Target Regions:** ACMG SF v3.3 BED

## Region & Variant Filtering

- Variants restricted to ACMG SF regions
- Only `PASS` variants retained
- GIAB high-confidence regions used during benchmarking

## Consensus Variant Filtering Strategy

Variants are filtered using caller-aware thresholds to balance sensitivity and precision during consensus integration.

### Shared Calls (DeepVariant + HaplotypeCaller)

Variants detected by both callers are retained using relaxed thresholds, as concordance between callers increases confidence.

Criteria:
- DP ≥ 10
- GQ ≥ 10

---

### HaplotypeCaller-Only Calls

Variants detected only by GATK HaplotypeCaller require stricter filtering due to higher false positive risk.

Criteria:
- QUAL ≥ 30
- DP ≥ 10
- GQ ≥ 20

---

### DeepVariant-Only Calls

Variants detected only by DeepVariant are retained using moderate quality thresholds reflecting DeepVariant's higher precision model.

Criteria:
- QUAL ≥ 10
- DP ≥ 10
- GQ ≥ 20

---

### Rationale

- Consensus calls allow relaxed genotype thresholds due to cross-caller support.
- Caller-specific thresholds reduce false positives while preserving sensitivity.
- DP ≥ 10 ensures sufficient local coverage while avoiding excessive loss of true variants in lower-coverage regions.

## Benchmarking Notes

- Validation performed using GIAB HG003 and HG002 samples.
- Precision metrics were calculated after ACMG region filtering.
- Performance may vary if caller combinations or filtering thresholds are modified.

---

## Tips and Best Practices

### File Path Guidelines

- Use absolute paths in the samplesheet when possible
- Ensure FASTQ files are accessible from the execution environment

### Performance Optimization

- For large cohorts, consider running samples in batches
- Use `--max_cpus` and `--max_memory` to control resource usage
- Enable `-resume` to restart from the last completed step if interrupted

### Quality Control

- Review coverage statistics in the QC directory
- Check for coverage gaps in critical regions
- Verify sex determination matches expected values
- Review Ti/Tv ratios for data quality indicators

---


## Getting Help

If you encounter issues:

1. Check the Nextflow logs
2. Review `.command.log` files in failed task directories
3. Contact your bioinformatics support team
4. Email: khatib@bionl.ai


