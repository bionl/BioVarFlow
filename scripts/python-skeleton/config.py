# ACMG SF Variants Report - Configuration
# This file contains all configurable parameters for the report generation

# Template Configuration
TEMPLATE_CONFIG = {
    "template_file": "index.html",
    "css_file": "style.css",
    "output_encoding": "utf-8"
}

# Excel/CSV Column Mappings
# These should match the actual column names in your Excel/CSV files
COLUMN_MAPPINGS = {
    "sample_summary": {
        "sample_id": "Sample_ID",
        "assay": "Assay",
        "build": "Build",
        "generated": "generated",
        "mean_coverage": "Mean_Coverage",
        "mapped_percent": "Mapped_Percent",
        "duplicate_percent": "Duplicate_Percent",
        "mean_insert_size": "Mean_Insert_Size",
        "acmg_covered_percent": "ACMG_PctRegionsCovered"
    },
    "variants": {
        "gene": "Gene",
        "hgvsc": "HGVSc",
        "hgvsp": "HGVSp",
        "clinvar": "ClinVar",
        "clinvar_stars": "ClinVar_Stars"
    },
    "coverage": {
        "gene": "Gene",
        "transcript": "MANE_ID",
        "coverage_20x": "Pct>=20x"
    }
}

# Filter Criteria
FILTER_CRITERIA = {
    "primary_variants": {
        "clinvar_pathogenic": ["Pathogenic", "Likely pathogenic"],
        "min_stars": 2  # Updated from 3 to 2 stars
    },
    "secondary_variants": {
        "clinvar_pathogenic": ["Pathogenic", "Likely pathogenic"],
        "max_stars": 1  # Less than 2 stars
    },
    "coverage_gaps": {
        "min_coverage_threshold": 100  # Show genes with <100% coverage at 20x
    }
}

# File Name Patterns
FILE_PATTERNS = {
    "excel_pattern": "{sample_id}_variants_lean_v1.xlsx",
    "csv_patterns": {
        "sample_summary": "Sample Summary-Table 1.csv",
        "variants": "ACMG SF (P-LP)-Table 1.csv",
        "coverage": "ACMG Genes Coverage-Table 1.csv",
        "coverage_gaps": "Coverage gaps-Table 1.csv",
        "pass_variants": "PASS variants-Table 1.csv"
    }
}

# Data Formatting Rules
FORMATTING_RULES = {
    "coverage_suffix": "x",
    "percentage_suffix": "%",
    "size_suffix": "bp",
    "date_format": "%Y-%m-%d %H:%M",
    "missing_data_placeholder": "—",
    "no_data_message": "No data available."
}

# ACMG SF v3.3 Gene List (for validation)
ACMG_GENES = [
    "SDHB", "MUTYH", "PCSK9", "RPE65", "CASQ2", "LMNA", "SDHC", "CACNA1S",
    "TNNT2", "RYR2", "RET", "BMPR1A", "ACTA2", "RBM20", "BAG3", "KCNQ1",
    "WT1", "MYBPC3", "SDHAF2", "MEN1", "SDHD", "PKP2", "ACVRL1", "MYL2",
    "HNF1A", "BRCA2", "RB1", "ATP7B", "MYH7", "MAX", "CALM1", "ACTC1",
    "FBN1", "TPM1", "SMAD3", "TSC2", "MYH11", "PALB2", "TP53", "BRCA1",
    "GAA", "DSC2", "DSG2", "TTR", "SMAD4", "STK11", "LDLR", "RYR1",
    "CALM3", "TNNI3", "APOB", "CALM2", "MSH2", "MSH6", "TMEM127", "TTN",
    "COL3A1", "CYP27A1", "DES", "NF2", "VHL", "TMEM43", "BTD", "TGFBR2",
    "MLH1", "SCN5A", "MYL3", "TNNC1", "APC", "DSP", "HFE", "PLN",
    "TRDN", "PMS2", "FLNC", "KCNH2", "PRKAG2", "TGFBR1", "ENG", "TSC1",
    "OTC", "GLA", "ABCD1"
]

# Pipeline Information
PIPELINE_INFO = {
    "default_pipeline": "Bionl_Lean_call v1.0",
    "default_databases": "ClinVar (2025-01), gnomAD v4.1, Ensembl VEP Release 115, REVEL (latest release), AlphaMissense (Science 2023, updated 2025-05)"
}

# Debug Settings
DEBUG_CONFIG = {
    "verbose_logging": True,
    "save_intermediate_data": False,
    "validate_gene_list": True
}