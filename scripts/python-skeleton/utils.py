"""
Utility functions for data processing and validation.
These helper functions support the main report generation pipeline.
"""

import re
from datetime import datetime
from typing import Any, List, Dict, Optional
from config import FORMATTING_RULES, ACMG_GENES

def format_coverage(value: Any) -> str:
    """
    Format coverage values with 'x' suffix.

    Args:
        value: Coverage value (number or string)

    Returns:
        Formatted coverage string (e.g., "62.08x")
    """
    if value is None or str(value).strip() == "":
        return FORMATTING_RULES["missing_data_placeholder"]

    try:
        # Convert to float and format
        coverage_num = float(str(value).replace('x', '').strip())
        return f"{coverage_num:.2f}x"
    except (ValueError, TypeError):
        return str(value)

def format_percentage(value: Any) -> str:
    """
    Format percentage values with '%' suffix.

    Args:
        value: Percentage value (number or string)

    Returns:
        Formatted percentage string (e.g., "99.97%")
    """
    if value is None or str(value).strip() == "":
        return FORMATTING_RULES["missing_data_placeholder"]

    try:
        # Convert to float and format
        pct_num = float(str(value).replace('%', '').strip())
        return f"{pct_num:.2f}%"
    except (ValueError, TypeError):
        return str(value)

def format_size(value: Any) -> str:
    """
    Format size values with 'bp' suffix.

    Args:
        value: Size value (number or string)

    Returns:
        Formatted size string (e.g., "400.1bp")
    """
    if value is None or str(value).strip() == "":
        return FORMATTING_RULES["missing_data_placeholder"]

    try:
        # Convert to float and format
        size_num = float(str(value).replace('bp', '').strip())
        return f"{size_num:.1f}bp"
    except (ValueError, TypeError):
        return str(value)

def format_date(value: Any) -> str:
    """
    Format date values consistently.

    Args:
        value: Date value (string or datetime)

    Returns:
        Formatted date string (e.g., "2025-09-04 08:58")
    """
    if value is None or str(value).strip() == "":
        return FORMATTING_RULES["missing_data_placeholder"]

    try:
        if isinstance(value, datetime):
            return value.strftime(FORMATTING_RULES["date_format"])

        # Try to parse string dates
        date_str = str(value).strip()

        # Common date patterns
        patterns = [
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%d %H:%M",
            "%Y-%m-%d",
            "%m/%d/%Y %H:%M",
            "%m/%d/%Y"
        ]

        for pattern in patterns:
            try:
                parsed_date = datetime.strptime(date_str, pattern)
                return parsed_date.strftime(FORMATTING_RULES["date_format"])
            except ValueError:
                continue

        # If no pattern matches, return as-is
        return date_str

    except Exception:
        return str(value)

def clean_gene_name(gene: str) -> str:
    """
    Clean and validate gene names.

    Args:
        gene: Gene name string

    Returns:
        Cleaned gene name
    """
    if not gene or str(gene).strip() == "":
        return FORMATTING_RULES["missing_data_placeholder"]

    # Remove extra whitespace and convert to uppercase
    cleaned = str(gene).strip().upper()

    # Remove any non-alphanumeric characters except hyphens
    cleaned = re.sub(r'[^A-Z0-9\-]', '', cleaned)

    return cleaned

def validate_gene_list(genes: List[str]) -> Dict[str, List[str]]:
    """
    Validate gene list against ACMG SF v3.3 genes.

    Args:
        genes: List of gene names to validate

    Returns:
        Dictionary with 'valid', 'invalid', and 'missing' gene lists
    """
    cleaned_genes = [clean_gene_name(gene) for gene in genes if gene]
    acmg_set = set(ACMG_GENES)
    input_set = set(cleaned_genes)

    return {
        'valid': sorted(list(input_set.intersection(acmg_set))),
        'invalid': sorted(list(input_set - acmg_set)),
        'missing': sorted(list(acmg_set - input_set))
    }

def safe_get_value(data: Dict, key: str, default: str = None) -> str:
    """
    Safely get value from dictionary with fallback.

    Args:
        data: Dictionary to search
        key: Key to look for
        default: Default value if key not found

    Returns:
        Value or default/placeholder
    """
    if default is None:
        default = FORMATTING_RULES["missing_data_placeholder"]

    value = data.get(key, default)

    if value is None or str(value).strip() == "":
        return default

    return str(value).strip()

def filter_variants_by_criteria(variants: List[Dict], criteria: Dict) -> List[Dict]:
    """
    Filter variants based on criteria (pathogenicity and stars).

    Args:
        variants: List of variant dictionaries
        criteria: Filter criteria dictionary

    Returns:
        Filtered list of variants
    """
    filtered = []

    for variant in variants:
        # Check ClinVar pathogenicity
        clinvar = variant.get('ClinVar', '').strip()
        if clinvar not in criteria.get('clinvar_pathogenic', []):
            continue

        # Check star rating if specified
        if 'min_stars' in criteria:
            try:
                stars = float(variant.get('ClinVar_Stars', 0))
                if stars < criteria['min_stars']:
                    continue
            except (ValueError, TypeError):
                continue

        if 'max_stars' in criteria:
            try:
                stars = float(variant.get('ClinVar_Stars', 0))
                if stars > criteria['max_stars']:
                    continue
            except (ValueError, TypeError):
                # If we can't parse stars, include in secondary analysis
                pass

        filtered.append(variant)

    return filtered

def generate_coverage_grid_html(coverage_data: List[Dict]) -> str:
    """
    Generate HTML for the coverage grid section.

    Args:
        coverage_data: List of coverage dictionaries

    Returns:
        HTML string for coverage grid
    """
    if not coverage_data:
        return f'<p class="no-data">{FORMATTING_RULES["no_data_message"]}</p>'

    html_items = []

    for item in coverage_data:
        gene = item.get('Gene', FORMATTING_RULES["missing_data_placeholder"])
        transcript = item.get('ENSG_ENST_ID', FORMATTING_RULES["missing_data_placeholder"])
        coverage = format_percentage(item.get('Pct>=20x', ''))

        html_items.append(f'''
          <div class="coverage-item">
            <span class="gene-name">{gene}</span>
            <span class="transcript">{transcript}</span>
            <span class="coverage">{coverage}</span>
          </div>''')

    return ''.join(html_items)

def generate_variants_table_html(variants: List[Dict]) -> str:
    """
    Generate HTML for variants table rows.

    Args:
        variants: List of variant dictionaries

    Returns:
        HTML string for table rows
    """
    if not variants:
        return '''
            <tr>
                <td colspan="3" class="no-data">No data available.</td>
            </tr>'''

    html_rows = []

    for variant in variants:
        gene = variant.get('Gene', FORMATTING_RULES["missing_data_placeholder"])
        hgvsc = variant.get('HGVSc', FORMATTING_RULES["missing_data_placeholder"])
        hgvsp = variant.get('HGVSp', FORMATTING_RULES["missing_data_placeholder"])

        html_rows.append(f'''
            <tr>
                <td class="font-semibold">{gene}</td>
                <td class="font-mono text-sm">{hgvsc}</td>
                <td class="font-mono text-sm">{hgvsp}</td>
            </tr>''')

    return ''.join(html_rows)

def log_processing_info(message: str, data: Any = None) -> None:
    """
    Log processing information for debugging.

    Args:
        message: Log message
        data: Optional data to include
    """
    from config import DEBUG_CONFIG

    if DEBUG_CONFIG.get('verbose_logging', False):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")

        if data is not None:
            print(f"           Data: {data}")