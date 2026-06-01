# Python Implementation Notes

## Overview

This directory contains a complete Python implementation skeleton for generating ACMG SF Variants Reports from Excel/CSV data.

## Architecture

### Core Components

1. **`generate_report.py`** - Main entry point and orchestrator
2. **`excel_reader.py`** - Data extraction from Excel/CSV files
3. **`template_processor.py`** - HTML template processing and placeholder replacement
4. **`report_generator.py`** - High-level report generation logic
5. **`utils.py`** - Utility functions for data formatting and validation
6. **`config.py`** - Configuration and mappings

### Data Flow

```
Excel/CSV Files → excel_reader.py → Structured Data
                                         ↓
Template Files ← template_processor.py ← report_generator.py
                                         ↓
                 generate_report.py → Final HTML Report
```

## Implementation Status

### ✅ Completed Components

- Complete Python skeleton with all modules
- Configuration system with mappings
- Utility functions for data formatting
- Test framework
- Documentation

### 🔧 Requires Implementation

- Actual data extraction logic in `excel_reader.py`
- Template placeholder replacement in `template_processor.py`
- Integration between components
- Error handling and validation

## Key Design Decisions

### No Third-Party Dependencies

The implementation uses only Python standard library:

- `pandas` → Use `csv` module
- `openpyxl` → Use `csv` exports from Excel
- `jinja2` → Use string replacement
- `reportlab` → Generate HTML, convert via browser

### Template Processing Strategy

Instead of template engines, the code uses:

1. Load HTML template as string
2. Define placeholder markers (e.g., `{{SAMPLE_ID}}`)
3. Replace placeholders with actual data
4. Save processed HTML

### Data Extraction Approach

1. Read CSV files exported from Excel tabs
2. Map columns using configuration
3. Apply filters and transformations
4. Structure data for template processing

## File Structure

```
python-skeleton/
├── generate_report.py      # Main entry point
├── excel_reader.py         # Data extraction
├── template_processor.py   # Template processing
├── report_generator.py     # Report orchestration
├── utils.py               # Utility functions
├── config.py              # Configuration
├── test_generation.py     # Test suite
└── output/                # Generated reports
```

## Configuration System

The `config.py` file contains:

- Column name mappings
- Filter criteria (updated to 2-star threshold)
- File naming patterns
- Formatting rules
- ACMG gene list for validation

## Testing

Run the test suite:

```bash
python test_generation.py
```

This will:

- Validate all required files exist
- Test the complete pipeline
- Generate a sample report
- Provide detailed feedback

## Integration with Nextflow

### Input Requirements

- Excel file: `{sample_id}_variants_lean_v1.xlsx`
- Or CSV files exported from Excel tabs

### Output

- HTML report: `{sample_id}_acmg_sf_report.html`
- Optional: PDF via headless browser conversion

### Nextflow Process Example

```nextflow
process GENERATE_ACMG_REPORT {
    input:
    path excel_file

    output:
    path "*.html"

    script:
    """
    python generate_report.py ${excel_file} ${excel_file.baseName}_report.html
    """
}
```

## Next Steps for Data Scientist

1. **Implement `excel_reader.py`**:

   - Read CSV files or Excel tabs
   - Map columns according to `config.py`
   - Return structured dictionaries

2. **Implement `template_processor.py`**:

   - Load HTML template
   - Replace placeholders with data
   - Handle conditional sections

3. **Test with sample data**:

   - Use provided sample CSV files
   - Run `test_generation.py`
   - Compare output with static template

4. **Integrate error handling**:

   - Missing files
   - Invalid data formats
   - Empty datasets

5. **Optimize for pipeline**:
   - Command-line argument parsing
   - Logging and monitoring
   - Performance considerations

## Questions for Data Scientist

1. **Excel vs CSV**: Do you prefer working with Excel files directly or CSV exports?
2. **Error Handling**: How should the system handle missing or invalid data?
3. **Performance**: Any specific performance requirements for the pipeline?
4. **Output Format**: Just HTML, or also PDF generation needed?
5. **Validation**: How strict should data validation be?

## Contact

For questions about the HTML template or UI requirements, contact the UI team.
For questions about the data structure or Excel files, contact the pipeline team.
