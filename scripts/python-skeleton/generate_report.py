#!/usr/bin/env python3
"""
ACMG SF Variants Report Generator

This script generates HTML reports from Excel data using the provided template.
Designed to be integrated into Nextflow pipelines.

Usage:
    python generate_report.py input.xlsx output_dir/ --sample-id HG003

Requirements:
    - Python 3.7+
    - pandas
    - openpyxl (for Excel reading)
    - No additional dependencies required
"""

import argparse
import sys
import os
from pathlib import Path
from datetime import datetime

# Import our modules
from excel_reader import ExcelDataReader
from template_processor import TemplateProcessor
from report_generator import ReportGenerator


def main():
    parser = argparse.ArgumentParser(
        description='Generate ACMG SF Variants Report from Excel data'
    )

    parser.add_argument(
        'input_excel',
        help='Path to input Excel file (e.g., HG003_variants_lean_v1.xlsx)'
    )

    parser.add_argument(
        'output_dir',
        help='Output directory for generated reports'
    )

    parser.add_argument(
        '--sample-id',
        help='Sample ID (optional, will be extracted from filename if not provided)'
    )

    parser.add_argument(
        '--assay',
        help='Assay type (e.g., WES, WGS, Panel)',
        default='WES'
    )

    parser.add_argument(
        '--template-dir',
        default='./template-files',
        help='Directory containing HTML template files'
    )

    parser.add_argument(
        '--format',
        choices=['html', 'pdf', 'both'],
        default='html',
        help='Output format(s)'
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )

    args = parser.parse_args()

    try:
        # Initialize components
        excel_reader = ExcelDataReader(args.input_excel, debug=args.debug)
        template_processor = TemplateProcessor(args.template_dir, debug=args.debug)
        report_generator = ReportGenerator(args.output_dir, debug=args.debug)

        # Extract sample ID if not provided
        sample_id = args.sample_id
        if not sample_id:
            sample_id = extract_sample_id_from_filename(args.input_excel)

        print(f"Processing sample: {sample_id}")

        # Step 1: Read and validate Excel data
        print("Reading Excel data...")
        data = excel_reader.read_all_tabs()
        excel_reader.validate_data(data)

        # Step 2: Transform data for template
        print("Processing data...")
        processed_data = excel_reader.process_for_template(data, sample_id, assay=args.assay)

        # Step 3: Generate HTML from template
        print("Generating HTML...")
        html_content = template_processor.process_template(processed_data)

        # Step 4: Save outputs
        print("Saving reports...")
        output_files = report_generator.save_reports(
            html_content,
            sample_id,
            format=args.format
        )

        print("Success!")
        for file_path in output_files:
            print(f"Generated: {file_path}")

    except Exception as e:
        print(f"Error: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def extract_sample_id_from_filename(filepath):
    """Extract sample ID from filename like 'HG003_variants_lean_v1.xlsx'"""
    filename = Path(filepath).stem
    # Assume format: {SAMPLE_ID}_variants_lean_v1
    parts = filename.split('_')
    if len(parts) >= 1:
        return parts[0]
    else:
        raise ValueError(f"Could not extract sample ID from filename: {filename}")


if __name__ == '__main__':
    main()