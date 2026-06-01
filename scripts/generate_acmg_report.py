#!/usr/bin/env python3
"""
ACMG SF Variants Report Generator - Excel Input Version
Generates HTML reports from Excel files produced by the lean report step
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime

try:
    import pandas as pd
except ImportError:
    print("Error: pandas is required. Install with: pip install pandas openpyxl")
    sys.exit(1)


class ACMGReportGenerator:
    """Complete report generator handling Excel input to HTML output"""
    
    def __init__(self, excel_path, template_path, output_dir, sample_id=None, debug=False):
        self.excel_path = Path(excel_path)
        self.template_path = Path(template_path)
        self.output_dir = Path(output_dir)
        self.sample_id = sample_id or self._extract_sample_id()
        self.debug = debug
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def _extract_sample_id(self):
        """Extract sample ID from Excel filename"""
        # Expected format: {SAMPLE_ID}_variants_lean_v1.xlsx
        filename = self.excel_path.stem
        return filename.split('_')[0] if '_' in filename else filename
    
    def _log(self, message):
        """Debug logging"""
        if self.debug:
            print(f"[DEBUG] {message}")
    
    def read_excel_data(self):
        """Read all relevant tabs from Excel file"""
        self._log(f"Reading Excel file: {self.excel_path}")
        
        try:
            excel_data = pd.read_excel(self.excel_path, sheet_name=None)
            self._log(f"Available tabs: {list(excel_data.keys())}")
            
            # Check for required tabs
            required_tabs = ['Sample Summary', 'ACMG SF (P-LP)']
            missing = [tab for tab in required_tabs if tab not in excel_data]
            if missing:
                raise ValueError(f"Missing required tabs: {missing}")
            
            return excel_data
        except Exception as e:
            raise Exception(f"Error reading Excel file: {e}")
    
    def extract_sample_info(self, data):
        """Extract sample information from Sample Summary tab"""
        df = data.get('Sample Summary')
        if df is None or len(df) == 0:
            return {}
        
        row = df.iloc[0]
        
        def safe_get(col_name, default='—', suffix=''):
            """Safely get column value with optional suffix"""
            try:
                val = row.get(col_name, default)
                if val == default or pd.isna(val):
                    return default
                return f"{float(val):.2f}{suffix}" if suffix else str(val)
            except (ValueError, TypeError):
                return default
        
        return {
            'sample_id': safe_get('Sample_ID', self.sample_id, ''),
            'assay': safe_get('Assay', 'WES', ''),
            'build': safe_get('Build', 'GRCh38', ''),
            'mean_coverage': safe_get('Mean_Coverage', '—', 'x'),
            'mapped_percent': safe_get('Mapped_Percent', '—', '%'),
            'duplicate_percent': safe_get('Duplicate_Percent', '—', '%'),
            'mean_insert_size': safe_get('Mean_Insert_Size', '—', 'bp'),
            'acmg_covered_percent': safe_get('ACMG_Covered_Percent', '—', '%'),
        }
    
    def extract_variants(self, data):
        """Extract and filter variants by star rating"""
        df = data.get('ACMG SF (P-LP)')
        if df is None or len(df) == 0:
            return {'primary': [], 'secondary': []}
        
        # Filter for pathogenic/likely pathogenic
        pathogenic_mask = df['ClinVar'].isin(['Pathogenic', 'Likely pathogenic'])
        variants_df = df[pathogenic_mask].copy()
        
        primary_variants = []
        secondary_variants = []
        
        for _, row in variants_df.iterrows():
            variant = {
                'gene': str(row.get('Gene', '—')),
                'hgvsc': str(row.get('HGVSc', '—')),
                'hgvsp': str(row.get('HGVSp', '—')),
            }
            
            # Get star rating
            try:
                stars = int(row.get('ClinVar_Stars', 0)) if pd.notna(row.get('ClinVar_Stars')) else 0
            except (ValueError, TypeError):
                stars = 0
            
            # Split by star rating: >=2 stars = primary, <2 stars = secondary
            if stars >= 2:
                primary_variants.append(variant)
            else:
                secondary_variants.append(variant)
        
        self._log(f"Found {len(primary_variants)} primary variants (>=2 stars)")
        self._log(f"Found {len(secondary_variants)} secondary variants (<2 stars)")
        
        return {'primary': primary_variants, 'secondary': secondary_variants}
    
    def extract_coverage(self, data):
        """Extract coverage data for genes"""
        df = data.get('ACMG Genes Coverage')
        if df is None or len(df) == 0:
            return {'all': [], 'gaps': []}
        
        all_coverage = []
        gaps = []
        
        for _, row in df.iterrows():
            # Try different possible column names for transcript ID
            transcript = (row.get('MANE_ID') or 
                         row.get('ENSG_ENST_ID') or 
                         row.get('Transcript_ID') or '—')
            
            try:
                coverage_val = float(row.get('Pct>=20x', 100))
            except (ValueError, TypeError):
                coverage_val = 100.0
            
            cov_item = {
                'gene': str(row.get('Gene', '—')),
                'transcript': str(transcript),
                'coverage': f"{coverage_val:.0f}%"
            }
            
            all_coverage.append(cov_item)
            
            # Track genes with <100% coverage
            if coverage_val < 100:
                gaps.append(cov_item)
        
        self._log(f"Found {len(all_coverage)} genes with coverage data")
        self._log(f"Found {len(gaps)} genes with coverage gaps (<100%)")
        
        return {'all': all_coverage, 'gaps': gaps}
    
    def generate_variant_rows(self, variants):
        """Generate HTML table rows for variants"""
        if not variants:
            return '              <tr><td colspan="3" class="text-center">No variants found</td></tr>'
        
        rows = []
        for v in variants:
            row = f'''              <tr>
                <td class="font-semibold">{v['gene']}</td>
                <td class="font-mono text-sm">{v['hgvsc']}</td>
                <td class="font-mono text-sm">{v['hgvsp']}</td>
              </tr>'''
            rows.append(row)
        
        return '\n'.join(rows)
    
    def generate_coverage_grid(self, coverage_items):
        """Generate HTML for coverage grid"""
        if not coverage_items:
            return '          <div class="coverage-item"><span class="gene-name">No coverage data available</span></div>'
        
        items = []
        for item in coverage_items:
            html = f'''          <div class="coverage-item">
            <span class="gene-name">{item['gene']}</span>
            <span class="transcript">{item['transcript']}</span>
            <span class="coverage">{item['coverage']}</span>
          </div>'''
            items.append(html)
        
        return '\n'.join(items)
    
    def generate_coverage_gaps_section(self, gaps):
        """Generate HTML for coverage gaps section"""
        if not gaps:
            return '<p class="no-data">No data available.</p>'
        
        # Simple list format
        items = [f"{g['gene']} ({g['coverage']})" for g in gaps]
        return '<div class="coverage-gaps">\n  <ul>\n' + '\n'.join(
            f'    <li>{item}</li>' for item in items
        ) + '\n  </ul>\n</div>'
    
    def generate_html_report(self, excel_data):
        """Generate complete HTML report"""
        self._log("Processing Excel data...")
        
        # Extract all data sections
        sample_info = self.extract_sample_info(excel_data)
        variants = self.extract_variants(excel_data)
        coverage = self.extract_coverage(excel_data)
        
        # Read template
        self._log(f"Reading template: {self.template_path}")
        with open(self.template_path, 'r', encoding='utf-8') as f:
            html = f.read()
        
        # Generate timestamp
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M')
        
        # Replace header information
        replacements = {
            'HG003': sample_info.get('sample_id', self.sample_id),
            'WES': sample_info.get('assay', 'WES'),
            'GRCh38': sample_info.get('build', 'GRCh38'),
            '2025-09-04 08:58': timestamp,
        }
        
        for old, new in replacements.items():
            html = html.replace(old, new)
        
        # Replace QC metrics
        qc_replacements = {
            '62.08x': sample_info.get('mean_coverage', '—'),
            '99.97%': sample_info.get('mapped_percent', '—'),
            '400.1bp': sample_info.get('mean_insert_size', '—'),
            '99.57%': sample_info.get('acmg_covered_percent', '—'),
        }
        
        for old, new in qc_replacements.items():
            html = html.replace(old, new)
        
        # Handle duplicate percent (special case with "—" default)
        dup_pct = sample_info.get('duplicate_percent', '—')
        html = html.replace(
            '<td class="metric-value">—</td>',
            f'<td class="metric-value">{dup_pct}</td>',
            1  # Only replace first occurrence
        )
        
        # Replace variant tables using multiple tbody approach
        import re
        
        # Generate variant rows
        primary_rows = self.generate_variant_rows(variants['primary'])
        secondary_rows = self.generate_variant_rows(variants['secondary'])
        
        # Find all tbody tags
        tbody_pattern = r'(<tbody>)(.*?)(</tbody>)'
        matches = list(re.finditer(tbody_pattern, html, re.DOTALL))
        
        if len(matches) >= 2:
            # Replace first tbody (Page 1 - primary variants)
            html = html[:matches[0].start(2)] + f'\n{primary_rows}\n            ' + html[matches[0].end(2):]
            
            # Re-find matches after first replacement
            matches = list(re.finditer(tbody_pattern, html, re.DOTALL))
            
            # Replace second tbody (Page 2 - secondary variants)
            if len(matches) >= 2:
                html = html[:matches[1].start(2)] + f'\n{secondary_rows}\n            ' + html[matches[1].end(2):]
        
        # Replace coverage gaps section
        gaps_html = self.generate_coverage_gaps_section(coverage['gaps'])
        html = html.replace(
            '<p class="no-data">No data available.</p>',
            gaps_html,
            1  # Only first occurrence
        )
        
        # Replace coverage grid
        coverage_grid_html = self.generate_coverage_grid(coverage['all'])
        grid_pattern = r'(<div class="coverage-grid">)(.*?)(</div>)'
        html = re.sub(
            grid_pattern,
            f'\\1\n{coverage_grid_html}\n        \\3',
            html,
            count=1,
            flags=re.DOTALL
        )
        
        return html
    
    def save_report(self, html_content):
        """Save HTML report to file"""
        output_file = self.output_dir / f"{self.sample_id}_acmg_report.html"
        
        self._log(f"Writing report to: {output_file}")
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        return output_file
    
    def generate(self):
        """Main generation workflow"""
        print(f"Processing sample: {self.sample_id}")
        
        # Read Excel data
        excel_data = self.read_excel_data()
        
        # Generate HTML
        html_content = self.generate_html_report(excel_data)
        
        # Save report
        output_file = self.save_report(html_content)
        
        print(f"Report generated: {output_file}")
        return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Generate ACMG SF Variants HTML Report from Excel'
    )
    parser.add_argument('excel_file', help='Input Excel file from lean report step')
    parser.add_argument('template_file', help='HTML template file (index.html)')
    parser.add_argument('output_dir', help='Output directory for report')
    parser.add_argument('--sample-id', help='Sample ID (auto-detected if not provided)')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    args = parser.parse_args()
    
    try:
        generator = ACMGReportGenerator(
            excel_path=args.excel_file,
            template_path=args.template_file,
            output_dir=args.output_dir,
            sample_id=args.sample_id,
            debug=args.debug
        )
        
        generator.generate()
        print("Success!")
        
    except Exception as e:
        print(f"Error: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()