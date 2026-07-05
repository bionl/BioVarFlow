"""
Template Processor Module

Handles processing of HTML templates and data substitution.
Uses string replacement instead of external templating libraries.
"""

import os
import re
from typing import Dict, Any, List
from pathlib import Path


class TemplateProcessor:
    """Processes HTML templates with data substitution"""

    def __init__(self, template_dir: str, debug: bool = False):
        self.template_dir = Path(template_dir)
        self.debug = debug

        # Template file paths
        self.html_template = self.template_dir / 'index.html'
        self.css_template = self.template_dir / 'style.css'

        if not self.html_template.exists():
            raise FileNotFoundError(f"HTML template not found: {self.html_template}")

    def process_template(self, data: Dict[str, Any]) -> str:
        """Process the HTML template with the provided data"""

        # Read the template
        with open(self.html_template, 'r', encoding='utf-8') as f:
            html_content = f.read()

        # Read CSS if it exists
        css_content = ""
        if self.css_template.exists():
            with open(self.css_template, 'r', encoding='utf-8') as f:
                css_content = f.read()

        # Replace CSS link with inline styles for standalone HTML
        if css_content:
            css_inline = f'<style>\n{css_content}\n</style>'
            html_content = re.sub(
                r'<link[^>]*href=["\']style\.css["\'][^>]*>',
                css_inline,
                html_content
            )

        # Process data substitutions
        html_content = self._substitute_header_data(html_content, data)
        html_content = self._substitute_variant_tables(html_content, data)
        html_content = self._substitute_qc_metrics(html_content, data)
        html_content = self._substitute_coverage_data(html_content, data)
        html_content = self._substitute_methods_limitations(html_content, data)

        return html_content

    def _substitute_header_data(self, html: str, data: Dict[str, Any]) -> str:
        """Replace header information in the template"""

        # Sample ID
        sample_id = data.get('sample_id', 'UNKNOWN')
        assay = (data.get("assay") or "WES").strip()   # expected: WES/WGS/NA
        reference_build = data.get("reference_build", "GRCh38")

        # 1) Replace sample id 
        html = html.replace("HG003", sample_id)

        # 2) Build assay labels
        assay_short_map = {"WES": "WES", "WGS": "WGS", "NA": "NA"}
        assay_long_map = {"WES": "Whole Exome Sequencing (WES)", "WGS": "Whole Genome Sequencing (WGS)", "NA":  "Not specified"}
        assay_short = assay_short_map.get(assay, assay)
        assay_long  = assay_long_map.get(assay, assay)

        # 3) Replace LONG assay phrase first (prevents "Whole Exome..." staying for WGS)
        html = html.replace("Whole Exome Sequencing (WES)", assay_long)
        html = html.replace("Whole Genome Sequencing (WGS)", assay_long)
        html = re.sub(r"\bWES\b", assay_short, html)
        html = re.sub(r"\bWGS\b", assay_short, html)


        # Other header fields
        replacements = {
            'GRCh38': reference_build,
            '2025-09-04 08:58': data.get('report_generated', ''),
            'Sarek 3.5.1 - ClinLEAN Reporting Workflow v1': data.get('pipeline', 'BioVarFlow v1.1.0'),
            'ClinVar, gnomAD (v4.1), VEP/Ensembl (Release 115)': data.get('databases', 'ClinVar (2025-01), gnomAD v4.1, Ensembl VEP Release 115, REVEL (latest release), AlphaMissense (Science 2023, updated 2025-05) SpliceAI (version 1.3.1), BayesDel (version 1.0)'),
        }

        for old_value, new_value in replacements.items():
            html = html.replace(old_value, str(new_value))

        return html

    def _substitute_variant_tables(self, html: str, data: Dict[str, Any]) -> str:
        """Replace variant tables for Page 1 and 2"""

        # Page 1 variants (high confidence)
        page1_variants = data.get('page1_variants', [])
        page1_html = self._generate_variant_rows(page1_variants)

        # Find and replace Page 1 table body
        page1_pattern = r'(<tbody>\s*<tr>\s*<td class="font-semibold">HNF1A</td>.*?</tr>\s*</tbody>)'
        if page1_html:
            html = re.sub(page1_pattern, f'<tbody>\n{page1_html}\n</tbody>', html, flags=re.DOTALL)
        else:
            # No variants found
            no_data_row = '<tr><td colspan="3" class="text-center">No variants meeting these criteria were identified in this sample.</td></tr>'
            html = re.sub(page1_pattern, f'<tbody>\n{no_data_row}\n</tbody>', html, flags=re.DOTALL)

        # Page 2 variants (lower confidence)
        page2_variants = data.get('page2_variants', [])
        page2_html = self._generate_variant_rows(page2_variants)

        # Find and replace Page 2 table body (has TP53 and GAA)
        page2_pattern = r'(<tbody>\s*<tr>\s*<td class="font-semibold">TP53</td>.*?</tr>\s*</tbody>)'
        if page2_html:
            html = re.sub(page2_pattern, f'<tbody>\n{page2_html}\n</tbody>', html, flags=re.DOTALL)
        else:
            # No variants found
            no_data_row = '<tr><td colspan="3" class="text-center">No variants meeting these criteria were identified in this sample.</td></tr>'
            html = re.sub(page2_pattern, f'<tbody>\n{no_data_row}\n</tbody>', html, flags=re.DOTALL)

        return html

    def _generate_variant_rows(self, variants: List[Dict[str, Any]]) -> str:
        """Generate HTML table rows for variants"""
        if not variants:
            return ""

        rows = []
        for variant in variants:
            gene = variant.get('gene', '—')
            hgvsc = variant.get('hgvsc', '—')
            hgvsp = variant.get('hgvsp', '—')
            clinvar_id = variant.get('clinvar_id', '')

            for key in ['hgvsc', 'hgvsp']:
                val = variant.get(key)
                if val is None or str(val).lower() == 'nan' or str(val).strip() == '':
                    variant[key] = '-'

            hgvsc = variant['hgvsc']
            hgvsp = variant['hgvsp']

            # Create gene cell with optional ClinVar link
            if clinvar_id:
                gene_cell = f'<td class="font-semibold"><a href="{clinvar_id}" target="_blank">{gene}</a></td>'
            else:
                gene_cell = f'<td class="font-semibold">{gene}</td>'

            row = f'''        <tr>
          {gene_cell}
          <td class="font-mono text-sm">{hgvsc}</td>
          <td class="font-mono text-sm">{hgvsp}</td>
        </tr>'''

            rows.append(row)

        return '\n'.join(rows)

    def _substitute_qc_metrics(self, html: str, data: Dict[str, Any]) -> str:
        """Replace QC metrics in the template"""

        # Static values to replace with dynamic data
        replacements = {
            '62.08x': data.get('mean_coverage', '—'),
            '99.97%': data.get('mapped_percent', '—'),
            '400.1bp': data.get('mean_insert_size', '—'),
            '99.57%': data.get('acmg_covered_percent', '—'),
        }

        # Handle duplicate percent (shows as "—" in template)
        #duplicate_percent = data.get('duplicate_percent', '—')
        #if duplicate_percent != '—':
        #    html = html.replace('<td class="metric-value">—</td>', f'<td class="metric-value">{duplicate_percent}</td>', 1)

        for old_value, new_value in replacements.items():
            html = html.replace(old_value, str(new_value))

        return html

    def _substitute_coverage_data(self, html: str, data: Dict[str, Any]) -> str:
        """Replace coverage data for both sections on Page 3"""

        # Coverage gaps section
        coverage_gaps = data.get('coverage_gaps', [])
        print(coverage_gaps)
        if coverage_gaps:
            # Generate coverage gaps content
            gaps_html = self._generate_coverage_gaps_html(coverage_gaps)
            # Replace the "No data available" message
            html = html.replace(
                '<p class="no-data">No exons below 20x coverage.</p>',
                gaps_html
            )

        # Overall coverage section
        overall_coverage = data.get('overall_coverage', [])
        if overall_coverage:
            coverage_grid_html = self._generate_coverage_grid_html(overall_coverage)
            # Replace the coverage grid content
            grid_pattern = r'(<div class="coverage-grid">.*?</div>)'
            html = re.sub(
                grid_pattern,
                f'<div class="coverage-grid">\n{coverage_grid_html}\n        </div>',
                html,
                flags=re.DOTALL
            )

        return html

    def _generate_coverage_gaps_html(self, gaps: List[Dict[str, Any]]) -> str:
        """Generate HTML for coverage gaps section"""
        if not gaps:
            return '<p class="no-data">No exons below 20x coverage.</p>'

        # Create a simple table for coverage gaps
        html = '<table class="coverage-gaps-table">\n'
        #html += '  <thead>\n    <tr>\n      <th>Gene</th>\n      <th>Transcript</th>\n      <th>Coverage</th>\n    </tr>\n  </thead>\n'
        html += '  <thead>\n    <tr>\n      <th>Gene</th>\n      <th>ExonStart</th>\n      <th>ExonEnd</th>\n      <th>Coverage</th>\n    </tr>\n  </thead>\n'
        html += '  <tbody>\n'

        for gap in gaps:
            html += f'    <tr>\n'
            html += f'      <td>{gap.get("gene", "—")}</td>\n'
            #html += f'      <td>{gap.get("transcript", "—")}</td>\n'
            html += f'      <td>{gap.get("ExonStart", "—")}</td>\n'
            html += f'      <td>{gap.get("ExonEnd", "—")}</td>\n'
            html += f'      <td>{gap.get("coverage", "—")}</td>\n'
            html += f'    </tr>\n'

        html += '  </tbody>\n</table>'
        return html

    def _generate_coverage_grid_html(self, coverage_data: List[Dict[str, Any]]) -> str:
        """Generate HTML for the overall coverage grid"""
        if not coverage_data:
            return ""

        items = []
        for item in coverage_data:
            gene = item.get('gene', '—')
            transcript = item.get('transcript', '—')
            coverage = item.get('coverage', '—')

            item_html = f'''          <div class="coverage-item">
            <span class="gene-name">{gene}</span>
            <span class="transcript">{transcript}</span>
            <span class="coverage">{coverage}</span>
          </div>'''
            print(item_html)
            items.append(item_html)

        return '\n'.join(items)

    def _substitute_methods_limitations(self, html: str, data: Dict[str, Any]) -> str:
        """Replace methods and limitations in the template"""
        # Sample ID
        #sample_id = data.get('sample_id', 'UNKNOWN')
        #html = html.replace('HG003', 'sample_id')

        # Methods and limitations
        replacements = {
            '2025-10-05 08:58': data.get('report_generated', ''),
            'BioVarFlow v1.1.0': data.get('pipeline', 'BioVarFlow v1.1.0'),
        }

        for old_value, new_value in replacements.items():
            html = html.replace(old_value, str(new_value))
        return html