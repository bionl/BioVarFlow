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
        html_content = self._substitute_hemonc_section(html_content, data)
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
        # NOTE: count=1 — the HemOnc page has its own identical fallback
        # sentence, filled separately by _substitute_hemonc_section. Without
        # this limit, ACMG gaps content would bleed into the HemOnc page.
        coverage_gaps = data.get('coverage_gaps', [])
        print(coverage_gaps)
        if coverage_gaps:
            gaps_html = self._generate_coverage_gaps_html(coverage_gaps)
            html = html.replace(
                '<p class="no-data">No exons below 20x coverage.</p>',
                gaps_html,
                1,
            )

        # Overall coverage section
        # NOTE: count=1 — same reason as above; the HemOnc page also has a
        # <div class="coverage-grid">…</div> block, filled by
        # _substitute_hemonc_section using its own marker regions.
        overall_coverage = data.get('overall_coverage', [])
        if overall_coverage:
            coverage_grid_html = self._generate_coverage_grid_html(overall_coverage)
            grid_pattern = r'(<div class="coverage-grid">.*?</div>)'
            html = re.sub(
                grid_pattern,
                f'<div class="coverage-grid">\n{coverage_grid_html}\n        </div>',
                html,
                count=1,
                flags=re.DOTALL,
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

    # ---------------------------------------------------------------
    # HemOnc pages (BioVarFlow_HemOnc branch)
    # ---------------------------------------------------------------
    def _substitute_hemonc_section(self, html: str, data: Dict[str, Any]) -> str:
        """
        Fill the four marker regions on the three HemOnc pages of the template:
          - HEMONC_VARIANTS_HIGH_{START,END}   → variants with ClinVar ≥2 stars
          - HEMONC_VARIANTS_LOW_{START,END}    → variants with ClinVar <2 stars / VUS
          - HEMONC_COVERAGE_GAPS_{START,END}   → exons with %20x < 100
          - HEMONC_COVERAGE_GRID_{START,END}   → per-gene overall coverage grid
          - HEMONC_COVERED_PCT_{START,END}     → Sample Summary metric cell

        Each marker's default content (a "no data" placeholder) is preserved
        when the corresponding data list is empty, so ACMG-only runs render
        cleanly.
        """

        def replace_between(html: str, start_marker: str, end_marker: str, replacement: str) -> str:
            """Replace everything between two HTML-comment markers with `replacement`.
            Markers themselves are preserved so the block can be re-substituted
            during re-rendering / debugging."""
            pattern = re.compile(
                re.escape(f"<!-- {start_marker} -->") + r".*?" + re.escape(f"<!-- {end_marker} -->"),
                flags=re.DOTALL,
            )
            wrapped = f"<!-- {start_marker} -->\n{replacement}\n              <!-- {end_marker} -->"
            # Use a lambda so any backreference-like sequences in `wrapped`
            # (e.g. \1, \g<...>) are treated literally, not as replacement
            # references. Matters for URLs and gene symbols with backslashes.
            return pattern.sub(lambda _m: wrapped, html, count=1)

        # ---- High-confidence variant table ----
        page1 = data.get('hemonc_page1_variants', [])
        rows = self._generate_variant_rows(page1)
        if not rows:
            rows = '<tr><td colspan="3" class="text-center">No variants meeting these criteria were identified in this sample.</td></tr>'
        html = replace_between(html, 'HEMONC_VARIANTS_HIGH_START', 'HEMONC_VARIANTS_HIGH_END', rows)

        # ---- Lower-confidence variant table ----
        page2 = data.get('hemonc_page2_variants', [])
        rows = self._generate_variant_rows(page2)
        if not rows:
            rows = '<tr><td colspan="3" class="text-center">No variants meeting these criteria were identified in this sample.</td></tr>'
        html = replace_between(html, 'HEMONC_VARIANTS_LOW_START', 'HEMONC_VARIANTS_LOW_END', rows)

        # ---- Coverage gaps list ----
        gaps = data.get('hemonc_coverage_gaps', [])
        gaps_html = (self._generate_coverage_gaps_html(gaps)
                     if gaps
                     else '<p class="no-data">No exons below 20x coverage.</p>')
        html = replace_between(html, 'HEMONC_COVERAGE_GAPS_START', 'HEMONC_COVERAGE_GAPS_END', gaps_html)

        # ---- Overall coverage grid ----
        grid = data.get('hemonc_overall_coverage', [])
        grid_html = self._generate_coverage_grid_html(grid) if grid else ''
        html = replace_between(html, 'HEMONC_COVERAGE_GRID_START', 'HEMONC_COVERAGE_GRID_END', grid_html)

        # ---- HemOnc-scoped Sample Summary metric cell ----
        pct = data.get('hemonc_covered_percent', '—')
        cell_html = f'<td class="metric-value">{pct}</td>'
        html = replace_between(html, 'HEMONC_COVERED_PCT_START', 'HEMONC_COVERED_PCT_END', cell_html)

        return html

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