"""
Report Generator Module

Handles saving HTML reports and optional PDF generation.
"""

import os
import subprocess
from pathlib import Path
from typing import List


class ReportGenerator:
    """Generates and saves reports in various formats"""

    def __init__(self, output_dir: str, debug: bool = False):
        self.output_dir = Path(output_dir)
        self.debug = debug

        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def save_reports(self, html_content: str, sample_id: str, format: str = 'html') -> List[str]:
        """Save reports in the specified format(s)"""

        output_files = []

        # Generate HTML
        if format in ['html', 'both']:
            html_path = self._save_html(html_content, sample_id)
            output_files.append(html_path)

        # Generate PDF (if requested)
        if format in ['pdf', 'both']:
            pdf_path = self._save_pdf(html_content, sample_id)
            if pdf_path:
                output_files.append(pdf_path)

        return output_files

    def _save_html(self, html_content: str, sample_id: str) -> str:
        """Save HTML report"""
        filename = f"{sample_id}_clinical_report.html"
        filepath = self.output_dir / filename

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(html_content)

        if self.debug:
            print(f"HTML saved to: {filepath}")

        return str(filepath)

    def _save_pdf(self, html_content: str, sample_id: str) -> str:
        """Save PDF report using system tools"""

        # First save HTML temporarily
        temp_html = self.output_dir / f"{sample_id}_temp.html"
        with open(temp_html, 'w', encoding='utf-8') as f:
            f.write(html_content)

        # PDF output path
        pdf_filename = f"{sample_id}_clinical_report.pdf"
        pdf_path = self.output_dir / pdf_filename

        try:
            # Try different PDF generation methods
            success = (
                self._try_wkhtmltopdf(temp_html, pdf_path) or
                self._try_weasyprint(temp_html, pdf_path) or
                self._try_chrome_headless(temp_html, pdf_path)
            )

            if success:
                if self.debug:
                    print(f"PDF saved to: {pdf_path}")
                return str(pdf_path)
            else:
                print("Warning: Could not generate PDF. No suitable PDF generator found.")
                print("Consider installing: wkhtmltopdf, weasyprint, or Google Chrome/Chromium")
                return None

        except Exception as e:
            print(f"Warning: PDF generation failed: {e}")
            return None

        finally:
            # Clean up temporary HTML
            if temp_html.exists():
                temp_html.unlink()

    def _try_wkhtmltopdf(self, html_path: Path, pdf_path: Path) -> bool:
        """Try using wkhtmltopdf for PDF generation"""
        try:
            subprocess.run([
                'wkhtmltopdf',
                '--page-size', 'A4',
                '--margin-top', '20mm',
                '--margin-bottom', '20mm',
                '--margin-left', '15mm',
                '--margin-right', '15mm',
                '--print-media-type',
                str(html_path),
                str(pdf_path)
            ], check=True, capture_output=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def _try_weasyprint(self, html_path: Path, pdf_path: Path) -> bool:
        """Try using WeasyPrint for PDF generation"""
        try:
            import weasyprint
            weasyprint.HTML(filename=str(html_path)).write_pdf(str(pdf_path))
            return True
        except (ImportError, Exception):
            return False

    def _try_chrome_headless(self, html_path: Path, pdf_path: Path) -> bool:
        """Try using Chrome/Chromium headless for PDF generation"""
        chrome_commands = [
            'google-chrome',
            'google-chrome-stable',
            'chromium',
            'chromium-browser',
            'chrome'
        ]

        for chrome_cmd in chrome_commands:
            try:
                subprocess.run([
                    chrome_cmd,
                    '--headless',
                    '--disable-gpu',
                    '--print-to-pdf=' + str(pdf_path),
                    '--print-to-pdf-no-header',
                    f'file://{html_path.absolute()}'
                ], check=True, capture_output=True)
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue

        return False