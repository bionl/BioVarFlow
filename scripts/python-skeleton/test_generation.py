#!/usr/bin/env python3
"""
Test script to validate the report generation process.
Run this script to test the complete pipeline with sample data.
"""

import os
import sys
from pathlib import Path
from generate_report import main

def test_report_generation():
    """Test the complete report generation pipeline."""

    # Get paths
    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    sample_data_dir = project_root / "sample-data" / "HG003_variants_lean_v1"

    # Check if sample data exists
    if not sample_data_dir.exists():
        print(f"❌ Sample data directory not found: {sample_data_dir}")
        return False

    print("🧪 Testing Report Generation Pipeline")
    print("=" * 50)

    # List available CSV files
    csv_files = list(sample_data_dir.glob("*.csv"))
    print(f"📁 Found {len(csv_files)} CSV files in sample data:")
    for csv_file in csv_files:
        print(f"   - {csv_file.name}")

    # Test the main generation function
    try:
        print("\n🚀 Starting report generation...")

        # Create output directory
        output_dir = current_dir / "output"
        output_dir.mkdir(exist_ok=True)

        # Run the main generation
        excel_path = str(sample_data_dir)  # Directory containing CSV files
        output_path = str(output_dir / "test_report.html")

        success = main(excel_path, output_path)

        if success:
            print(f"✅ Report generated successfully!")
            print(f"📄 Output saved to: {output_path}")

            # Check if file was created
            if Path(output_path).exists():
                file_size = Path(output_path).stat().st_size
                print(f"📊 File size: {file_size:,} bytes")
                return True
            else:
                print("❌ Output file was not created")
                return False
        else:
            print("❌ Report generation failed")
            return False

    except Exception as e:
        print(f"❌ Error during generation: {str(e)}")
        print(f"🔍 Error type: {type(e).__name__}")
        import traceback
        traceback.print_exc()
        return False

def validate_template_files():
    """Validate that all required template files exist."""

    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    template_dir = project_root / "template-files"

    required_files = [
        "index.html",
        "style.css"
    ]

    print("\n🔍 Validating Template Files")
    print("=" * 30)

    all_valid = True
    for filename in required_files:
        file_path = template_dir / filename
        if file_path.exists():
            size = file_path.stat().st_size
            print(f"✅ {filename} ({size:,} bytes)")
        else:
            print(f"❌ {filename} - MISSING")
            all_valid = False

    return all_valid

def validate_sample_data():
    """Validate that all required sample data files exist."""

    current_dir = Path(__file__).parent
    project_root = current_dir.parent
    sample_dir = project_root / "sample-data" / "HG003_variants_lean_v1"

    required_files = [
        "Sample Summary-Table 1.csv",
        "ACMG SF (Reportable)-Table 1.csv",
        "ACMG Genes Coverage-Table 1.csv",
        "Coverage gaps-Table 1.csv",
        "PASS variants-Table 1.csv"
    ]

    print("\n🔍 Validating Sample Data Files")
    print("=" * 35)

    all_valid = True
    for filename in required_files:
        file_path = sample_dir / filename
        if file_path.exists():
            size = file_path.stat().st_size
            print(f"✅ {filename} ({size:,} bytes)")
        else:
            print(f"❌ {filename} - MISSING")
            all_valid = False

    return all_valid

if __name__ == "__main__":
    print("🎯 ACMG SF Variants Report - Test Suite")
    print("=" * 50)

    # Validate prerequisites
    template_valid = validate_template_files()
    sample_valid = validate_sample_data()

    if not template_valid or not sample_valid:
        print("\n❌ Prerequisites not met. Please ensure all required files exist.")
        sys.exit(1)

    # Run the test
    success = test_report_generation()

    if success:
        print("\n🎉 All tests passed! The report generation pipeline is working correctly.")
        sys.exit(0)
    else:
        print("\n❌ Tests failed. Please check the error messages above.")
        sys.exit(1)