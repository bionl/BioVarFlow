#!/usr/bin/env python3
"""
Test script for the sample summary functionality
"""

import os
import tempfile
import pandas as pd

def create_test_files():
    """Create test files to verify the sample summary parsing"""
    
    # Create test flagstat file
    flagstat_content = """1000000 + 0 in total (QC-passed reads + QC-failed reads)
950000 + 0 mapped (95.00% : N/A)
50000 + 0 duplicates (5.00% : N/A)
"""
    
    # Create test stats file
    stats_content = """# This file was produced by samtools stats (1.17) using the command: samtools stats
SN	total length:	1000000000
SN	bases mapped:	950000000
SN	Q30 bases:	900000000
"""
    
    # Create test sex check file
    sex_check_content = """Sample: HG003
X chromosome depth: 15000
Y chromosome depth: 5000
X/Y ratio: 3.000
Predicted sex: Male
"""
    
    # Create test bcftools stats file
    bcftools_stats_content = """# This file was produced by bcftools stats (1.17) using the command: bcftools stats
TSTV	1000	2000	2.000
HET	500	1000	2.000
"""
    
    # Create test coverage summary file
    coverage_content = """Sample: HG003
Regions: 1000
Mean_coverage: 45.2
Pct_>=10x: 98.5%
Pct_>=20x: 95.2%
Pct_>=30x: 90.1%
"""
    
    # Create test files
    test_files = {
        'flagstat.txt': flagstat_content,
        'stats.txt': stats_content,
        'sex_check.txt': sex_check_content,
        'bcftools_stats.txt': bcftools_stats_content,
        'coverage_summary.txt': coverage_content
    }
    
    temp_dir = tempfile.mkdtemp()
    file_paths = {}
    
    for filename, content in test_files.items():
        filepath = os.path.join(temp_dir, filename)
        with open(filepath, 'w') as f:
            f.write(content)
        file_paths[filename] = filepath
    
    return temp_dir, file_paths

def test_sample_summary_parsing():
    """Test the sample summary parsing functionality"""
    
    # Import the function from the main script
    import sys
    sys.path.append('.')
    
    # Create test files
    temp_dir, file_paths = create_test_files()
    
    try:
        # Import the parsing function
        from generate_lean_report_org import parse_sample_summary_data
        
        # Temporarily set the global variables
        import generate_lean_report_org
        generate_lean_report_org.flagstat_file = file_paths['flagstat.txt']
        generate_lean_report_org.stats_file = file_paths['stats.txt']
        generate_lean_report_org.sex_check_file = file_paths['sex_check.txt']
        generate_lean_report_org.bcftools_stats_file = file_paths['bcftools_stats.txt']
        generate_lean_report_org.coverage_summary = file_paths['coverage_summary.txt']
        
        # Parse the data
        summary_data = parse_sample_summary_data()
        
        # Print results
        print("Sample Summary Test Results:")
        print("=" * 40)
        for key, value in summary_data.items():
            print(f"{key}: {value}")
        
        # Verify expected values
        expected_values = {
            'Total_Reads': 1000000,
            'Mapped_Percent': 95.0,
            'Duplicate_Percent': 5.0,
            'Yield_Gb': 1.0,
            'Q30_Percent': 90.0,
            'X_Depth': 15000,
            'Y_Depth': 5000,
            'X_Y_Ratio': 3.0,
            'Predicted_Sex': 'Male',
            'Ti_Tv_Ratio': 2.0,
            'Het_Hom_Ratio': 2.0,
            'Mean_Coverage': 45.2,
            'Pct_10x': 98.5,
            'Pct_20x': 95.2,
            'Pct_30x': 90.1
        }
        
        print("\nVerification:")
        print("=" * 40)
        all_passed = True
        for key, expected in expected_values.items():
            actual = summary_data.get(key)
            if actual == expected:
                print(f"✓ {key}: {actual}")
            else:
                print(f"✗ {key}: expected {expected}, got {actual}")
                all_passed = False
        
        if all_passed:
            print("\n✓ All tests passed!")
        else:
            print("\n✗ Some tests failed!")
            
    finally:
        # Clean up
        import shutil
        shutil.rmtree(temp_dir)

if __name__ == "__main__":
    test_sample_summary_parsing() 