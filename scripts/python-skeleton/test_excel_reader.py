#!/usr/bin/env python3
"""
Test script for the updated Excel reader
"""

import sys
import os
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from excel_reader import ExcelDataReader

def test_excel_reader(excel_path):
    """Test the Excel reader with a sample file"""
    try:
        print(f"Testing Excel reader with: {excel_path}")
        
        # Initialize reader
        reader = ExcelDataReader(excel_path, debug=True)
        
        # Read all tabs
        print("\n1. Reading Excel data...")
        data = reader.read_all_tabs()
        
        # Validate data
        print("\n2. Validating data...")
        reader.validate_data(data)
        
        # Process for template
        print("\n3. Processing data for template...")
        sample_id = "HG003"  # Default sample ID
        processed_data = reader.process_for_template(data, sample_id)
        
        # Print summary
        print("\n4. Processed data summary:")
        print(f"  Sample ID: {processed_data.get('sample_id', 'N/A')}")
        print(f"  Assay: {processed_data.get('assay', 'N/A')}")
        print(f"  Reference Build: {processed_data.get('reference_build', 'N/A')}")
        print(f"  Page 1 variants: {len(processed_data.get('page1_variants', []))}")
        print(f"  Page 2 variants: {len(processed_data.get('page2_variants', []))}")
        print(f"  Coverage gaps: {len(processed_data.get('coverage_gaps', []))}")
        print(f"  Overall coverage entries: {len(processed_data.get('overall_coverage', []))}")
        
        # Print QC metrics
        print("\n5. QC Metrics:")
        for key, value in processed_data.items():
            if key in ['mean_coverage', 'mapped_percent', 'duplicate_percent', 'mean_insert_size', 'acmg_covered_percent']:
                print(f"  {key}: {value}")
        
        print("\n✅ Excel reader test completed successfully!")
        return True
        
    except Exception as e:
        print(f"\n❌ Error testing Excel reader: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python test_excel_reader.py <excel_file>")
        sys.exit(1)
    
    excel_path = sys.argv[1]
    if not os.path.exists(excel_path):
        print(f"File not found: {excel_path}")
        sys.exit(1)
    
    success = test_excel_reader(excel_path)
    sys.exit(0 if success else 1)
