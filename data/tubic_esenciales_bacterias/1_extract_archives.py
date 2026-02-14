#!/usr/bin/env python3
"""
Step 1: Extract compressed archives
Unzip CSV files and decompress gzipped sequence files
"""
import os
import zipfile
import gzip
import shutil
from pathlib import Path

def extract_archives():
    """Extract all compressed files in the current directory"""
    script_dir = Path(__file__).parent
    os.chdir(script_dir)
    
    print("=" * 60)
    print("Step 1: Extracting Archives")
    print("=" * 60)
    
    # Unzip CSV files
    csv_zips = ['deg_bacteria.csv.zip', 'deg_annotation_p.csv.zip']
    for zip_name in csv_zips:
        if os.path.exists(zip_name):
            print(f"\n[1/2] Extracting {zip_name}...")
            with zipfile.ZipFile(zip_name, 'r') as zip_ref:
                zip_ref.extractall('.')
            print(f"✓ Extracted {zip_name}")
        else:
            print(f"⚠ Warning: {zip_name} not found")
    
    # Decompress gzipped sequence files
    gz_files = ['DEG10.aa.gz', 'DEG10.nt.gz']
    for gz_name in gz_files:
        out_name = gz_name.replace('.gz', '')
        if os.path.exists(gz_name):
            print(f"\n[2/2] Decompressing {gz_name}...")
            with gzip.open(gz_name, 'rb') as f_in:
                with open(out_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"✓ Decompressed to {out_name}")
            
            # Check file type
            with open(out_name, 'rb') as f:
                header = f.read(100).decode('utf-8', errors='ignore')
                if header.startswith('>'):
                    print(f"  → Detected multi-FASTA format")
                else:
                    print(f"  ⚠ Warning: Unexpected format (not FASTA)")
        else:
            print(f"⚠ Warning: {gz_name} not found")
    
    print("\n" + "=" * 60)
    print("Extraction complete!")
    print("=" * 60)

if __name__ == "__main__":
    extract_archives()
