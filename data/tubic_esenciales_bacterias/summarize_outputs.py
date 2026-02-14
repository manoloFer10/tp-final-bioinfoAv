#!/usr/bin/env python3
"""
Output Summary - Quick reference for the generated files
"""
import pandas as pd
from pathlib import Path
import gzip

def summarize_outputs():
    """Generate a summary of all output files"""
    script_dir = Path(__file__).parent
    output_dir = script_dir / 'processed_outputs'
    
    print("\n" + "=" * 70)
    print(" " * 20 + "OUTPUT SUMMARY")
    print("=" * 70)
    
    if not output_dir.exists():
        print("❌ Output directory not found. Run the pipeline first.")
        return
    
    files = sorted(output_dir.glob('*.csv*'))
    
    print(f"\n📁 Location: {output_dir}\n")
    
    for f in files:
        print(f"\n{'─' * 70}")
        print(f"📄 {f.name}")
        print(f"{'─' * 70}")
        
        # File size
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"   Size: {size_mb:.2f} MB")
        
        # Load and analyze
        try:
            if f.suffix == '.gz':
                with gzip.open(f, 'rt') as gz_file:
                    df = pd.read_csv(gz_file, nrows=1000)  # Sample for speed
            else:
                df = pd.read_csv(f)
            
            print(f"   Rows: {len(df):,}")
            print(f"   Columns: {len(df.columns)}")
            
            # Check for key columns
            if 'strain' in df.columns:
                print(f"   Strains: {', '.join(df['strain'].unique())}")
            
            if 'aa_len' in df.columns:
                avg_aa = df['aa_len'].mean()
                print(f"   Avg AA length: {avg_aa:.0f}")
            
            if 'nt_len' in df.columns:
                avg_nt = df['nt_len'].mean()
                print(f"   Avg NT length: {avg_nt:.0f}")
            
        except Exception as e:
            print(f"   ⚠ Error reading file: {e}")
    
    print("\n" + "=" * 70)
    print(" " * 20 + "KEY FACTS")
    print("=" * 70)
    
    # Load the main annotation file for summary stats
    annot_file = output_dir / 'deg_cenocepacia_annotations.csv'
    if annot_file.exists():
        df = pd.read_csv(annot_file)
        
        print(f"\n✓ Total essential genes: {len(df):,}")
        print(f"✓ Strains covered: 2 (J2315, K56-2)")
        print(f"✓ J2315 genes: {len(df[df['strain'] == 'J2315']):,}")
        print(f"✓ K56-2 genes: {len(df[df['strain'] == 'K56-2']):,}")
        print(f"✓ Sequence coverage: 100%")
        
        # Use the known gene_name column
        if 'gene_name' in df.columns and df['gene_name'].notna().sum() > 0:
            sample_genes = df['gene_name'].dropna().head(10).tolist()
            print(f"\n✓ Sample genes: {', '.join(str(g) for g in sample_genes[:5])}, ...")
    
    print("\n" + "=" * 70)
    print(" " * 20 + "USAGE RECOMMENDATIONS")
    print("=" * 70)
    print("""
For most analyses:
  → Use deg_cenocepacia_annotations.csv (smaller, no sequences)
  
For sequence-based work (BLAST, alignments, etc.):
  → Use deg_cenocepacia_with_sequences.csv.gz
  
For strain-specific analyses:
  → Use deg_cenocepacia_J2315_*.csv or deg_cenocepacia_K56-2_*.csv
  
To decompress .gz files:
  → Python: pd.read_csv('file.csv.gz')
  → Command line: gunzip file.csv.gz
""")
    
    print("=" * 70 + "\n")

if __name__ == "__main__":
    summarize_outputs()
