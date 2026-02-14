#!/usr/bin/env python3
"""
Step 3: Main filtering pipeline for Burkholderia cenocepacia
Filters DEG data, adds strain tags, joins sequences, and exports outputs
"""
import pandas as pd
import pickle
import gzip
from pathlib import Path
import re

def load_sequence_maps(script_dir):
    """Load the pickled sequence dictionaries"""
    aa_map = {}
    nt_map = {}
    
    aa_pickle = script_dir / 'deg_aa_map.pkl'
    nt_pickle = script_dir / 'deg_nt_map.pkl'
    
    if aa_pickle.exists():
        with open(aa_pickle, 'rb') as f:
            aa_map = pickle.load(f)
        print(f"  ✓ Loaded {len(aa_map):,} AA sequences")
    else:
        print("  ⚠ Warning: AA sequence map not found")
    
    if nt_pickle.exists():
        with open(nt_pickle, 'rb') as f:
            nt_map = pickle.load(f)
        print(f"  ✓ Loaded {len(nt_map):,} NT sequences")
    else:
        print("  ⚠ Warning: NT sequence map not found")
    
    return aa_map, nt_map

def assign_strain(organism_str):
    """
    Assign strain tag based on organism string
    Returns: 'J2315', 'K56-2', or 'unknown'
    """
    if pd.isna(organism_str):
        return 'unknown'
    
    org_lower = str(organism_str).lower()
    
    if 'j2315' in org_lower:
        return 'J2315'
    elif 'k56-2' in org_lower or 'k562' in org_lower:
        return 'K56-2'
    else:
        return 'unknown'

def main():
    """Main pipeline execution"""
    script_dir = Path(__file__).parent
    output_dir = script_dir / 'processed_outputs'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("Step 3: Filtering Pipeline for B. cenocepacia")
    print("=" * 60)
    
    # 1. Load the main bacteria table
    print("\n[1/6] Loading deg_bacteria.csv...")
    bacteria_file = script_dir / 'deg_bacteria.csv'
    if not bacteria_file.exists():
        print(f"❌ Error: {bacteria_file} not found. Run step 1 first.")
        return
    
    # The DEG bacteria file uses semicolons and has NO header row
    # Expected columns based on typical DEG format:
    column_names = [
        'organism', 'reference', 'pmid', 'accession', 'essential_genes', 
        'total_genes', 'secondary_reference', 'secondary_pmid', 
        'tertiary_reference', 'quaternary_reference', 'method', 
        'condition', 'dataset', 'date'
    ]
    
    try:
        df = pd.read_csv(
            bacteria_file, 
            sep=';', 
            header=None, 
            names=column_names,
            low_memory=False, 
            on_bad_lines='skip', 
            encoding='utf-8',
            dtype=str
        )
    except Exception as e:
        print(f"  ⚠ UTF-8 encoding failed, trying latin-1...")
        try:
            df = pd.read_csv(
                bacteria_file, 
                sep=';', 
                header=None,
                names=column_names,
                low_memory=False, 
                on_bad_lines='skip', 
                encoding='latin-1',
                dtype=str
            )
        except Exception as e2:
            print(f"  ❌ Error loading CSV: {e2}")
            return
    
    # Clean quotes from values
    for col in df.columns:
        if df[col].dtype == 'object':
            df[col] = df[col].str.strip('"')
    
    print(f"  ✓ Loaded {len(df):,} total records")
    print(f"  Columns: {list(df.columns)}")
    
    # 2. Filter to Burkholderia cenocepacia
    print("\n[2/6] Filtering to Burkholderia cenocepacia...")
    
    organism_col = 'organism'  # We set this in column_names
    
    if organism_col in df.columns:
        print(f"  Using column: '{organism_col}'")
        mask = df[organism_col].str.contains(
            'Burkholderia cenocepacia', 
            case=False, 
            na=False
        )
        df_filtered = df[mask].copy()
        print(f"  ✓ Found {len(df_filtered):,} Burkholderia cenocepacia records")
    else:
        print("  ❌ Error: Cannot find organism column")
        print(f"  Available columns: {list(df.columns)}")
        return
    
    if len(df_filtered) == 0:
        print("  ⚠ Warning: No records found for Burkholderia cenocepacia")
        print("  Showing sample organism values:")
        print(df[organism_col].value_counts().head(10))
        return
    
    # 3. Add strain tag
    print("\n[3/6] Adding strain tags...")
    df_filtered['strain'] = df_filtered[organism_col].apply(assign_strain)
    
    strain_counts = df_filtered['strain'].value_counts()
    print("  Strain distribution:")
    for strain, count in strain_counts.items():
        print(f"    {strain}: {count:,} records")
    
    # Optionally drop 'unknown' strains
    if 'unknown' in strain_counts.index and strain_counts['unknown'] > 0:
        print(f"  ⚠ Keeping {strain_counts['unknown']} unknown strain records")
    
    # 4. Load sequence maps
    print("\n[4/6] Loading sequence data...")
    aa_map, nt_map = load_sequence_maps(script_dir)
    
    # 5. Load annotation file to get individual gene records with DEG IDs
    print("\n[5/6] Loading annotation file with gene-level data...")
    annotation_file = script_dir / 'deg_annotation_p.csv'
    
    if not annotation_file.exists():
        print(f"  ❌ Error: {annotation_file} not found. This file contains gene-level records.")
        return
    
    # Define proper column names for the annotation file (no header row)
    annot_column_names = [
        'dataset', 'deg_id', 'gene_name', 'gene_id', 'cog', 
        'functional_category', 'function', 'organism', 'accession',
        'growth_condition', 'expression_level', 'gene_ontology', 
        'uniprot_id', 'extra'
    ]
    
    try:
        df_annot = pd.read_csv(
            annotation_file, 
            sep=';', 
            header=None,
            names=annot_column_names,
            low_memory=False, 
            on_bad_lines='skip', 
            dtype=str
        )
        # Clean quotes
        for col in df_annot.columns:
            if df_annot[col].dtype == 'object':
                df_annot[col] = df_annot[col].str.strip('"')
    except Exception as e:
        print(f"  ❌ Error loading annotation file: {e}")
        return
    
    print(f"  ✓ Loaded {len(df_annot):,} gene records")
    print(f"  Columns: {list(df_annot.columns)}")
    
    # Use the explicitly defined column names
    dataset_col = 'dataset'
    deg_id_col = 'deg_id'
    
    # Get datasets from filtered bacteria records
    filtered_datasets = df_filtered['dataset'].unique()
    print(f"  ✓ Filtering to {len(filtered_datasets)} B. cenocepacia datasets")
    
    # Filter annotation to matching datasets
    df_genes = df_annot[df_annot[dataset_col].isin(filtered_datasets)].copy()
    print(f"  ✓ Found {len(df_genes):,} gene records for B. cenocepacia")
    
    if len(df_genes) == 0:
        print("  ⚠ Warning: No gene records found for B. cenocepacia datasets")
        return
    
    # Merge strain information from  bacteria file
    strain_map = df_filtered.set_index('dataset')['strain'].to_dict()
    df_genes['strain'] = df_genes[dataset_col].map(strain_map)
    
    # Also merge other useful info from bacteria file
    for col in ['organism', 'method', 'condition', 'pmid', 'reference']:
        if col in df_filtered.columns:
            col_map = df_filtered.set_index('dataset')[col].to_dict()
            if col in df_genes.columns:
                df_genes[f'{col}_dataset'] = df_genes[dataset_col].map(col_map)
            else:
                df_genes[col] = df_genes[dataset_col].map(col_map)
    
    # 6. Join sequences
    print("\n[6/6] Joining sequences to gene records...")
    
    # We already have deg_id_col defined above
    print(f"  Using DEG ID column: '{deg_id_col}'")
    # We already have deg_id_col defined above
    print(f"  Using DEG ID column: '{deg_id_col}'")
    
    # Map sequences to df_genes
    df_genes['aa_seq'] = df_genes[deg_id_col].map(aa_map)
    df_genes['nt_seq'] = df_genes[deg_id_col].map(nt_map)
    
    # Add length columns
    df_genes['aa_len'] = df_genes['aa_seq'].apply(
        lambda x: len(x) if pd.notna(x) else None
    )
    df_genes['nt_len'] = df_genes['nt_seq'].apply(
        lambda x: len(x) if pd.notna(x) else None
    )
    
    # Statistics
    aa_coverage = df_genes['aa_seq'].notna().sum() / len(df_genes) * 100
    nt_coverage = df_genes['nt_seq'].notna().sum() / len(df_genes) * 100
    
    print(f"  ✓ AA sequence coverage: {aa_coverage:.1f}%")
    print(f"  ✓ NT sequence coverage: {nt_coverage:.1f}%")
    
    if aa_coverage < 50:
        print("  ⚠ Warning: Low AA sequence coverage. Check DEG ID extraction.")
    
    # Use df_genes (gene-level data) as our final filtered dataframe
    df_filtered = df_genes
    
    # 7. Export outputs
    print("\n" + "=" * 60)
    print("Exporting outputs...")
    print("=" * 60)
    
    # Output 1: Annotations only (no sequences)
    cols_no_seq = [col for col in df_filtered.columns 
                   if col not in ['aa_seq', 'nt_seq']]
    df_annotations = df_filtered[cols_no_seq]
    
    out_annot = output_dir / 'deg_cenocepacia_annotations.csv'
    df_annotations.to_csv(out_annot, index=False)
    print(f"\n[1/4] ✓ Saved: {out_annot.name}")
    print(f"       {len(df_annotations):,} records, {len(df_annotations.columns)} columns")
    
    # Output 2: Full data with sequences (compressed)
    out_full = output_dir / 'deg_cenocepacia_with_sequences.csv'
    df_filtered.to_csv(out_full, index=False)
    
    # Compress it
    with open(out_full, 'rb') as f_in:
        with gzip.open(str(out_full) + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
    out_full.unlink()  # Remove uncompressed version
    
    print(f"[2/4] ✓ Saved: {out_full.name}.gz")
    print(f"       {len(df_filtered):,} records, {len(df_filtered.columns)} columns")
    
    # Output 3 & 4: Per-strain splits
    for strain in ['J2315', 'K56-2']:
        df_strain = df_filtered[df_filtered['strain'] == strain]
        
        if len(df_strain) > 0:
            # Annotations
            out_strain_annot = output_dir / f'deg_cenocepacia_{strain}_annotations.csv'
            df_strain[cols_no_seq].to_csv(out_strain_annot, index=False)
            
            # With sequences
            out_strain_full = output_dir / f'deg_cenocepacia_{strain}_with_sequences.csv'
            df_strain.to_csv(out_strain_full, index=False)
            
            with open(out_strain_full, 'rb') as f_in:
                with gzip.open(str(out_strain_full) + '.gz', 'wb') as f_out:
                    f_out.writelines(f_in)
            out_strain_full.unlink()
            
            idx = 3 if strain == 'J2315' else 4
            print(f"[{idx}/4] ✓ Saved: {strain} strain files")
            print(f"       {len(df_strain):,} records")
    
    # Final sanity checks
    print("\n" + "=" * 60)
    print("Sanity Checks")
    print("=" * 60)
    
    print(f"\n✓ Total records: {len(df_filtered):,}")
    print(f"✓ Strains found: {df_filtered['strain'].nunique()}")
    print(f"✓ Strain breakdown:")
    for strain, count in df_filtered['strain'].value_counts().items():
        pct = count / len(df_filtered) * 100
        print(f"    {strain}: {count:,} ({pct:.1f}%)")
    
    duplicates = df_filtered['deg_id'].duplicated().sum()
    print(f"\n✓ Duplicate DEG IDs: {duplicates}")
    
    if 'aa_seq' in df_filtered.columns:
        missing_aa = df_filtered['aa_seq'].isna().sum()
        missing_aa_pct = missing_aa / len(df_filtered) * 100
        print(f"✓ Missing AA sequences: {missing_aa} ({missing_aa_pct:.1f}%)")
    
    if 'nt_seq' in df_filtered.columns:
        missing_nt = df_filtered['nt_seq'].isna().sum()
        missing_nt_pct = missing_nt / len(df_filtered) * 100
        print(f"✓ Missing NT sequences: {missing_nt} ({missing_nt_pct:.1f}%)")
    
    print("\n" + "=" * 60)
    print("Pipeline complete! ✓")
    print("=" * 60)
    print(f"\nOutputs saved to: {output_dir}")

if __name__ == "__main__":
    main()
