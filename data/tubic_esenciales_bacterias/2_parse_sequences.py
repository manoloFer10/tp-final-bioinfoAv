#!/usr/bin/env python3
"""
Step 2: Parse FASTA sequences into lookup dictionaries
Extracts DEG IDs and creates aa_map and nt_map
"""
import re
from pathlib import Path
from collections import defaultdict
import pickle

def parse_fasta_to_dict(fasta_file, verbose=True):
    """
    Parse a FASTA file and return a dict mapping deg_id -> sequence
    Handles multi-FASTA format with DEG ID extraction from headers
    """
    seq_dict = {}
    current_id = None
    current_seq = []
    
    # Regex to extract DEG IDs (e.g., DEG10530001, DEG10000123)
    deg_pattern = re.compile(r'DEG\d{4,}')
    
    duplicates = defaultdict(int)
    lines_processed = 0
    
    if verbose:
        print(f"\nParsing {fasta_file}...")
    
    with open(fasta_file, 'r') as f:
        for line in f:
            lines_processed += 1
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous sequence
                if current_id and current_seq:
                    sequence = ''.join(current_seq)
                    if current_id in seq_dict:
                        duplicates[current_id] += 1
                        # Keep first occurrence or append suffix
                        current_id = f"{current_id}_{duplicates[current_id]}"
                    seq_dict[current_id] = sequence
                
                # Extract new DEG ID from header
                match = deg_pattern.search(line)
                if match:
                    current_id = match.group(0)
                else:
                    # Fallback: use entire header (cleaned)
                    current_id = line[1:].split()[0]
                
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id and current_seq:
            sequence = ''.join(current_seq)
            if current_id in seq_dict:
                duplicates[current_id] += 1
                current_id = f"{current_id}_{duplicates[current_id]}"
            seq_dict[current_id] = sequence
    
    if verbose:
        print(f"  ✓ Processed {lines_processed:,} lines")
        print(f"  ✓ Extracted {len(seq_dict):,} sequences")
        if duplicates:
            print(f"  ⚠ Found {len(duplicates)} duplicate DEG IDs (appended suffixes)")
    
    return seq_dict

def main():
    """Parse both AA and NT sequence files"""
    script_dir = Path(__file__).parent
    
    print("=" * 60)
    print("Step 2: Parsing Sequence Files")
    print("=" * 60)
    
    # Parse amino acid sequences
    aa_file = script_dir / 'DEG10.aa'
    if aa_file.exists():
        aa_map = parse_fasta_to_dict(aa_file)
        
        # Save as pickle for fast loading
        pickle_file = script_dir / 'deg_aa_map.pkl'
        with open(pickle_file, 'wb') as f:
            pickle.dump(aa_map, f)
        print(f"  → Saved to {pickle_file.name}")
    else:
        print(f"⚠ Warning: {aa_file} not found. Run step 1 first.")
        aa_map = {}
    
    # Parse nucleotide sequences
    nt_file = script_dir / 'DEG10.nt'
    if nt_file.exists():
        nt_map = parse_fasta_to_dict(nt_file)
        
        # Save as pickle for fast loading
        pickle_file = script_dir / 'deg_nt_map.pkl'
        with open(pickle_file, 'wb') as f:
            pickle.dump(nt_map, f)
        print(f"  → Saved to {pickle_file.name}")
    else:
        print(f"⚠ Warning: {nt_file} not found. Run step 1 first.")
        nt_map = {}
    
    # Summary statistics
    if aa_map or nt_map:
        print("\n" + "=" * 60)
        print("Parsing complete!")
        print(f"  AA sequences: {len(aa_map):,}")
        print(f"  NT sequences: {len(nt_map):,}")
        
        # Sample a few sequences
        if aa_map:
            sample_ids = list(aa_map.keys())[:3]
            print(f"\n  Sample AA IDs: {', '.join(sample_ids)}")
            for sid in sample_ids[:1]:
                print(f"    {sid}: {len(aa_map[sid])} aa")
        
        print("=" * 60)

if __name__ == "__main__":
    main()
