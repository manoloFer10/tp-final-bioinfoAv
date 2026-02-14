# DEG Database Processing Pipeline for Burkholderia cenocepacia

## Overview
This pipeline processes the Database of Essential Genes (DEG) to extract and filter data specifically for *Burkholderia cenocepacia* strains J2315 and K56-2.

## Input Files (Required)
- `deg_bacteria.csv.zip` - Main DEG records with organism information
- `deg_annotation_p.csv.zip` - Dataset-level metadata and annotations
- `DEG10.aa.gz` - Protein (amino acid) sequences
- `DEG10.nt.gz` - Nucleotide sequences

## Pipeline Scripts

### `run_pipeline.py` - Master Script
Runs all pipeline steps in sequence. **Use this to execute the entire pipeline.**

```bash
python run_pipeline.py
```

### Individual Steps

#### Step 1: `1_extract_archives.py`
- Unzips CSV files from `.zip` archives
- Decompresses sequence files from `.gz` format
- Validates file formats (multi-FASTA detection)

#### Step 2: `2_parse_sequences.py`
- Parses FASTA files (DEG10.aa and DEG10.nt)
- Extracts DEG IDs using regex patterns
- Creates lookup dictionaries (deg_id → sequence)
- Saves as pickle files for fast loading

#### Step 3: `3_filter_pipeline.py`
- Loads main bacteria table
- Filters to *Burkholderia cenocepacia* records
- Assigns strain tags (J2315, K56-2, or unknown)
- Joins amino acid and nucleotide sequences
- Merges annotation metadata
- Exports multiple output formats

## Output Files

All outputs are saved to: `processed_outputs/`

### Generated Files:
1. **`deg_cenocepacia_annotations.csv`**
   - All B. cenocepacia records without sequences
   - Smaller file for quick analysis
   
2. **`deg_cenocepacia_with_sequences.csv.gz`**
   - Complete data including AA/NT sequences
   - Compressed for storage efficiency

3. **`deg_cenocepacia_J2315_annotations.csv`**
   - J2315 strain only (no sequences)

4. **`deg_cenocepacia_J2315_with_sequences.csv.gz`**
   - J2315 strain with sequences (compressed)

5. **`deg_cenocepacia_K56-2_annotations.csv`**
   - K56-2 strain only (no sequences)

6. **`deg_cenocepacia_K56-2_with_sequences.csv.gz`**
   - K56-2 strain with sequences (compressed)

### Additional Files (Intermediate):
- `deg_aa_map.pkl` - Pickled AA sequence dictionary
- `deg_nt_map.pkl` - Pickled NT sequence dictionary
- `deg_bacteria.csv` - Extracted main table
- `deg_annotation_p.csv` - Extracted annotation table
- `DEG10.aa` - Decompressed AA sequences
- `DEG10.nt` - Decompressed NT sequences

## Output Columns

### Core Columns:
- Original columns from DEG database (organism, gene names, etc.)
- `strain` - Added tag: "J2315", "K56-2", or "unknown"
- `aa_seq` - Amino acid sequence (if available)
- `nt_seq` - Nucleotide sequence (if available)
- `aa_len` - Length of amino acid sequence
- `nt_len` - Length of nucleotide sequence

## Sanity Checks

The pipeline performs automatic validation:
- ✓ Both strains present and tagged
- ✓ Record counts per strain
- ✓ Sequence coverage percentage
- ✓ Duplicate DEG ID detection
- ✓ Missing sequence statistics

## Usage Example

```bash
# Run the complete pipeline
python run_pipeline.py

# Or run individual steps
python 1_extract_archives.py
python 2_parse_sequences.py
python 3_filter_pipeline.py
```

## Requirements

```python
pandas
pickle
gzip
zipfile
pathlib
re
```

## Notes

- The pipeline automatically handles missing data
- Duplicate DEG IDs are suffixed with `_1`, `_2`, etc.
- Sequences are matched using DEG ID extraction from FASTA headers
- Low sequence coverage (<50%) triggers a warning

## Troubleshooting

### No records found for B. cenocepacia
- Check organism column names in the CSV
- Verify organism string format matches "Burkholderia cenocepacia"

### Low sequence coverage
- Verify DEG ID regex pattern matches your FASTA headers
- Check that DEG IDs in CSV match those in FASTA files

### File not found errors
- Ensure all input .zip and .gz files are present
- Run scripts in sequence (or use run_pipeline.py)
