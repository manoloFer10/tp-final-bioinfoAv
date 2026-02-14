#!/usr/bin/env python3
"""
Master Pipeline Runner
Executes all steps in sequence to process DEG data for Burkholderia cenocepacia
"""
import sys
import subprocess
from pathlib import Path
import time

def run_script(script_name, description):
    """Run a Python script and handle errors"""
    print("\n" + "█" * 70)
    print(f"█  {description}")
    print("█" * 70 + "\n")
    
    script_path = Path(__file__).parent / script_name
    
    if not script_path.exists():
        print(f"❌ Error: Script not found: {script_name}")
        return False
    
    start_time = time.time()
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=False,
            text=True,
            check=True
        )
        
        elapsed = time.time() - start_time
        print(f"\n✓ Completed in {elapsed:.1f}s")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Error running {script_name}")
        print(f"Exit code: {e.returncode}")
        return False
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        return False

def main():
    """Run the complete pipeline"""
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  DEG Database Processing Pipeline for B. cenocepacia  ".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "═" * 68 + "╝")
    
    pipeline_steps = [
        ("1_extract_archives.py", "Step 1: Extract compressed archives"),
        ("2_parse_sequences.py", "Step 2: Parse FASTA sequences"),
        ("3_filter_pipeline.py", "Step 3: Filter and export B. cenocepacia data"),
    ]
    
    overall_start = time.time()
    
    for script, description in pipeline_steps:
        success = run_script(script, description)
        
        if not success:
            print("\n" + "=" * 70)
            print("❌ Pipeline stopped due to error")
            print("=" * 70)
            sys.exit(1)
    
    # Success summary
    overall_elapsed = time.time() - overall_start
    
    print("\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  ✓ PIPELINE COMPLETED SUCCESSFULLY  ".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("║" + f"  Total time: {overall_elapsed:.1f}s  ".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "═" * 68 + "╝")
    
    # Show output location
    output_dir = Path(__file__).parent / 'processed_outputs'
    print(f"\n📁 Output files saved to:")
    print(f"   {output_dir}")
    
    if output_dir.exists():
        print(f"\n📄 Generated files:")
        for f in sorted(output_dir.iterdir()):
            size_mb = f.stat().st_size / (1024 * 1024)
            print(f"   • {f.name} ({size_mb:.2f} MB)")
    
    print("\n")

if __name__ == "__main__":
    main()
