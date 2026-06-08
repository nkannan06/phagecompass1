import os
import subprocess
import glob
import pandas as pd

def run_padloc(fasta_path: str, out_dir: str) -> str:
    """
    Runs PADLOC on a given FASTA file and returns the path to the output CSV.
    """
    os.makedirs(out_dir, exist_ok=True)
    accession = os.path.splitext(os.path.basename(fasta_path))[0]
    
    # The padloc command
    cmd = ["padloc", "--fna", fasta_path, "--outdir", out_dir]
    
    try:
        # Run PADLOC programmatically
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running PADLOC on {fasta_path}:\n{e.stderr}")
        return ""
    
    # PADLOC generates output CSVs in the form: <fasta_prefix>_padloc.csv
    csv_pattern = os.path.join(out_dir, f"{accession}*.csv")
    csv_files = glob.glob(csv_pattern)
    
    if not csv_files:
        print(f"PADLOC complete for {accession} — 0 defense systems detected (No CSV generated)")
        return ""
        
    csv_path = csv_files[0]
    
    # Count defense systems by reading the CSV
    try:
        df = pd.read_csv(csv_path)
        n_systems = len(df)
    except pd.errors.EmptyDataError:
        n_systems = 0
        
    print(f"PADLOC complete for {accession} — {n_systems} defense systems detected")
    
    return csv_path
