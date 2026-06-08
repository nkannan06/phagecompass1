import os
from Bio import Entrez
import urllib.error

# Note: Entrez requires an email. Replacing with a placeholder UVA email format.
Entrez.email = "phagecompass_dev@virginia.edu"

def fetch_genome(accession: str, out_dir: str) -> str:
    """
    Fetches a genome by accession from NCBI and saves it as a FASTA file.
    """
    os.makedirs(out_dir, exist_ok=True)
    filepath = os.path.join(out_dir, f"{accession}.fasta")
    
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            fasta_data = handle.read()
            
        with open(filepath, "w") as f:
            f.write(fasta_data)
            
        print(f"Downloaded {accession} -> {filepath}")
        return filepath
        
    except urllib.error.URLError as e:
        print(f"Network error while downloading {accession}: {e}")
        return ""
    except Exception as e:
        print(f"Error downloading {accession}: {e}")
        return ""
