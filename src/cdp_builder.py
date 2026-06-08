import os
import time
import json
from Bio import Entrez, SeqIO
from rich.console import Console
from rich.table import Table

Entrez.email = "phagecompass_dev@virginia.edu"

PHAGE_PANEL = {
    "phiKZ":   "NC_004629",
    "LUZ19":   "NC_010821",
    "LUZ24":   "NC_010780",
    "PAK_P1":  "NC_022746",
    "LUZ7":    "NC_013691",
    "phiPA3":  "NC_014141",
    "14-1":    "NC_012418",
    "LMA2":    "NC_019526",
    "E79":     "NC_007458", # EL phage
    "Ab31":    "NC_027399"
}

ANTI_DEFENSE_KEYWORDS = {
    "anti-RM": [
        "ocr", "arn", "ral", "dmd1", "stp",
        "restriction alleviation",
        "anti-restriction"
    ],
    "anti-CRISPR": [
        "acr", "anti-crispr", "anti-cas",
        "crispr inhibitor"
    ],
    "anti-CBASS": [
        "anti-cbass", "cap5", "nd1",
        "cbass inhibitor", "cyclic nucleotide",
        "phosphodiesterase"
    ],
    "anti-Thoeris": [
        "tad1", "anti-thoeris",
        "thoeris inhibitor"
    ],
    "anti-Retron": [
        "retron inhibitor", "anti-retron"
    ],
    "anti-Lamassu": [
        "lamassu inhibitor", "anti-lamassu"
    ]
}

def fetch_phage_genomes(panel, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for name, accession in panel.items():
        filepath = os.path.join(out_dir, f"{name}.gb")
        try:
            with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
                data = handle.read()
            with open(filepath, "w") as f:
                f.write(data)
            time.sleep(1) # Rate limit

            record = SeqIO.read(filepath, "genbank")
            length = len(record.seq)
            print(f"Fetched {name} ({accession}) — {length}bp")
            paths[name] = filepath
        except Exception as e:
            print(f"Error fetching/parsing {name}: {e}")
    return paths

def scan_genbank_for_antidefense(gb_path: str) -> dict:
    hits = {k: False for k in ANTI_DEFENSE_KEYWORDS.keys()}
    try:
        record = SeqIO.read(gb_path, "genbank")
    except Exception as e:
        print(f"Error parsing {gb_path}: {e}")
        return hits

    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0].lower()
            function_ = feature.qualifiers.get("function", [""])[0].lower()
            note = feature.qualifiers.get("note", [""])[0].lower()

            combined_text = f"{product} {function_} {note}"

            for category, keywords in ANTI_DEFENSE_KEYWORDS.items():
                for kw in keywords:
                    if kw in combined_text:
                        hits[category] = True
                        break
    return hits

def main():
    out_dir = "phagecompass/data/phage_panel"
    results_dir = "phagecompass/results"
    os.makedirs(results_dir, exist_ok=True)

    print("Fetching genomes...")
    paths = fetch_phage_genomes(PHAGE_PANEL, out_dir)

    cdps = []
    console = Console()
    table = Table(title="Phage Counter-Defense Profiles")
    table.add_column("Phage", style="cyan")
    table.add_column("Length (bp)", justify="right", style="green")
    table.add_column("CDS Count", justify="right", style="blue")
    table.add_column("Anti-Defense Found", style="magenta")
    table.add_column("Coverage", justify="right", style="yellow")
    table.add_column("Note", style="white")

    for name, acc in PHAGE_PANEL.items():
        gb_path = paths.get(name)
        if not gb_path:
            continue

        try:
            record = SeqIO.read(gb_path, "genbank")
            genome_length = len(record.seq)
            cds_count = sum(1 for f in record.features if f.type == "CDS")
        except Exception:
            genome_length = 0
            cds_count = 0

        if genome_length < 10000:
            print(f"ERROR: {name} genome too short ({genome_length}bp) — bad accession, skipping.")
            continue

        raw_hits = scan_genbank_for_antidefense(gb_path)
        anti_defense_systems = [cat for cat, detected in raw_hits.items() if detected]
        cdp_coverage = len(anti_defense_systems) / 6.0

        annotation_note = "Standard annotation scan."
        if genome_length > 200000 and cdp_coverage == 0.0:
            annotation_note = "Jumbo phage: anti-defense genes likely present but unannotated. Manual curation recommended."

        cdp = {
            "phage_name": name,
            "accession": acc,
            "genome_length": genome_length,
            "cds_count": cds_count,
            "anti_defense_systems": anti_defense_systems,
            "cdp_coverage": round(cdp_coverage, 2),
            "raw_hits": raw_hits,
            "annotation_note": annotation_note
        }
        cdps.append(cdp)

        sys_str = ", ".join(cdp["anti_defense_systems"]) if cdp["anti_defense_systems"] else "None"
        table.add_row(
            cdp["phage_name"],
            str(cdp["genome_length"]),
            str(cdp["cds_count"]),
            sys_str,
            f"{cdp['cdp_coverage']:.2f}",
            cdp["annotation_note"]
        )

    json_path = os.path.join(results_dir, "phage_cdps.json")
    with open(json_path, "w") as f:
        json.dump(cdps, f, indent=4)

    print("\n")
    console.print(table)
    print(f"\nSaved CDPs to {json_path}")

if __name__ == "__main__":
    main()
