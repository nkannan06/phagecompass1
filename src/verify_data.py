from Bio import Entrez, SeqIO
import json
import time
import os
from rich.console import Console
from rich.table import Table
import subprocess
from datetime import date

Entrez.email = "nitishkannan2024@gmail.com"

BACTERIAL_GENOMES = {
    "PAO1":          "AE004091",
    "PA7":           "NC_011770",
    "K_pneumoniae":  "CP012001"
}

PHAGE_PANEL = {
    "phiKZ":   "NC_004629",
    "LUZ19":   "NC_010821",
    "LUZ24":   "NC_010779",
    "PAK_P1":  "NC_022746",
    "LUZ7":    "NC_013691",
    "phiPA3":  "JN699037",
    "14-1":    "NC_014169",
    "LMA2":    "NC_019526",
    "E79":     "NC_007042",
    "Ab31":    "NC_027399"
}

def verify_accession(accession: str) -> dict:
    try:
        handle = Entrez.esummary(
            db="nucleotide",
            id=accession,
            retmode="json"
        )
        summary = json.loads(handle.read())
        uid = list(summary["result"].keys())[1]
        record = summary["result"][uid]
        return {
            "accession":   accession,
            "title":       record["title"],
            "organism":    record.get("organism", ""),
            "length":      record["slen"],
            "update_date": record["updatedate"]
        }
    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return {}

def check_local_gb(filepath: str) -> dict:
    try:
        record = SeqIO.read(filepath, "genbank")
        return {
            "local_id":          record.id,
            "local_length":      len(record.seq),
            "local_organism":    record.annotations.get("organism", "unknown"),
            "local_cds_count":   sum(1 for f in record.features if f.type == "CDS")
        }
    except Exception as e:
        return {"error": str(e)}

def run_padloc_version():
    try:
        version_out = subprocess.check_output(["padloc", "--version"], text=True).strip()
        db_list_out = subprocess.check_output(["padloc", "--db-list"], text=True).strip()
        return version_out, db_list_out
    except Exception:
        return "unknown", "unknown"

def main():
    console = Console()
    report = {
        "padloc_version": "",
        "padloc_db_version": "",
        "run_date": date.today().isoformat(),
        "bacterial_genomes": [],
        "phage_panel": [],
        "local_file_checks": []
    }

    # Step 2: Bacteria
    bact_table = Table(title="Bacterial Genome Verification")
    bact_table.add_column("Accession")
    bact_table.add_column("Title")
    bact_table.add_column("Organism")
    bact_table.add_column("Length")

    for name, acc in BACTERIAL_GENOMES.items():
        meta = verify_accession(acc)
        if meta:
            report["bacterial_genomes"].append(meta)
            bact_table.add_row(meta["accession"], meta["title"][:30]+"...", meta["organism"], str(meta["length"]))
        time.sleep(1)

    console.print(bact_table)

    # Step 3: Phage
    phage_table = Table(title="Phage Panel Verification")
    phage_table.add_column("Name")
    phage_table.add_column("Accession")
    phage_table.add_column("Title")
    phage_table.add_column("Organism")
    phage_table.add_column("Length")

    red_flags = []

    for name, acc in PHAGE_PANEL.items():
        meta = verify_accession(acc)
        if meta:
            report["phage_panel"].append({"name": name, **meta})
            phage_table.add_row(name, meta["accession"], meta["title"][:30]+"...", meta["organism"], str(meta["length"]))
            
            # Check for red flags
            org_lower = meta["organism"].lower()
            title_lower = meta["title"].lower()
            if "phage" not in org_lower and "virus" not in org_lower:
                red_flags.append(f"{name} ({acc}): Organism '{meta['organism']}' might not be a phage/virus.")
            if "chromosome" in title_lower or "plasmid" in title_lower:
                red_flags.append(f"{name} ({acc}): Title suggests chromosome/plasmid.")
            if meta["length"] < 10000:
                red_flags.append(f"{name} ({acc}): Length {meta['length']}bp is < 10000bp.")
            if meta["length"] > 1000000:
                red_flags.append(f"{name} ({acc}): Length {meta['length']}bp is > 1000000bp (likely a bacterium).")

        time.sleep(1)

    console.print(phage_table)
    
    if red_flags:
        console.print("[red bold]\nRED FLAGS:[/red bold]")
        for flag in red_flags:
            console.print(f"[red]- {flag}[/red]")
    else:
        console.print("\n[green]No Red Flags Triggered![/green]")

    # Step 4: Local File checks
    local_dir = "phagecompass/data/phage_panel"
    local_comp_table = Table(title="Local vs NCBI Comparison")
    local_comp_table.add_column("Phage")
    local_comp_table.add_column("NCBI Length")
    local_comp_table.add_column("Local Length")
    local_comp_table.add_column("Mismatch")

    for p in report["phage_panel"]:
        name = p["name"]
        ncbi_len = p["length"]
        gb_path = os.path.join(local_dir, f"{name}.gb")
        
        mismatch = "File not found"
        if os.path.exists(gb_path):
            local_meta = check_local_gb(gb_path)
            if "error" not in local_meta:
                loc_len = local_meta["local_length"]
                diff = abs(ncbi_len - loc_len)
                mismatch = f"Diff: {diff}" if diff > 100 else "OK"
                report["local_file_checks"].append({
                    "name": name,
                    "ncbi_length": ncbi_len,
                    "local_length": loc_len,
                    "mismatch": mismatch
                })
                local_comp_table.add_row(name, str(ncbi_len), str(loc_len), mismatch)
            else:
                local_comp_table.add_row(name, str(ncbi_len), "ERROR", local_meta["error"])
        else:
            local_comp_table.add_row(name, str(ncbi_len), "N/A", "File not found")

    console.print(local_comp_table)

    # Step 5: PADLOC
    padloc_ver, padloc_db_list = run_padloc_version()
    print(f"\nPADLOC Version:\n{padloc_ver}")
    print(f"\nPADLOC DB List:\n{padloc_db_list}")
    
    report["padloc_version"] = padloc_ver
    report["padloc_db_version"] = padloc_db_list

    # Step 6: Save provenance report
    os.makedirs("phagecompass/results", exist_ok=True)
    with open("phagecompass/results/provenance_report.json", "w") as f:
        json.dump(report, f, indent=4)
        
    print("\nSaved provenance report to phagecompass/results/provenance_report.json")

if __name__ == "__main__":
    main()
