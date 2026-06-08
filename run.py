import argparse
import sys
import os
import json
from datetime import date
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
import subprocess

from src.fetch_genome import fetch_genome
from src.run_padloc import run_padloc
from src.dbs_calculator import parse_padloc_output, calculate_dbs
from src.matching_engine import match_panel

def main():
    parser = argparse.ArgumentParser(description="PhageCompass Pipeline")
    parser.add_argument("--accession", required=True, help="Bacterial genome accession")
    parser.add_argument("--name", help="Human-readable name for output files")
    parser.add_argument("--output", default="results/", help="Output directory")
    parser.add_argument("--panel", default="results/phage_cdps.json", help="Path to phage CDPs JSON")
    parser.add_argument("--skip-fetch", action="store_true", help="Skip fetching and running PADLOC if files exist")

    args = parser.parse_args()

    name = args.name if args.name else args.accession
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    console = Console()

    fasta_path = f"data/genomes/{args.accession}.fasta"
    padloc_out_dir = "data/padloc_output"
    csv_path = os.path.join(padloc_out_dir, f"{args.accession}.fasta_padloc.csv")

    try:
        # Step 1
        console.print(f"[1/5] Fetching genome {args.accession} from NCBI...")
        if args.skip_fetch and os.path.exists(fasta_path):
            console.print("[green]✔[/green] Skipped fetch, file exists.")
        else:
            fasta_path = fetch_genome(args.accession, "data/genomes")
            if not fasta_path or not os.path.exists(fasta_path):
                raise Exception("Failed to fetch genome.")
            console.print("[green]✔[/green] Genome fetched.")

        # Step 2
        console.print("[2/5] Running PADLOC defense system detection...")
        if args.skip_fetch and os.path.exists(csv_path):
            console.print("[green]✔[/green] Skipped PADLOC, output exists.")
        else:
            csv_path = run_padloc(fasta_path, padloc_out_dir)
            if not csv_path or not os.path.exists(csv_path):
                raise Exception("Failed to run PADLOC.")
            console.print("[green]✔[/green] PADLOC complete.")

        # Step 3
        console.print("[3/5] Calculating Defense Burden Score...")
        profile = parse_padloc_output(csv_path)
        dbs_result = calculate_dbs(profile)
        console.print("[green]✔[/green] DBS calculated.")

        # Step 4
        console.print("[4/5] Loading phage panel...")
        if not os.path.exists(args.panel):
            raise Exception(f"Panel JSON not found at {args.panel}")
        with open(args.panel, "r") as f:
            phage_cdps = json.load(f)
        console.print("[green]✔[/green] Panel loaded.")

        # Step 5
        console.print("[5/5] Scoring phage-host compatibility...")
        rankings = match_panel(dbs_result, phage_cdps)
        console.print("[green]✔[/green] Scoring complete.")

        # Print table
        table = Table(title=f"PhageCompass Results: {name} (DBS: {dbs_result['DBS']})")
        table.add_column("Rank", justify="right")
        table.add_column("Phage", style="cyan")
        table.add_column("Score", justify="right")
        table.add_column("Recommendation")
        table.add_column("Countered", style="green")
        table.add_column("Uncountered", style="red")

        for i, rank in enumerate(rankings):
            score = rank["phagecompass_score"]
            if score >= 75:
                score_str = f"[green]{score}[/green]"
            elif score >= 55:
                score_str = f"[yellow]{score}[/yellow]"
            elif score >= 35:
                score_str = f"[color(208)]{score}[/color(208)]"
            else:
                score_str = f"[red]{score}[/red]"

            countered = [s["system"] for s in rank["system_breakdown"] if s["status"] == "countered"]
            uncountered = [s["system"] for s in rank["system_breakdown"] if s["status"] == "uncountered"]

            table.add_row(
                str(i+1),
                rank["phage_name"],
                score_str,
                rank["recommendation"],
                ", ".join(countered) if countered else "None",
                ", ".join(uncountered) if uncountered else "None"
            )

        console.print(table)
        print("\n")

        # Save JSON
        try:
            padloc_ver = subprocess.check_output(["padloc", "--version"], text=True).strip()
        except:
            padloc_ver = "unknown"

        top_phage = rankings[0] if rankings else {}
        top_countered = [s["system"] for s in top_phage.get("system_breakdown", []) if s["status"] == "countered"]

        report = {
            "query_accession": args.accession,
            "query_name": name,
            "run_date": date.today().isoformat(),
            "padloc_version": padloc_ver,
            "DBS": dbs_result["DBS"],
            "dominant_defense_systems": dbs_result["dominant_systems"],
            "system_count": dbs_result["system_count"],
            "phage_rankings": rankings,
            "top_recommendation": {
                "phage": top_phage.get("phage_name"),
                "score": top_phage.get("phagecompass_score"),
                "recommendation": top_phage.get("recommendation"),
                "countered_systems": top_countered
            } if top_phage else {}
        }
        
        report_path = os.path.join(output_dir, f"{name}_phagecompass_report.json")
        with open(report_path, "w") as f:
            json.dump(report, f, indent=4)

        summary = (
            f"Query:      {name} ({args.accession})\n"
            f"DBS:        {dbs_result['DBS']} / 100\n"
            f"Top phage:  {top_phage.get('phage_name')} (score: {top_phage.get('phagecompass_score')})\n"
            f"Saved to:   {report_path}"
        )
        console.print(Panel(summary, title="PhageCompass Report Summary", expand=False))

    except Exception as e:
        console.print(f"[bold red][ERROR] Step failed: {str(e)}[/bold red]")
        console.print("Pipeline aborted. Check inputs and try again.")
        sys.exit(1)

if __name__ == "__main__":
    main()
