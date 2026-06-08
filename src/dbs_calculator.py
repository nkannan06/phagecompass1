import pandas as pd
import math
from rich.console import Console
from rich.table import Table

SYSTEM_WEIGHTS = {
    "CBASS":            8.0,
    "CRISPR":           9.0,
    "RM_type_I":        7.0,
    "RM_type_II":       6.0,
    "RM_type_III":      6.5,
    "Thoeris":          7.0,
    "Zorya":            6.5,
    "Druantia":         6.0,
    "Gabija":           5.5,
    "Pycsar":           7.0,
    "retron":           5.0,
    "Lamassu":          5.5,
    "Avs":              6.0,
    "SoFic":            4.5,
    "Helicase":         3.5,
    "PDC":              4.0,
    "PD-":              4.0,
    "DEFAULT":          3.5
}

def get_system_weight(system_name: str) -> float:
    """Maps system names to weight categories using prefix matching."""
    system_name_lower = system_name.lower()
    
    # Explicit check for cas systems to map to CRISPR
    if system_name_lower.startswith("cas"):
        return SYSTEM_WEIGHTS["CRISPR"]
        
    for key, weight in SYSTEM_WEIGHTS.items():
        if key != "DEFAULT" and system_name_lower.startswith(key.lower()):
            return weight
    return SYSTEM_WEIGHTS["DEFAULT"]

def parse_padloc_output(csv_path: str) -> list[dict]:
    """Parses a PADLOC CSV output into defense system records."""
    try:
        df = pd.read_csv(csv_path)
    except pd.errors.EmptyDataError:
        return []
    except FileNotFoundError:
        print(f"File not found: {csv_path}")
        return []

    if df.empty:
        return []

    profile = []
    # Group rows by system.number
    grouped = df.groupby('system.number')
    for sys_num, group in grouped:
        system_name = group['system'].iloc[0]

        # Completeness: mean of target.coverage and hmm.coverage for all rows
        target_cov = group['target.coverage'].tolist()
        hmm_cov = group['hmm.coverage'].tolist()
        all_covs = target_cov + hmm_cov
        completeness = sum(all_covs) / len(all_covs) if all_covs else 0.0

        component_count = len(group)

        profile.append({
            "system_name": system_name,
            "system_number": int(sys_num),
            "completeness": completeness,
            "component_count": component_count
        })

    return profile

def calculate_dbs(defense_profile: list[dict]) -> dict:
    """Calculates DBS and annotates the defense profile with weights/contributions."""
    annotated_profile = []
    raw_score = 0.0

    for sys in defense_profile:
        weight = get_system_weight(sys["system_name"])
        comp_count = sys["component_count"]
        completeness = sys["completeness"]

        contribution = weight * completeness * math.log(1 + comp_count)
        raw_score += contribution

        annotated_sys = sys.copy()
        annotated_sys["weight"] = weight
        annotated_sys["contribution"] = round(contribution, 2)
        annotated_profile.append(annotated_sys)

    # Reduced scalar from 5.0 to 2.2
    dbs = min(100.0, raw_score * 2.2)

    # Sort by contribution descending to get dominant systems
    annotated_profile.sort(key=lambda x: x["contribution"], reverse=True)
    dominant_systems = [sys["system_name"] for sys in annotated_profile[:3]]

    return {
        "DBS": round(dbs, 2),
        "system_count": len(defense_profile),
        "dominant_systems": dominant_systems,
        "profile": annotated_profile
    }

def pretty_print_dbs(result: dict, genome_name: str):
    """Prints a visually appealing summary using rich."""
    console = Console()
    table = Table(title=f"Defense Systems for {genome_name}")

    table.add_column("System", justify="left", style="cyan", no_wrap=True)
    table.add_column("Number", justify="right", style="magenta")
    table.add_column("Completeness", justify="right", style="green")
    table.add_column("Components", justify="right", style="blue")
    table.add_column("Weight", justify="right", style="yellow")
    table.add_column("Contribution", justify="right")

    for sys in result["profile"]:
        contrib = sys["contribution"]
        # Color code contribution
        if contrib >= 10:
            contrib_str = f"[red]{contrib}[/red]"
        elif contrib >= 5:
            contrib_str = f"[yellow]{contrib}[/yellow]"
        else:
            contrib_str = f"[green]{contrib}[/green]"

        table.add_row(
            sys["system_name"],
            str(sys["system_number"]),
            f"{sys['completeness']:.2f}",
            str(sys["component_count"]),
            f"{sys['weight']:.1f}",
            contrib_str
        )

    console.print(table)
    console.print(f"[bold]Defense Burden Score (DBS): {result['DBS']}/100[/bold]")
    console.print(f"Dominant systems: {', '.join(result['dominant_systems'])}\n")
