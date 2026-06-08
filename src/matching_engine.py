import json
import sys
import os
from rich.console import Console
from rich.table import Table

COUNTER_MAP = {
    "RM_type_I":    ["anti-RM"],
    "RM_type_II":   ["anti-RM"],
    "RM_type_III":  ["anti-RM"],
    "cbass":        ["anti-CBASS"],
    "CBASS":        ["anti-CBASS"],
    "CRISPR":       ["anti-CRISPR"],
    "cas_":         ["anti-CRISPR"],
    "Thoeris":      ["anti-Thoeris"],
    "retron":       ["anti-Retron"],
    "Lamassu":      ["anti-Lamassu"],
    "Zorya":        ["anti-RM"],      # Zorya is partially countered by anti-RM machinery
    "SoFic":        [],               # no known counter
    "Helicase":     [],               # no known counter
    "PDC":          [],               # candidate systems, counters unknown
    "PD-":          [],               # same
    "DMS":          [],               # same
    "Tiamat":       [],               # same
    "PrrC":         ["anti-RM"],      # PrrC is partially inhibited by anti-RM proteins
    "DEFAULT":      []
}

def get_counters(system_name: str) -> list[str]:
    system_name_lower = system_name.lower()
    
    if system_name_lower.startswith("cas"):
        return COUNTER_MAP["cas_"]
        
    for key, counters in COUNTER_MAP.items():
        if key != "DEFAULT" and system_name_lower.startswith(key.lower()):
            return counters
    return COUNTER_MAP["DEFAULT"]

def get_recommendation(score: float) -> str:
    if score >= 75:
        return "High viability — prioritize for experimental testing"
    elif score >= 55:
        return "Moderate viability — consider as part of cocktail"
    elif score >= 35:
        return "Low viability — resistance likely, use only if no better option"
    else:
        return "Not recommended — significant resistance predicted"

def score_phage(bacterial_profile: list[dict], phage_cdp: dict) -> dict:
    base_score = 50.0
    system_breakdown = []
    
    for S in bacterial_profile:
        potency = S["weight"] * S["completeness"]
        counters = get_counters(S["system_name"])
        
        if not counters:
            base_score -= potency * 0.3
            status = "no_known_counter"
            active_counter = None
        elif any(c in phage_cdp.get("anti_defense_systems", []) for c in counters):
            base_score += potency * 0.8
            status = "countered"
            active_counter = next(c for c in counters if c in phage_cdp["anti_defense_systems"])
        else:
            base_score -= potency * 0.5
            status = "uncountered"
            active_counter = None
            
        system_breakdown.append({
            "system": S["system_name"],
            "potency": round(potency, 2),
            "status": status,
            "counter": active_counter
        })
        
    phagecompass_score = max(0.0, min(100.0, base_score))
    
    return {
        "phage_name": phage_cdp["phage_name"],
        "phagecompass_score": round(phagecompass_score, 2),
        "recommendation": get_recommendation(phagecompass_score),
        "system_breakdown": system_breakdown
    }

def match_panel(dbs_result: dict, phage_cdps: list[dict]) -> list[dict]:
    rankings = []
    for cdp in phage_cdps:
        score_result = score_phage(dbs_result["profile"], cdp)
        rankings.append(score_result)
        
    rankings.sort(key=lambda x: x["phagecompass_score"], reverse=True)
    return rankings
