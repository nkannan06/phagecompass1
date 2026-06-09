from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import json, os, glob
from src.fetch_genome import fetch_genome
from src.run_padloc import run_padloc
from src.dbs_calculator import parse_padloc_output, calculate_dbs
from src.matching_engine import match_panel

app = FastAPI(
    title="PhageCompass",
    description="Defense-aware phage therapy matching engine",
    version="1.0.0"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"]
)

class QueryRequest(BaseModel):
    accession: str
    name: str = None

@app.get("/")
def root():
    return {
        "tool": "PhageCompass",
        "version": "1.0.0",
        "description": "Defense-aware phage therapy matching engine",
        "usage": "POST /analyze with accession number",
        "docs": "/docs"
    }

@app.get("/health")
def health():
    return {"status": "ok"}

@app.post("/analyze")
def analyze(query: QueryRequest):
    try:
        name = query.name or query.accession
        os.makedirs("data/genomes", exist_ok=True)
        os.makedirs(f"data/padloc_output/{name}", exist_ok=True)

        fasta = fetch_genome(query.accession, "data/genomes")
        csv_path = run_padloc(fasta, f"data/padloc_output/{name}")
        profile = parse_padloc_output(csv_path)
        dbs_result = calculate_dbs(profile)

        with open("results/phage_cdps.json") as f:
            cdps = json.load(f)

        rankings = match_panel(dbs_result, cdps)

        return {
            "accession": query.accession,
            "name": name,
            "DBS": dbs_result["DBS"],
            "system_count": dbs_result["system_count"],
            "dominant_systems": dbs_result["dominant_systems"],
            "phage_rankings": rankings,
            "top_recommendation": rankings[0]
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/panel")
def get_panel():
    with open("results/phage_cdps.json") as f:
        return json.load(f)