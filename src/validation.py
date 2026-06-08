import os
import sys
import json
import pandas as pd
os.environ.pop('MPLBACKEND', None)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix, matthews_corrcoef

# Add src directory to path
sys.path.append(os.path.join(os.path.dirname(__file__)))
from dbs_calculator import parse_padloc_output, calculate_dbs
from matching_engine import score_phage

GROUND_TRUTH = [
    ("phiKZ",  "AE004091", "PAO1",   1),
    ("phiKZ",  "NC_011770","PA7",    0),
    ("LUZ7",   "AE004091", "PAO1",   1),
    ("LUZ7",   "NC_011770","PA7",    0),
    ("LUZ24",  "AE004091", "PAO1",   0),
    ("PAK_P1", "AE004091", "PAO1",   1),
    ("PAK_P1", "NC_011770","PA7",    0),
    ("phiKZ",  "CP012001", "Kpn",    1),
    ("LUZ24",  "CP012001", "Kpn",    0),
    ("Ab31",   "CP012001", "Kpn",    1)
]

DECISION_THRESHOLD = 48.0

def main():
    os.makedirs("phagecompass/results", exist_ok=True)

    with open("phagecompass/results/phage_cdps.json", "r") as f:
        cdps = json.load(f)

    cdp_map = {c["phage_name"]: c for c in cdps}

    results = []
    for phage, acc, bac_name, observed in GROUND_TRUTH:
        csv_path = f"phagecompass/data/padloc_output/{acc}.fasta_padloc.csv"
        profile = parse_padloc_output(csv_path)
        dbs_res = calculate_dbs(profile)

        phage_cdp = cdp_map.get(phage)
        if not phage_cdp:
            print(f"Warning: {phage} CDP not found.")
            continue

        score_res = score_phage(dbs_res["profile"], phage_cdp)
        pc_score = score_res["phagecompass_score"]
        predicted = 1 if pc_score >= DECISION_THRESHOLD else 0

        results.append({
            "phage": phage,
            "bacteria": bac_name,
            "accession": acc,
            "pc_score": pc_score,
            "observed": observed,
            "predicted": predicted,
            "correct": bool(observed == predicted)
        })

    df = pd.DataFrame(results)

    y_obs = df["observed"]
    y_pred = df["predicted"]
    y_score = df["pc_score"]

    acc_score = accuracy_score(y_obs, y_pred)
    try:
        auc = roc_auc_score(y_obs, y_score)
    except ValueError:
        auc = 0.5
    mcc = matthews_corrcoef(y_obs, y_pred)
    cm = confusion_matrix(y_obs, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel() if len(cm.ravel()) == 4 else (0,0,0,0)

    print("=== VALIDATION METRICS ===")
    print(f"Accuracy:  {sum(df['correct'])}/{len(df)} correct ({acc_score*100:.1f}%)")
    print(f"AUC-ROC:   {auc:.2f}")
    print(f"MCC:       {mcc:.2f}")
    print("Confusion Matrix:")
    print(f"  True Positives:  {tp} (correctly predicted sensitive)")
    print(f"  True Negatives:  {tn} (correctly predicted resistant)")
    print(f"  False Positives: {fp} (predicted sensitive, was resistant)")
    print(f"  False Negatives: {fn} (predicted resistant, was sensitive)")
    print("\n=== PAIR RESULTS ===")
    print(df[["phage", "bacteria", "pc_score", "observed", "predicted", "correct"]].to_string(index=False))

    # ROC Curve
    from sklearn.metrics import roc_curve
    fpr, tpr, thresholds = roc_curve(y_obs, y_score)
    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, color="steelblue", lw=2, label=f"PhageCompass (AUC = {auc:.2f})")
    plt.plot([0,1],[0,1], color="gray", linestyle="--", label="Random")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("PhageCompass ROC Curve\nPhage-Host Infectivity Prediction")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig("phagecompass/results/roc_curve.png", dpi=300)
    plt.close()
    print("\nROC curve saved to phagecompass/results/roc_curve.png")

    # Score Distribution
    plt.figure(figsize=(8, 5))
    sensitive_scores = df[df["observed"] == 1]["pc_score"]
    resistant_scores = df[df["observed"] == 0]["pc_score"]

    plt.hist(sensitive_scores, bins=8, alpha=0.6, color="steelblue", label="Sensitive (infects)")
    plt.hist(resistant_scores, bins=8, alpha=0.6, color="firebrick", label="Resistant")
    plt.axvline(x=DECISION_THRESHOLD, color="black", linestyle="--", label=f"Decision threshold ({DECISION_THRESHOLD})")
    plt.xlabel("PhageCompass Score")
    plt.ylabel("Count")
    plt.title("Score Distribution by Infectivity Outcome")
    plt.legend()
    plt.tight_layout()
    plt.savefig("phagecompass/results/score_distribution.png", dpi=300)
    plt.close()
    print("Score distribution saved to phagecompass/results/score_distribution.png")

    report = {
        "n_pairs": len(df),
        "accuracy": round(acc_score, 2),
        "auc_roc": round(auc, 2),
        "mcc": round(mcc, 2),
        "threshold": DECISION_THRESHOLD,
        "confusion_matrix": {
            "TP": int(tp), "TN": int(tn), "FP": int(fp), "FN": int(fn)
        },
        "pair_results": df.to_dict(orient="records")
    }

    with open("phagecompass/results/validation_report.json", "w") as f:
        json.dump(report, f, indent=4)
    print("Validation report saved to phagecompass/results/validation_report.json")

if __name__ == "__main__":
    main()
