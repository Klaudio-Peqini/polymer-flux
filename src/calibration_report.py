
from __future__ import annotations
import csv, json, os, sys
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(__file__))
from reduced_polymer_model import ReducedPolymerFlood1D

ROOT = os.path.dirname(os.path.dirname(__file__))
CFG = os.path.join(ROOT, "configs", "presets.json")
OUT = os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)

with open(CFG, "r", encoding="utf-8") as f:
    data = json.load(f)

targets = data["paper_targets"]
cases = [
    ("driza1_base", "water", "driza1_water"),
    ("driza1_base", "polymer", "driza1_polymer"),
    ("driza1_base", "hybrid", "driza1_hybrid"),
    ("driza2_base", "water", "driza2_water"),
    ("driza2_base", "polymer", "driza2_polymer"),
]

rows = []
labels = []
target_rf = []
model_rf = []

for preset, mode, target_key in cases:
    model = ReducedPolymerFlood1D(data["presets"][preset])
    res = model.run(mode=mode, n_steps=360)
    trf = targets[target_key]["rf_percent"]
    toil = targets[target_key]["oil_ml"]
    mrf = float(res.rf_percent[-1])
    moil = float(res.cum_oil_ml[-1])
    rows.append({
        "preset": preset,
        "mode": mode,
        "target_rf_percent": trf,
        "model_rf_percent": mrf,
        "abs_rf_error": abs(mrf - trf),
        "target_oil_ml": toil,
        "model_oil_ml": moil,
        "abs_oil_error_ml": abs(moil - toil),
        "final_delta_p_bar": float(res.delta_p_bar[-1]),
        "mean_cp_ppm": float(res.cp.mean()),
        "mean_ads_ug_g": float(res.ads.mean()),
    })
    labels.append(f"{preset}\n{mode}")
    target_rf.append(trf)
    model_rf.append(mrf)

with open(os.path.join(OUT, "calibration_report.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

with open(os.path.join(OUT, "calibration_report.md"), "w", encoding="utf-8") as f:
    f.write("# Calibration report\n\n")
    f.write("| Preset | Mode | Target RF [%] | Model RF [%] | |RF error| | Target oil [ml] | Model oil [ml] | |Oil error| [ml] | Final dP [bar] |\n")
    f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
    for r in rows:
        f.write(f"| {r['preset']} | {r['mode']} | {r['target_rf_percent']:.2f} | {r['model_rf_percent']:.2f} | {r['abs_rf_error']:.2f} | {r['target_oil_ml']:.2f} | {r['model_oil_ml']:.2f} | {r['abs_oil_error_ml']:.2f} | {r['final_delta_p_bar']:.3f} |\n")

plt.figure(figsize=(10,5))
x = range(len(labels))
plt.plot(list(x), target_rf, marker="o", label="Paper target")
plt.plot(list(x), model_rf, marker="s", label="Reduced model")
plt.xticks(list(x), labels, rotation=25, ha="right")
plt.ylabel("Recovery factor [%]")
plt.title("Paper targets vs reduced-model predictions")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "calibration_comparison.png"), dpi=180)
print("Wrote calibration report.")
