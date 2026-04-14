
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))
from reduced_polymer_model import ReducedPolymerFlood1D

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT = os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)

cases = [
    ("driza1_base", "water"),
    ("driza1_base", "polymer"),
    ("driza1_base", "hybrid"),
    ("driza2_base", "water"),
    ("driza2_base", "polymer"),
    ("patos_base", "polymer"),
]

plt.figure(figsize=(10,6))
for preset, mode in cases:
    params = ReducedPolymerFlood1D.load_preset(preset)
    res = ReducedPolymerFlood1D(params).run(mode=mode, n_steps=360)
    plt.plot(res.time_h, res.rf_percent, label=f"{preset}:{mode}")
plt.xlabel("Time [h]")
plt.ylabel("RF [%]")
plt.title("Demo comparison across presets and schedules")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "demo_rf_overview.png"), dpi=180)
print("Wrote demo figure.")
