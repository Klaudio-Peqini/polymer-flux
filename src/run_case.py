
from __future__ import annotations
import argparse, csv, json, os
import matplotlib.pyplot as plt
from reduced_polymer_model import ReducedPolymerFlood1D

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", required=True)
    ap.add_argument("--mode", required=True, choices=["water","polymer","hybrid"])
    ap.add_argument("--steps", type=int, default=320)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.dirname(__file__)), "outputs"))
    args = ap.parse_args()

    params = ReducedPolymerFlood1D.load_preset(args.preset)
    model = ReducedPolymerFlood1D(params)
    result = model.run(mode=args.mode, n_steps=args.steps)

    os.makedirs(args.outdir, exist_ok=True)
    stem = f"{args.preset}_{args.mode}"

    with open(os.path.join(args.outdir, stem + ".csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time_h","rf_percent","cum_oil_ml","delta_p_bar"])
        for row in zip(result.time_h, result.rf_percent, result.cum_oil_ml, result.delta_p_bar):
            w.writerow(row)

    with open(os.path.join(args.outdir, stem + "_summary.json"), "w", encoding="utf-8") as f:
        json.dump({
            "preset": args.preset,
            "mode": args.mode,
            "final_rf_percent": float(result.rf_percent[-1]),
            "final_cum_oil_ml": float(result.cum_oil_ml[-1]),
            "final_delta_p_bar": float(result.delta_p_bar[-1]),
            "mean_ads_ug_g": float(result.ads.mean()),
            "mean_cp_ppm": float(result.cp.mean()),
        }, f, indent=2)

    plt.figure(figsize=(11,7))
    plt.subplot(2,2,1)
    plt.plot(result.time_h, result.rf_percent)
    plt.xlabel("Time [h]"); plt.ylabel("RF [%]"); plt.title("Recovery factor")
    plt.subplot(2,2,2)
    plt.plot(result.time_h, result.cum_oil_ml)
    plt.xlabel("Time [h]"); plt.ylabel("Cumulative oil [ml]"); plt.title("Cumulative oil")
    plt.subplot(2,2,3)
    plt.plot(model.x_cm, result.sw, label="Sw")
    plt.plot(model.x_cm, result.cp / max(model.c_inj_ppm,1e-12), label="c/cinj")
    plt.xlabel("x [cm]"); plt.title("Final profiles"); plt.legend()
    plt.subplot(2,2,4)
    plt.plot(model.x_cm, result.ads, label="Adsorption [ug/g]")
    plt.plot(result.time_h, result.delta_p_bar, label="dP [bar]")
    plt.title("Adsorption + pressure trend"); plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, stem + ".png"), dpi=180)
    print("Done.")

if __name__ == "__main__":
    main()
