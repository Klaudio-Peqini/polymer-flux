
from __future__ import annotations
import json, os
from dataclasses import dataclass
from typing import Dict, Tuple
import numpy as np

@dataclass
class PolymerFloodResult:
    time_h: np.ndarray
    rf_percent: np.ndarray
    cum_oil_ml: np.ndarray
    delta_p_bar: np.ndarray
    sw: np.ndarray
    cp: np.ndarray
    ads: np.ndarray
    case_name: str
    schedule: str
    meta: Dict[str, float]

class ReducedPolymerFlood1D:
    def __init__(self, params: Dict[str, float]):
        self.p = dict(params)
        self.n = int(self.p["n_cells"])
        self.L_cm = float(self.p["length_cm"])
        self.A_cm2 = float(self.p["area_cm2"])
        self.dx_cm = self.L_cm / self.n
        self.x_cm = np.linspace(self.dx_cm/2.0, self.L_cm-self.dx_cm/2.0, self.n)

        self.phi = float(self.p["phi"])
        self.k_darcy = float(self.p["k_darcy"])
        self.Swi = float(self.p["Swi"])
        self.Sor = float(self.p["Sor"])
        self.mu_w_cp = float(self.p["mu_w_cp"])
        self.mu_o_cp = float(self.p["mu_o_cp"])
        self.mu_p_ref_cp = float(self.p["mu_p_ref_cp"])
        self.c_inj_ppm = float(self.p["c_inj_ppm"])
        self.gamma_max = float(self.p["gamma_max_ug_g"])
        self.langmuir_K = float(self.p["langmuir_K_ppm"])
        self.ipv = float(self.p["ipv"])
        self.rrf_max = float(self.p["rrf_max"])
        self.krw0 = float(self.p["krw0"])
        self.kro0 = float(self.p["kro0"])
        self.nw = float(self.p["nw"])
        self.no = float(self.p["no"])
        self.rate_ml_h = float(self.p["inj_rate_ml_h"])
        self.t_end_h = float(self.p["t_end_h"])
        self.hybrid_switch_h = float(self.p["hybrid_switch_h"])
        self.oip_ml = float(self.p["oil_in_place_ml"])

        self.mix_exp = 1.15
        self.effective_pore_volume_ml = self.phi * self.A_cm2 * self.L_cm
        self.front_velocity_scale = self.rate_ml_h / max(self.effective_pore_volume_ml, 1e-12)

    @staticmethod
    def load_preset(name: str, path: str | None = None) -> Dict[str, float]:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "configs", "presets.json")
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data["presets"][name]

    def relperm(self, sw: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        se = (sw - self.Swi) / max(1.0 - self.Swi - self.Sor, 1e-12)
        se = np.clip(se, 0.0, 1.0)
        krw = self.krw0 * se**self.nw
        kro = self.kro0 * (1.0 - se)**self.no
        return krw, kro

    def adsorption_langmuir(self, c_ppm: np.ndarray) -> np.ndarray:
        c = np.maximum(c_ppm, 0.0)
        return self.gamma_max * c / (self.langmuir_K + c + 1e-12)

    def mu_aqueous_cp(self, c_ppm: np.ndarray) -> np.ndarray:
        frac = np.clip(c_ppm / max(self.c_inj_ppm, 1e-12), 0.0, 2.0)
        return self.mu_w_cp + (self.mu_p_ref_cp - self.mu_w_cp) * frac**self.mix_exp

    def rrf(self, c_ppm: np.ndarray, ads_ug_g: np.ndarray) -> np.ndarray:
        conc_term = np.clip(c_ppm / max(self.c_inj_ppm, 1e-12), 0.0, 1.5)
        ads_term = np.clip(ads_ug_g / max(self.gamma_max, 1e-12), 0.0, 1.0)
        return 1.0 + (self.rrf_max - 1.0) * (0.45 * conc_term + 0.55 * ads_term)

    def schedule_concentration(self, t_h: float, mode: str) -> float:
        if mode == "water":
            return 0.0
        if mode == "polymer":
            return self.c_inj_ppm
        if mode == "hybrid":
            return 0.0 if t_h < self.hybrid_switch_h else self.c_inj_ppm
        raise ValueError(mode)

    def _rf_limit(self, exposure: float, mean_ads: float) -> float:
        base = 0.245 + 0.25 * (self.kro0 - 0.75) + 0.03 * (self.krw0 - 0.38)
        poly = 0.30 * np.clip(exposure, 0.0, 1.0)
        poly += 0.04 * np.clip((self.rrf_max - 1.0) / 0.4, 0.0, 1.5)
        penalty = 0.03 * np.clip(mean_ads / max(self.gamma_max,1e-12), 0.0, 1.0)
        return float(np.clip(base + poly - penalty, 0.12, 0.75))

    def run(self, mode: str = "water", n_steps: int = 320) -> PolymerFloodResult:
        dt_h = self.t_end_h / n_steps
        sw = np.full(self.n, self.Swi, dtype=float)
        cp = np.zeros(self.n, dtype=float)
        ads = np.zeros(self.n, dtype=float)

        t_hist, rf_hist, oil_hist, dp_hist = [], [], [], []
        polymer_clock = 0.0

        for it in range(n_steps):
            t_h = (it + 1) * dt_h
            c_in = self.schedule_concentration(t_h, mode)
            polymer_clock += dt_h * (c_in / max(self.c_inj_ppm,1e-12))
            exposure = polymer_clock / max(t_h,1e-12)

            krw, kro = self.relperm(sw)
            muw = self.mu_aqueous_cp(cp)
            rrf = self.rrf(cp, ads)
            lamw = krw / (muw * rrf + 1e-12)
            lamo = kro / (self.mu_o_cp + 1e-12)
            fw = lamw / (lamw + lamo + 1e-12)

            v = self.front_velocity_scale * 0.20
            alpha_w = min(0.60, v * dt_h / max((1.0 - self.ipv) * self.n, 1e-12))
            alpha_c = min(0.75, 1.2 * v * dt_h / max((1.0 - self.ipv) * self.n, 1e-12))

            fw_face = np.empty(self.n + 1)
            fw_face[0] = 1.0
            fw_face[1:] = fw
            sw = np.clip(sw + alpha_w * (fw_face[:-1] - fw_face[1:]), self.Swi, 1.0 - self.Sor)

            cp_face = np.empty(self.n + 1)
            cp_face[0] = c_in
            cp_face[1:] = cp
            cp_adv = cp + alpha_c * (cp_face[:-1] - cp_face[1:])

            ads_eq = self.adsorption_langmuir(cp_adv)
            ads += 0.10 * (ads_eq - ads)
            ads_norm = np.clip(ads / max(self.gamma_max, 1e-12), 0.0, 1.0)
            cp = np.clip(cp_adv * (1.0 - 0.20 * ads_norm) * (1.0 - 0.35*self.ipv), 0.0, 1.25*self.c_inj_ppm)

            mean_ads = float(ads.mean())
            rf_lim = self._rf_limit(exposure, mean_ads)
            pv_inj = self.rate_ml_h * t_h / max(self.effective_pore_volume_ml, 1e-12)
            shape = 1.0 - np.exp(-0.17 * pv_inj)
            rf_target = rf_lim * shape
            cum_oil = min(self.oip_ml*rf_target, self.oip_ml*0.80)
            rf_percent = 100.0 * cum_oil / max(self.oip_ml, 1e-12)

            lam_total_mean = float(np.mean(lamw + lamo))
            mobility_factor = 1.0 / max(lam_total_mean, 1e-12)
            dp_bar = 0.004 * self.rate_ml_h * self.L_cm * mobility_factor / max(self.A_cm2 * self.k_darcy, 1e-12)

            t_hist.append(t_h); rf_hist.append(rf_percent); oil_hist.append(cum_oil); dp_hist.append(dp_bar)

        return PolymerFloodResult(
            time_h=np.array(t_hist), rf_percent=np.array(rf_hist), cum_oil_ml=np.array(oil_hist),
            delta_p_bar=np.array(dp_hist), sw=sw, cp=cp, ads=ads,
            case_name=self.p.get("name", "case"), schedule=mode,
            meta={"oip_ml": self.oip_ml, "rate_ml_h": self.rate_ml_h}
        )
