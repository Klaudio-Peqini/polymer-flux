# Polymer Flooding Reduced-Order Model                                                                                                                                            

---

## Overview

This repository provides a **research-grade reduced-order simulator** for polymer flooding in porous media.

It is designed to bridge the gap between:

- **Laboratory coreflood experiments (Driza 1 & Driza 2)**
- **Physical understanding of polymer EOR mechanisms**
- **Field-scale intuition (Patos-Marinza)**

---

## Objectives

This model enables:

- Rapid **screening of polymer injection strategies**
- Insight into **mobility control and sweep efficiency**
- Transparent **physics-driven interpretation**
- A **lightweight alternative** to full simulators like OPM Flow

---

## Physical Model Description

### Governing Variables

The model evolves:

- Pressure: $ p(x,t) $
- Water saturation: $ S_w(x,t) $
- Polymer concentration: $ c(x,t) $

---

### Polymer Transport

Polymer is transported in the aqueous phase:

$$
\frac{\partial}{\partial t}(\phi_{eff} S_w c) + \nabla \cdot (u_w c) = - \rho_r \frac{\partial \Gamma(c)}{\partial t}
$$

---

### Effective Porosity

$$
\phi_{eff} = \phi (1 - s_{ipv})
$$

---

### Adsorption (Langmuir)

$$
\Gamma(c) = \Gamma_{max} \frac{c}{K + c}
$$

---

### Effective Viscosity

$$
\mu_{w,eff}(c) = \mu_w + f(c)(\mu_{wp} - \mu_w)
$$

---

### Mobility

$$
\lambda_w = \frac{k_{rw}}{\mu_{w,eff}(c)}, \quad \lambda_o = \frac{k_{ro}}{\mu_o}
$$

---

### 🪨 Residual Resistance Factor

$$
k_{eff} = \frac{k}{R_k(c)}
$$

---

## Numerical Method

- Finite Volume Method (TPFA)
- Upwind flux discretization
- Semi-implicit pressure solve
- Explicit transport update (v2)

---

## Repository Structure

```
polymer_v2_package/
├── src/
│   ├── reduced_polymer_model.py
│   ├── run_case.py
│   ├── calibration_report.py
├── configs/
│   └── presets.json
├── examples/
│   └── run_demo.py
├── outputs/
├── notes/
├── README.md
└── requirements.txt
```

---

## Installation

```bash
pip install -r requirements.txt
```

---

## Running Simulations

```bash
python src/run_case.py --preset driza1_base --mode polymer
```

### Modes

- `water`
- `polymer`
- `hybrid`

---

## Preset Library

### Driza 1

| Preset | Description |
|-------|------------|
| driza1_base | Best estimate |
| driza1_conservative | Lower polymer efficiency |
| driza1_strong | Aggressive polymer |

---

### Driza 2

| Preset | Description |
|-------|------------|
| driza2_base | Validation case |

---

### Patos-Marinza

| Preset | Description |
|-------|------------|
| patos_base | Representative |
| patos_heavy_oil | Extreme viscosity |
| patos_high_perm | High permeability |

---

## 📊 Experimental Targets

| Case | RF |
|------|----|
| Driza 1 water | 25% |
| Driza 1 polymer | 53% |
| Driza 1 hybrid | 48% |
| Driza 2 water | 25.8% |
| Driza 2 polymer | 52% |

---

## Model Performance

✔ Captures correct ordering  
✔ Approximates magnitude  
✔ Reproduces polymer benefit  

---

## Extended Study Cases

### Polymer concentration scan
1000–3000 ppm

### Viscosity scan
10–50 cP

### Adsorption scan
50–200 µg/g

### Permeability reduction
0.7–0.95

---

## Outputs

- recovery.csv
- profiles.csv
- pressure.csv
- summary.json
- plots.png

---

## Calibration Workflow

1. Fit waterflood
2. Match breakthrough
3. Enable polymer
4. Tune adsorption & viscosity
5. Validate hybrid

---

## Limitations

- 1D only
- No heterogeneity
- No gas phase
- Simplified numerics

---

## Future Work

- Fully implicit solver
- AD Jacobian
- 2D/3D extension
- Shear-thinning rheology
- OPM coupling

---

## Discussion Points

- Adsorption vs mobility control
- Polymer efficiency limits
- Scaling lab → field
- Numerical vs physical uncertainty

---

## Final Insight

This framework provides:

**Understanding > black-box simulation**

It is ideal for:
- teaching
- research prototyping
- parameter sensitivity studies

---

## Contact / Contribution

Open for extension toward:
- HPC integration
- multi-phase models
- machine learning calibration

