
# Notes for the professor

## What this package tries to accomplish
This package turns the earlier toy script into a more paper-aligned reduced model with:

- polymer concentration transport,
- Langmuir adsorption,
- inaccessible pore volume,
- concentration-dependent aqueous viscosity,
- residual resistance factor,
- water / polymer / hybrid schedules,
- calibration comparison against the five paper targets.

## What it does not claim
It is not a substitute for:
- OPM Flow,
- fully implicit black-oil simulation,
- 3D core geometry,
- automatic differentiation / exact Jacobian Newton solves.

## Suggested next scientific step
Use this reduced model to:
1. calibrate waterflood first,
2. calibrate polymer uplift second,
3. compare sensitivity of adsorption, viscosity, and mobility reduction,
4. only then decide which mechanisms deserve full OPM-scale study.
