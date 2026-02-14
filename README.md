# Glucosamine QSP Model for Knee Osteoarthritis

A Quantitative Systems Pharmacology (QSP) model for evaluating the efficacy of glucosamine sulfate in slowing structural progression and reducing pain in knee osteoarthritis.

## Overview

This repository contains the MATLAB implementation of a mechanistic QSP model that integrates:

- **Pharmacokinetics:** Steady-state synovial fluid concentrations for sulfate and HCl formulations
- **Cartilage dynamics:** GAG synthesis/degradation with dual drug effect mechanisms (anti-catabolic + anabolic)
- **Clinical endpoints:** Joint space width (JSW) narrowing and WOMAC pain scores
- **Placebo modeling:** Exponential onset placebo response

The model was calibrated against published clinical trial data (Reginster 2001, Pavelka 2002, GAIT 2006) and externally validated against two independent trials (GUIDE 2007, MOVES 2015).

## Repository Structure

```
glucosamine-qsp/
├── run_all.m                    # Master script — runs the full pipeline
├── src/                         # All source code
│   ├── simulate_glucosamine.m   # Core simulation function
│   ├── M01_model_structure.m    # Parameter definitions and bounds
│   ├── M02_model_equations.m    # Model equations test
│   ├── M03_calibration_data.m   # Calibration data assembly
│   ├── M04_set_R_parameters.m   # Set calibrated parameter values
│   ├── M05_sensitivity_analysis.m  # Local sensitivity analysis
│   ├── M06_profile_likelihood.m    # Profile likelihood for identifiability
│   ├── M07_bootstrap_ci.m         # Bootstrap confidence intervals
│   ├── M08_residual_diagnostics.m  # Goodness-of-fit assessment
│   ├── M09_virtual_population.m    # Virtual population (N=500)
│   ├── M10_external_validation.m   # GUIDE 2007 & MOVES 2015 validation
│   ├── M11_translational_simulations.m  # Dose-response, bioavailability, etc.
│   ├── M12_figures.m               # Main figures (Figures 1-5)
│   └── M13_supplementary_figures.m # Supplementary figures (S1-S5)
├── data_clean/                  # Generated: processed calibration data
├── model_output/                # Generated: saved model results (.mat)
├── figures/                     # Generated: publication figures
├── LICENSE                      # MIT License
└── README.md                    # This file
```

## System Requirements

### Software

- **MATLAB R2019b** or later (tested on R2023b)
- **Optimization Toolbox** — required for `fmincon` in bootstrap analysis (M07)
- **Statistics and Machine Learning Toolbox** — required for `prctile`, `ksdensity`, `normpdf`

### Hardware

- No special hardware required
- Bootstrap analysis (M07, N=1000 replicates) may take 2–8 hours depending on CPU
- All other steps complete in minutes

## Installation and Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/hamedgilzadkohan/glucosamine-qsp.git
   cd glucosamine-qsp
   ```

2. **Open MATLAB** and navigate to the repository root folder.

3. **Verify toolboxes** (in MATLAB Command Window):
   ```matlab
   ver('optim')    % Should show Optimization Toolbox
   ver('stats')    % Should show Statistics and Machine Learning Toolbox
   ```

## Running the Pipeline

### Full pipeline (all steps)

```matlab
run('run_all.m')
```

This executes all 13 analysis steps sequentially. Each step saves its outputs to `data_clean/` or `model_output/`, so downstream steps can load intermediate results.

### Individual steps

Each script can be run independently once its dependencies have been generated:

```matlab
addpath('src');

% Example: run calibration data assembly
calibration_data = M03_calibration_data();

% Example: generate figures (requires model_output/ to exist)
run('src/M12_figures.m');
```

### Dependency chain

| Step | Script | Depends on |
|------|--------|------------|
| 1 | M01 | — |
| 2 | M02 | M01, simulate_glucosamine |
| 3 | M03 | — |
| 4 | M04 | M01, simulate_glucosamine |
| 5 | M05 | M04 output |
| 6 | M06 | M04, M03 output |
| 7 | M07 | M04, M03 output |
| 8 | M08 | M04, M03 output |
| 9 | M09 | M04 output |
| 10 | M10 | M09 output |
| 11 | M11 | M09 output |
| 12 | M12 | M04, M03, M05, M09, M10, M11 output |
| 13 | M13 | M04, M03, M07 output |

## Expected Outputs

### Key numerical results

| Endpoint | Value |
|----------|-------|
| 3-year JSW treatment effect (GS 1500 mg) | ~0.24 mm |
| 3-year placebo JSW change | ~−0.21 mm |
| Week 24 WOMAC pain (placebo, HCl trial) | ~29.7 points |

### Generated figures

**Main figures** (in `figures/`):
- `Fig1_calibration.png` — Model calibration (JSW, GAG, Pain)
- `Fig2_sensitivity.png` — Tornado sensitivity analysis
- `Fig3_virtual_population.png` — Virtual population with 90% prediction intervals
- `Fig4_validation.png` — External validation (GUIDE 2007, MOVES 2015)
- `Fig5_translational.png` — Dose-response, bioavailability, duration, stratification

**Supplementary figures:**
- `FigS1_profile_likelihood.png` — Parameter identifiability
- `FigS2_residual_diagnostics.png` — Residuals and goodness-of-fit
- `FigS3_bootstrap_distributions.png` — Bootstrap parameter distributions
- `FigS4_bioavailability_impact.png` — Bioavailability sensitivity
- `FigS5_kpl_sensitivity.png` — Placebo rate constant sensitivity

### Saved data

All intermediate results are saved as `.mat` files in `model_output/`:
- `calibrated_params.mat`, `calibration_results.mat`
- `sensitivity_analysis_results.mat`
- `profile_likelihood_results.mat`
- `bootstrap_results.mat`
- `residual_diagnostics.mat`
- `virtual_population_results.mat`
- `validation_results.mat`
- `translational_simulations.mat`

## Key Model Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| k_syn | 0.1753 /year | GAG synthesis rate |
| k_deg | 0.2073 /year | GAG degradation rate |
| Imax_deg | 1.0 | Maximum degradation inhibition |
| Pain₀ | 46.65 | Baseline WOMAC pain (0–100) |
| Pmax | 17.62 | Maximum placebo response |
| kpl | 0.02 /day | Placebo effect rate constant |
| F_sulfate | 0.22 | Bioavailability (crystalline sulfate) |
| F_HCl | 0.11 | Bioavailability (HCl formulation) |
| Emax_syn | 0.15 | Maximum synthesis stimulation |
| EC50_syn | 2.0 ug/mL | EC50 for synthesis stimulation |

## Citation

If you use this code, please cite:

> Soheili M, Gilzad Kohan H. A Quantitative Systems Pharmacology Model for Glucosamine Sulfate in Knee Osteoarthritis: 
Mechanistic Integration of Pharmacokinetics, Cartilage Homeostasis, and Clinical Outcomes.
*Journal of Pharmacy and Pharmaceutical Sciences*, 2026. DOI: pending

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

For questions about the model or code, please open an issue on this repository.
