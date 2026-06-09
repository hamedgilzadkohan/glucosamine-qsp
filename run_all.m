%% =============================================================================
%% run_all.m
%% Master Script - Reproduces All Analyses
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%%
%% This script runs the complete analysis pipeline in the correct order.
%% A random seed of rng(42) is used in each script for reproducibility.
%%
%% Steps 14-17 are reviewer-response analyses:
%%   Step 14    - equal-bioavailability sensitivity (Reviewer 1)
%%   Steps 15-17 - wide-range IC50,deg sensitivity / re-calibration (Reviewer 2)
%%
%% Requirements: MATLAB R2023b+ with Optimization Toolbox
%% Estimated runtime: ~30-60 minutes (bootstrap is the bottleneck)
%%
%% =============================================================================
fprintf('=============================================================\n');
fprintf('GLUCOSAMINE QSP MODEL - FULL ANALYSIS PIPELINE\n');
fprintf('=============================================================\n\n');
tic;
%% Step 1: Model Structure and Parameter Definitions
fprintf('Step 1/17: Model structure...\n');
[param_info, default_params, param_bounds] = M01_model_structure();
%% Step 2: Model Equations Test
fprintf('\nStep 2/17: Equations test...\n');
M02_model_equations;
%% Step 3: Calibration Data Preparation
fprintf('\nStep 3/17: Calibration data...\n');
calibration_data = M03_calibration_data();
%% Step 4: Set Calibrated Parameters (matched to R manuscript)
fprintf('\nStep 4/17: Set R parameters...\n');
calibration_results = M04_set_R_parameters();
%% Step 5: Sensitivity Analysis
fprintf('\nStep 5/17: Sensitivity analysis...\n');
sensitivity_results = M05_sensitivity_analysis();
%% Step 6: Profile Likelihood
fprintf('\nStep 6/17: Profile likelihood...\n');
profile_results = M06_profile_likelihood();
%% Step 7: Bootstrap Confidence Intervals
fprintf('\nStep 7/17: Bootstrap CIs (this may take ~20-40 min)...\n');
bootstrap_results = M07_bootstrap_ci();
%% Step 8: Residual Diagnostics
fprintf('\nStep 8/17: Residual diagnostics...\n');
diagnostic_results = M08_residual_diagnostics();
%% Step 9: Virtual Population Simulation
fprintf('\nStep 9/17: Virtual population...\n');
vpop_results = M09_virtual_population();
%% Step 10: External Validation
fprintf('\nStep 10/17: External validation...\n');
validation_results = M10_external_validation();
%% Step 11: Translational Simulations
fprintf('\nStep 11/17: Translational simulations...\n');
trans_results = M11_translational_simulations();
%% Step 12: Publication Figures
fprintf('\nStep 12/17: Generating main figures...\n');
M12_figures;
%% Step 13: Supplementary Figures
fprintf('\nStep 13/17: Generating supplementary figures...\n');
M13_supplementary_figures;
%% Step 14: Equal-Bioavailability Sensitivity Analysis (Reviewer Response)
fprintf('\nStep 14/17: Equal-bioavailability sensitivity analysis...\n');
equal_ba_results = M14_equal_bioavailability_sensitivity();
%% Step 15: Wide-Range IC50,deg Re-Calibration (Reviewer 2, Points 7-8)
%% Re-fits the 5 free parameters at each IC50,deg; writes
%% model_output/ic50_recalibration_results.mat (consumed by M17).
fprintf('\nStep 15/17: IC50 wide-range re-calibration...\n');
ic50_recal = M15_IC50_recalibration();
%% Step 16: Forward IC50,deg Sweep (Reviewer 2)
%% Holds calibrated parameters fixed; writes
%% model_output/ic50_forward_sweep_results.mat (consumed by M17).
fprintf('\nStep 16/17: IC50 forward sweep...\n');
fwd = M16_IC50_forward_sweep();
%% Step 17: Combined 4-Panel Figure S6 (Reviewer 2)
%% Reads the M15 and M16 results above; must run after both.
fprintf('\nStep 17/17: Building combined Figure S6...\n');
M17_make_figS6_combined();
%% Done
elapsed = toc;
fprintf('\n=============================================================\n');
fprintf('COMPLETE PIPELINE FINISHED\n');
fprintf('Total runtime: %.1f minutes\n', elapsed/60);
fprintf('=============================================================\n');
fprintf('\nOutputs:\n');
fprintf('  data_clean/    - Processed calibration data\n');
fprintf('  model_output/  - All analysis results (.mat)\n');
fprintf('  figures/       - All publication figures (.png, .fig)\n');
fprintf('=============================================================\n');
