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
%% Requirements: MATLAB R2023b+ with Optimization Toolbox
%% Estimated runtime: ~30-60 minutes (bootstrap is the bottleneck)
%%
%% =============================================================================

fprintf('=============================================================\n');
fprintf('GLUCOSAMINE QSP MODEL - FULL ANALYSIS PIPELINE\n');
fprintf('=============================================================\n\n');

tic;

%% Step 1: Model Structure and Parameter Definitions
fprintf('Step 1/14: Model structure...\n');
[param_info, default_params, param_bounds] = M01_model_structure();

%% Step 2: Model Equations Test
fprintf('\nStep 2/14: Equations test...\n');
M02_model_equations;

%% Step 3: Calibration Data Preparation
fprintf('\nStep 3/14: Calibration data...\n');
calibration_data = M03_calibration_data();

%% Step 4: Set Calibrated Parameters (matched to R manuscript)
fprintf('\nStep 4/14: Set R parameters...\n');
calibration_results = M04_set_R_parameters();

%% Step 5: Sensitivity Analysis
fprintf('\nStep 5/14: Sensitivity analysis...\n');
sensitivity_results = M05_sensitivity_analysis();

%% Step 6: Profile Likelihood
fprintf('\nStep 6/14: Profile likelihood...\n');
profile_results = M06_profile_likelihood();

%% Step 7: Bootstrap Confidence Intervals
fprintf('\nStep 7/14: Bootstrap CIs (this may take ~20-40 min)...\n');
bootstrap_results = M07_bootstrap_ci();

%% Step 8: Residual Diagnostics
fprintf('\nStep 8/14: Residual diagnostics...\n');
diagnostic_results = M08_residual_diagnostics();

%% Step 9: Virtual Population Simulation
fprintf('\nStep 9/14: Virtual population...\n');
vpop_results = M09_virtual_population();

%% Step 10: External Validation
fprintf('\nStep 10/14: External validation...\n');
validation_results = M10_external_validation();

%% Step 11: Translational Simulations
fprintf('\nStep 11/14: Translational simulations...\n');
trans_results = M11_translational_simulations();

%% Step 12: Publication Figures
fprintf('\nStep 12/14: Generating main figures...\n');
M12_figures;

%% Step 13: Supplementary Figures
fprintf('\nStep 13/14: Generating supplementary figures...\n');
M13_supplementary_figures;

%% Step 14: Equal-Bioavailability Sensitivity Analysis (Reviewer Response)
fprintf('\nStep 14/14: Equal-bioavailability sensitivity analysis...\n');
equal_ba_results = M14_equal_bioavailability_sensitivity();

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
