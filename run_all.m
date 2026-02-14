%% =============================================================================
%% run_all.m
%% Master Script: Execute Full Glucosamine QSP Model Pipeline
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%%
%% USAGE:
%%   1. Open MATLAB
%%   2. Navigate to the repository root folder
%%   3. Run:  run('run_all.m')
%%
%% PIPELINE:
%%   Step 1:  M01 - Parameter definitions
%%   Step 2:  M02 - Model equations test
%%   Step 3:  M03 - Calibration data assembly
%%   Step 4:  M04 - Set calibrated parameters
%%   Step 5:  M05 - Local sensitivity analysis
%%   Step 6:  M06 - Profile likelihood analysis
%%   Step 7:  M07 - Bootstrap confidence intervals  ** SLOW (2-8 hours) **
%%   Step 8:  M08 - Residual diagnostics
%%   Step 9:  M09 - Virtual population simulation (N=500)
%%   Step 10: M10 - External validation (GUIDE 2007, MOVES 2015)
%%   Step 11: M11 - Translational simulations
%%   Step 12: M12 - Main figures (Figures 1-5)
%%   Step 13: M13 - Supplementary figures (Figures S1-S5)
%%
%% REQUIREMENTS:
%%   MATLAB R2019b or later
%%   Optimization Toolbox (for M07 bootstrap)
%%   Statistics and Machine Learning Toolbox (for M07, M09, M13)
%%
%% =============================================================================

clear; clc; close all;

%% =========================================================================
%% SET PATH - Automatically detect repository root
%% =========================================================================

% Detect the directory containing this script
repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root)
    repo_root = pwd;
end
cd(repo_root);

% Add source directory to path
addpath(fullfile(repo_root, 'src'));

fprintf('=============================================================\n');
fprintf('GLUCOSAMINE QSP MODEL - FULL PIPELINE\n');
fprintf('Working directory: %s\n', pwd);
fprintf('Started: %s\n', string(datetime('now')));
fprintf('=============================================================\n\n');

% Create output directories
if ~exist('data_clean', 'dir'),   mkdir('data_clean'); end
if ~exist('model_output', 'dir'), mkdir('model_output'); end
if ~exist('figures', 'dir'),      mkdir('figures'); end

%% =========================================================================
%% STEP 1: Verify all source files are present
%% =========================================================================

fprintf('--- STEP 1: Verifying source files ---\n');

required_files = {
    'simulate_glucosamine.m', ...
    'M01_model_structure.m', ...
    'M02_model_equations.m', ...
    'M03_calibration_data.m', ...
    'M04_set_R_parameters.m', ...
    'M05_sensitivity_analysis.m', ...
    'M06_profile_likelihood.m', ...
    'M07_bootstrap_ci.m', ...
    'M08_residual_diagnostics.m', ...
    'M09_virtual_population.m', ...
    'M10_external_validation.m', ...
    'M11_translational_simulations.m', ...
    'M12_figures.m', ...
    'M13_supplementary_figures.m'
};

all_found = true;
for i = 1:length(required_files)
    if exist(required_files{i}, 'file')
        fprintf('  [OK] %s\n', required_files{i});
    else
        fprintf('  [MISSING] %s\n', required_files{i});
        all_found = false;
    end
end

if ~all_found
    error('Missing source files. Ensure all .m files are in the src/ directory.');
end
fprintf('  All files found.\n\n');

%% =========================================================================
%% STEP 2: Load parameter definitions (M01)
%% =========================================================================

fprintf('--- STEP 2: Loading parameter definitions (M01) ---\n');

[param_info, default_params, param_bounds] = M01_model_structure();
fprintf('\n');

%% =========================================================================
%% STEP 3: Model equations test (M02)
%% =========================================================================

fprintf('--- STEP 3: Model equations test (M02) ---\n');

run(fullfile('src', 'M02_model_equations.m'));
fprintf('\n');

%% =========================================================================
%% STEP 4: Calibration data assembly (M03)
%% =========================================================================

fprintf('--- STEP 4: Calibration data assembly (M03) ---\n');

calibration_data = M03_calibration_data();
fprintf('\n');

%% =========================================================================
%% STEP 5: Set calibrated parameters (M04)
%% =========================================================================

fprintf('--- STEP 5: Setting calibrated parameters (M04) ---\n');

try
    calibration_results = M04_set_R_parameters();
    params = calibration_results.params;
    fprintf('  Parameters loaded from M04\n');
catch ME
    fprintf('  M04 returned error: %s\n', ME.message);
    fprintf('  Using default parameters from M01\n');
    params = default_params;
    calibration_results.params = params;
    calibration_results.param_info = param_info;
end

% Ensure calibration_results is saved for downstream scripts
save('model_output/calibration_results.mat', 'calibration_results');
fprintf('  Saved: model_output/calibration_results.mat\n\n');

%% =========================================================================
%% STEP 6: Verification simulation
%% =========================================================================

fprintf('--- STEP 6: Verification simulation ---\n');

sim_p = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
sim_t = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
te = sim_t.JSW_change(end) - sim_p.JSW_change(end);

fprintf('  3-Year Placebo JSW change:  %.4f mm\n', sim_p.JSW_change(end));
fprintf('  3-Year Treatment JSW change: %.4f mm\n', sim_t.JSW_change(end));
fprintf('  3-Year Treatment Effect:     %.4f mm\n\n', te);

%% =========================================================================
%% STEP 7: Sensitivity Analysis (M05)
%% =========================================================================

fprintf('--- STEP 7: Sensitivity Analysis (M05) ---\n');

try
    sensitivity_results = M05_sensitivity_analysis();
    fprintf('  Sensitivity analysis complete\n\n');
catch ME
    fprintf('  WARNING: Sensitivity analysis failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 8: Profile Likelihood (M06)
%% =========================================================================

fprintf('--- STEP 8: Profile Likelihood (M06) ---\n');

try
    profile_results = M06_profile_likelihood();
    fprintf('  Profile likelihood complete\n\n');
catch ME
    fprintf('  WARNING: Profile likelihood failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 9: Bootstrap Confidence Intervals (M07)
%% *** THIS IS THE SLOW STEP - may take 2-8 hours ***
%% =========================================================================

fprintf('--- STEP 9: Bootstrap CI (M07) - this may take HOURS ---\n');
fprintf('  Started: %s\n', string(datetime('now')));

try
    bootstrap_results = M07_bootstrap_ci();
    fprintf('  Bootstrap complete at: %s\n', string(datetime('now')));
    fprintf('  Valid replicates: %d / %d\n\n', ...
            bootstrap_results.n_valid, bootstrap_results.n_boot);
catch ME
    fprintf('  WARNING: Bootstrap failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 10: Residual Diagnostics (M08)
%% =========================================================================

fprintf('--- STEP 10: Residual Diagnostics (M08) ---\n');

try
    diagnostic_results = M08_residual_diagnostics();
    fprintf('  Residual diagnostics complete\n\n');
catch ME
    fprintf('  WARNING: Residual diagnostics failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 11: Virtual Population (M09)
%% =========================================================================

fprintf('--- STEP 11: Virtual Population (M09) ---\n');

try
    vpop_results = M09_virtual_population();
    fprintf('  Virtual population complete\n\n');
catch ME
    fprintf('  WARNING: Virtual population failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 12: External Validation (M10)
%% =========================================================================

fprintf('--- STEP 12: External Validation (M10) ---\n');

try
    validation_results = M10_external_validation();
    fprintf('  External validation complete\n\n');
catch ME
    fprintf('  WARNING: External validation failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 13: Translational Simulations (M11)
%% =========================================================================

fprintf('--- STEP 13: Translational Simulations (M11) ---\n');

try
    trans_results = M11_translational_simulations();
    fprintf('  Translational simulations complete\n\n');
catch ME
    fprintf('  WARNING: Translational simulations failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 14: Generate Main Figures (M12)
%% =========================================================================

fprintf('--- STEP 14: Generating Main Figures (M12) ---\n');

try
    run(fullfile('src', 'M12_figures.m'));
    fprintf('  All 5 main figures generated\n\n');
catch ME
    fprintf('  WARNING: Figure generation failed: %s\n\n', ME.message);
end

%% =========================================================================
%% STEP 15: Generate Supplementary Figures (M13)
%% =========================================================================

fprintf('--- STEP 15: Generating Supplementary Figures (M13) ---\n');

try
    run(fullfile('src', 'M13_supplementary_figures.m'));
    fprintf('  All 5 supplementary figures generated\n\n');
catch ME
    fprintf('  WARNING: Supplementary figure generation failed: %s\n\n', ME.message);
end

%% =========================================================================
%% FINAL SUMMARY
%% =========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('   PIPELINE COMPLETE\n');
fprintf('=============================================================\n');
fprintf('  Finished: %s\n\n', string(datetime('now')));

% Calibrated Parameters
fprintf('CALIBRATED PARAMETERS:\n');
fprintf('  k_syn     = %.6f /year\n', params.k_syn);
fprintf('  k_deg     = %.6f /year\n', params.k_deg);
fprintf('  Imax_deg  = %.4f\n', params.Imax_deg);
fprintf('  Pain_0    = %.2f WOMAC points\n', params.Pain_0);
fprintf('  Pmax      = %.2f WOMAC points\n\n', params.Pmax);

% Treatment Effect
sim_p_final = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
sim_t_final = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
te_final = sim_t_final.JSW_change(end) - sim_p_final.JSW_change(end);

fprintf('TREATMENT EFFECTS:\n');
fprintf('  3-Year JSW Treatment Effect: %.3f mm\n', te_final);
fprintf('  3-Year Placebo JSW change:   %.3f mm\n', sim_p_final.JSW_change(end));
fprintf('  3-Year Treatment JSW change: %.3f mm\n\n', sim_t_final.JSW_change(end));

% Pain
sim_p_pain = simulate_glucosamine(params, 1500, 24, 'HCl', false);
sim_t_pain = simulate_glucosamine(params, 1500, 24, 'HCl', true);

fprintf('PAIN RESULTS (24-week, HCl):\n');
fprintf('  Placebo Week 24 Pain:  %.1f\n', sim_p_pain.Pain(end));
fprintf('  Treatment Week 24 Pain: %.1f\n', sim_t_pain.Pain(end));
fprintf('  Pain difference:         %.1f WOMAC points\n\n', ...
        sim_p_pain.Pain(end) - sim_t_pain.Pain(end));

fprintf('=============================================================\n');
fprintf('Output locations:\n');
fprintf('  data_clean/    - Processed calibration data\n');
fprintf('  model_output/  - Saved model results (.mat files)\n');
fprintf('  figures/        - Publication figures (.png, .fig)\n');
fprintf('=============================================================\n');
