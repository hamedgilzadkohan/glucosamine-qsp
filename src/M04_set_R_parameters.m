%% =============================================================================
%% M04_set_R_parameters.m
%% Set Calibrated Parameters to Match R Manuscript Exactly
%% =============================================================================
%%
%% This script sets the calibrated parameters to EXACTLY match the R model
%% so that MATLAB-generated figures are consistent with manuscript text.
%%
%% R Calibrated Values (from calibrated_parameters_v3.xlsx):
%%   k_syn    = 0.175258 /year
%%   k_deg    = 0.207263 /year
%%   Imax_deg = 1.000000 (at boundary)
%%   Pain_0   = 46.651892
%%   Pmax     = 17.624509
%%
%% =============================================================================

function calibration_results = M04_set_R_parameters()

    fprintf('=============================================================\n');
    fprintf('SETTING PARAMETERS TO MATCH R MANUSCRIPT\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD MODEL STRUCTURE
    %% =============================================================================
    
    [param_info, base_params, ~] = M01_model_structure();
    
    %% =============================================================================
    %% 2. SET R CALIBRATED PARAMETERS
    %% =============================================================================
    
    % Start with base parameters
    calibrated_params = base_params;
    
    % Overwrite with R calibrated values
    % Optimized parameters from R
    calibrated_params.k_syn    = 0.175258;    % GAG synthesis rate (/year)
    calibrated_params.k_deg    = 0.207263;    % GAG degradation rate (/year)
    calibrated_params.Imax_deg = 1.000000;    % Max degradation inhibition
    calibrated_params.Pain_0   = 46.651892;   % Baseline pain (0-100)
    calibrated_params.Pmax     = 17.624509;   % Maximum placebo effect
    
    % Fixed parameters from R
    calibrated_params.F_sulfate   = 0.22;     % Bioavailability (crystalline sulfate)
    calibrated_params.F_HCl       = 0.11;     % Bioavailability (HCl formulation)
    calibrated_params.Vd_apparent = 200;      % Apparent volume of distribution (L)
    calibrated_params.R_syn       = 0.25;     % Synovial:plasma concentration ratio
    calibrated_params.GAG_0       = 0.95;     % Initial GAG content (normalized)
    calibrated_params.JSW_0       = 4.2;      % Baseline joint space width (mm)
    calibrated_params.gamma_JSW   = 1.0;      % GAG-JSW power coefficient
    calibrated_params.IC50_deg    = 3.0;      % IC50 for degradation inhibition (µug/mL)
    calibrated_params.alpha_struct = 150;     % Pain sensitivity to GAG loss
    calibrated_params.beta_drug   = 5.0;      % Direct drug effect on pain
    calibrated_params.EC50_pain   = 3.0;      % EC50 for pain relief (µug/mL)
    calibrated_params.kpl         = 0.02;     % Placebo effect rate constant (/day)
    
    %% =============================================================================
    %% 3. DISPLAY PARAMETERS
    %% =============================================================================
    
    fprintf('CALIBRATED PARAMETERS (matched to R):\n');
    fprintf('  k_syn     = %.6f /year\n', calibrated_params.k_syn);
    fprintf('  k_deg     = %.6f /year\n', calibrated_params.k_deg);
    fprintf('  Imax_deg  = %.6f\n', calibrated_params.Imax_deg);
    fprintf('  Pain_0    = %.6f\n', calibrated_params.Pain_0);
    fprintf('  Pmax      = %.6f\n\n', calibrated_params.Pmax);
    
    fprintf('FIXED PARAMETERS:\n');
    fprintf('  F_sulfate   = %.2f\n', calibrated_params.F_sulfate);
    fprintf('  F_HCl       = %.2f\n', calibrated_params.F_HCl);
    fprintf('  GAG_0       = %.2f\n', calibrated_params.GAG_0);
    fprintf('  JSW_0       = %.1f mm\n', calibrated_params.JSW_0);
    fprintf('  IC50_deg    = %.1f µug/mL\n', calibrated_params.IC50_deg);
    fprintf('  kpl         = %.3f /day\n', calibrated_params.kpl);
    
    %% =============================================================================
    %% 4. VERIFY MODEL PREDICTIONS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('VERIFYING MODEL PREDICTIONS\n');
    fprintf('=============================================================\n\n');
    
    % Simulate with R parameters
    sim_p = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', false);
    sim_t = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', true);
    
    % 3-year endpoints (week 156)
    y3_placebo = sim_p.JSW_change(157);
    y3_treat = sim_t.JSW_change(157);
    treatment_effect = y3_treat - y3_placebo;
    
    fprintf('JSW PREDICTIONS (3-year):\n');
    fprintf('  Placebo change:    %.4f mm\n', y3_placebo);
    fprintf('  Treatment change:  %.4f mm\n', y3_treat);
    fprintf('  Treatment effect:  %.4f mm\n', treatment_effect);
    fprintf('  R manuscript:      ~0.22 mm\n');
    
    % Pain predictions
    w24_pain_placebo = sim_p.Pain(sim_p.week == 24);
    w24_pain_treat = sim_t.Pain(sim_t.week == 24);
    
    fprintf('\nPAIN PREDICTIONS (Week 24):\n');
    fprintf('  Placebo:    %.1f points\n', w24_pain_placebo);
    fprintf('  Treatment:  %.1f points\n', w24_pain_treat);
    
    %% =============================================================================
    %% 5. CREATE RESULTS STRUCTURE (compatible with downstream scripts)
    %% =============================================================================
    
    calibration_results = struct();
    calibration_results.params = calibrated_params;
    calibration_results.param_info = param_info;
    calibration_results.objective_value = 0;  % Not optimized, set to R manuscript
    calibration_results.par_names = {'k_syn', 'k_deg', 'Imax_deg', 'Pain_0', 'Pmax'};
    calibration_results.final_params = [calibrated_params.k_syn, ...
                                         calibrated_params.k_deg, ...
                                         calibrated_params.Imax_deg, ...
                                         calibrated_params.Pain_0, ...
                                         calibrated_params.Pmax];
    calibration_results.optimization_results = {};
    calibration_results.method = 'Fixed to R manuscript values';
    
    % Store predictions
    calibration_results.predictions.jsw_y3_placebo = y3_placebo;
    calibration_results.predictions.jsw_y3_treat = y3_treat;
    calibration_results.predictions.treatment_effect = treatment_effect;
    calibration_results.predictions.pain_w24_placebo = w24_pain_placebo;
    calibration_results.predictions.pain_w24_treat = w24_pain_treat;
    
    %% =============================================================================
    %% 6. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    % Save in format expected by downstream scripts
    save('model_output/calibrated_params.mat', 'calibrated_params');
    save('model_output/calibration_results.mat', 'calibration_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to:\n');
    fprintf('  model_output/calibrated_params.mat\n');
    fprintf('  model_output/calibration_results.mat\n');
    fprintf('=============================================================\n');
    fprintf('\nPARAMETERS SET TO R MANUSCRIPT VALUES - READY FOR ANALYSIS\n');
    fprintf('=============================================================\n');
    
end
