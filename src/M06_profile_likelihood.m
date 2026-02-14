%% =============================================================================
%% M06_profile_likelihood.m
%% Profile Likelihood Analysis for Parameter Identifiability (FIXED)
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB Implementation
%%
%% Method: For each parameter, fix at range of values and compute objective
%%         95% CI from chi-squared threshold (delta = 1.92 for df=1)
%%
%% =============================================================================

function profile_results = M06_profile_likelihood()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('PROFILE LIKELIHOOD ANALYSIS\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD MODEL AND DATA
    %% =============================================================================
    
    fprintf('Loading calibrated model and data...\n');
    
    load('model_output/calibration_results.mat', 'calibration_results');
    calibrated_params = calibration_results.params;
    
    load('data_clean/calibration_data.mat', 'calibration_data');
    jsw_data = calibration_data.jsw_data;
    womac_data = calibration_data.womac_pain_data;
    
    jsw_targets = jsw_data(jsw_data.week > 0, :);
    womac_targets = womac_data(strcmp(womac_data.study, 'GAIT 2006'), :);
    
    % Compute baseline objective with calibrated parameters
    baseline_obj = compute_objective(calibrated_params, jsw_targets, womac_targets);
    fprintf('  Baseline objective: %.4f\n', baseline_obj);
    
    %% =============================================================================
    %% 2. SETUP
    %% =============================================================================
    
    % Parameters to profile
    profile_params = {'k_syn', 'k_deg', 'Imax_deg', 'Pain_0', 'Pmax'};
    
    % Define parameter ranges for profiling
    param_ranges.k_syn = linspace(0.08, 0.20, 15);
    param_ranges.k_deg = linspace(0.12, 0.25, 15);
    param_ranges.Imax_deg = linspace(0.4, 1.0, 15);
    param_ranges.Pain_0 = linspace(35, 55, 15);
    param_ranges.Pmax = linspace(10, 28, 15);
    
    n_points = 15;
    chi2_threshold = 1.92;  % 95% CI threshold for 1 df
    
    profile_data = struct();
    
    %% =============================================================================
    %% 3. PROFILE LIKELIHOOD COMPUTATION (SIMPLIFIED)
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('COMPUTING PROFILE LIKELIHOODS\n');
    fprintf('=============================================================\n\n');
    
    for p = 1:length(profile_params)
        prof_param = profile_params{p};
        fprintf('Profiling %s...', prof_param);
        
        profile_values = param_ranges.(prof_param);
        profile_obj = zeros(length(profile_values), 1);
        
        for j = 1:length(profile_values)
            % Create parameter set with fixed profiled parameter
            test_params = calibrated_params;
            test_params.(prof_param) = profile_values(j);
            
            % Compute objective
            profile_obj(j) = compute_objective(test_params, jsw_targets, womac_targets);
        end
        
        profile_data.(prof_param).param_value = profile_values';
        profile_data.(prof_param).objective = profile_obj;
        profile_data.(prof_param).delta_obj = profile_obj - min(profile_obj);
        
        fprintf(' Done (min obj = %.2f)\n', min(profile_obj));
    end
    
    %% =============================================================================
    %% 4. CALCULATE CONFIDENCE INTERVALS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('95%% CONFIDENCE INTERVALS\n');
    fprintf('=============================================================\n\n');
    
    ci_summary = table('Size', [length(profile_params), 5], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'parameter', 'estimate', 'ci_lower', 'ci_upper', 'identifiable'});
    
    for p = 1:length(profile_params)
        prof_param = profile_params{p};
        prof_data = profile_data.(prof_param);
        cal_val = calibrated_params.(prof_param);
        
        % Find where profile crosses threshold
        below_threshold = prof_data.delta_obj <= chi2_threshold;
        
        if sum(below_threshold) > 0
            valid_vals = prof_data.param_value(below_threshold);
            ci_lower = min(valid_vals);
            ci_upper = max(valid_vals);
            identifiable = true;
        else
            ci_lower = min(prof_data.param_value);
            ci_upper = max(prof_data.param_value);
            identifiable = false;
        end
        
        ci_summary.parameter(p) = prof_param;
        ci_summary.estimate(p) = cal_val;
        ci_summary.ci_lower(p) = ci_lower;
        ci_summary.ci_upper(p) = ci_upper;
        ci_summary.identifiable(p) = identifiable;
        
        if identifiable
            fprintf('  %s: %.4f [%.4f, %.4f]\n', prof_param, cal_val, ci_lower, ci_upper);
        else
            fprintf('  %s: %.4f [NOT IDENTIFIABLE]\n', prof_param, cal_val);
        end
    end
    
    %% =============================================================================
    %% 5. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    profile_results.profiles = profile_data;
    profile_results.ci_summary = ci_summary;
    profile_results.baseline_objective = baseline_obj;
    profile_results.chi2_threshold = chi2_threshold;
    profile_results.n_profile_points = n_points;
    
    save('model_output/profile_likelihood_results.mat', 'profile_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/profile_likelihood_results.mat\n');
    fprintf('=============================================================\n');
    
end

%% =============================================================================
%% HELPER FUNCTION: COMPUTE OBJECTIVE
%% =============================================================================

function sse = compute_objective(params, jsw_targets, womac_targets)
    % Compute sum of squared errors for given parameters
    
    total_sse = 0;
    
    % JSW component
    try
        sim_p = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
        sim_t = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
        
        for i = 1:height(jsw_targets)
            wk = jsw_targets.week(i);
            obs = jsw_targets.jsw_change_mm(i);
            arm = jsw_targets.arm{i};
            
            if contains(arm, 'Placebo')
                idx = find(sim_p.week == wk);
                if ~isempty(idx)
                    pred = sim_p.JSW_change(idx);
                else
                    continue;
                end
            else
                idx = find(sim_t.week == wk);
                if ~isempty(idx)
                    pred = sim_t.JSW_change(idx);
                else
                    continue;
                end
            end
            
            % Weighted residual (SE ~ 0.08 mm)
            total_sse = total_sse + ((pred - obs) / 0.08)^2;
        end
    catch
        total_sse = total_sse + 100;
    end
    
    % WOMAC component
    try
        sim_p24 = simulate_glucosamine(params, 1500, 24, 'HCl', false);
        sim_t24 = simulate_glucosamine(params, 1500, 24, 'HCl', true);
        
        for i = 1:height(womac_targets)
            wk = womac_targets.week(i);
            obs = womac_targets.womac_pain_normalized(i);
            arm = womac_targets.arm{i};
            
            if contains(arm, 'Placebo')
                idx = find(sim_p24.week == wk);
                if ~isempty(idx)
                    pred = sim_p24.Pain(idx);
                else
                    continue;
                end
            else
                idx = find(sim_t24.week == wk);
                if ~isempty(idx)
                    pred = sim_t24.Pain(idx);
                else
                    continue;
                end
            end
            
            % Lower weight for pain (SE ~ 10 points)
            total_sse = total_sse + 0.1 * ((pred - obs) / 10)^2;
        end
    catch
        total_sse = total_sse + 50;
    end
    
    sse = total_sse;
end
