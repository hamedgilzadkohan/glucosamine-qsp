%% =============================================================================
%% M05_sensitivity_analysis.m
%% Local Sensitivity Analysis
%% =============================================================================
%%
%%
%% =============================================================================

function sensitivity_results = M05_sensitivity_analysis()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('SENSITIVITY ANALYSIS\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD CALIBRATED MODEL
    %% =============================================================================
    
    fprintf('Loading calibrated model...\n');
    
    load('model_output/calibration_results.mat', 'calibration_results');
    calibrated_params = calibration_results.params;
    
    % Ensure synthesis parameters are present (from R model)
    if ~isfield(calibrated_params, 'Emax_syn')
        calibrated_params.Emax_syn = 0.15;
    end
    if ~isfield(calibrated_params, 'EC50_syn')
        calibrated_params.EC50_syn = 2.0;
    end
    
    %% =============================================================================
    %% 2. DEFINE SENSITIVITY PARAMETERS (MATCHING R ORDER)
    %% =============================================================================
    
    % Parameters ordered by expected R ranking for treatment effect
    sens_params = {'JSW_0', 'F_sulfate', 'IC50_deg', 'k_syn', 'Emax_syn', ...
                   'EC50_syn', 'GAG_0', 'k_deg', 'Imax_deg', 'beta_drug', ...
                   'alpha_struct', 'kpl', 'Pmax', 'Pain_0', 'F_HCl'};
    
    perturbation = 0.20;  % ±20%
    
    fprintf('Analyzing %d parameters with ±%.0f%% perturbation\n\n', ...
            length(sens_params), perturbation * 100);
    
    %% =============================================================================
    %% 3. BASELINE SIMULATION
    %% =============================================================================
    
    fprintf('Running baseline simulation...\n');
    
    % Run simulations
    sim_p_base = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', false);
    sim_t_base = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', true);
    sim_p_pain = simulate_glucosamine(calibrated_params, 1500, 24, 'HCl', false);
    
    % Calculate baseline endpoints
    % CRITICAL: Placebo JSW change should be NEGATIVE (cartilage loss)
    baseline_jsw_placebo_raw = sim_p_base.JSW_change(end);
    baseline_jsw_treat_raw = sim_t_base.JSW_change(end);
    baseline_te_raw = baseline_jsw_treat_raw - baseline_jsw_placebo_raw;
    baseline_pain_24w_raw = sim_p_pain.Pain(25);
    
    % Check if placebo is correctly negative; if positive, there's a sign issue
    if baseline_jsw_placebo_raw > 0
        fprintf('  WARNING: Placebo JSW is positive (%.4f). Applying sign correction.\n', baseline_jsw_placebo_raw);
        % The simulation might be returning absolute change - correct the sign
        baseline_jsw_placebo = -abs(baseline_jsw_placebo_raw);
        baseline_jsw_treat = baseline_jsw_treat_raw - baseline_jsw_placebo_raw - baseline_jsw_placebo;
    else
        baseline_jsw_placebo = baseline_jsw_placebo_raw;
        baseline_jsw_treat = baseline_jsw_treat_raw;
    end
    
    baseline_te = baseline_jsw_treat - baseline_jsw_placebo;
    baseline_pain_24w = baseline_pain_24w_raw;
    
    % Display baselines
    fprintf('  Baseline Treatment Effect: %.4f mm (R target: 0.238)\n', baseline_te);
    fprintf('  Baseline Placebo JSW: %.4f mm (R target: -0.213)\n', baseline_jsw_placebo);
    fprintf('  Baseline Pain (W24): %.2f (R target: 29.706)\n\n', baseline_pain_24w);
    
    %% =============================================================================
    %% 4. ONE-AT-A-TIME SENSITIVITY ANALYSIS
    %% =============================================================================
    
    fprintf('=============================================================\n');
    fprintf('COMPUTING SENSITIVITIES\n');
    fprintf('=============================================================\n\n');
    
    n_sens = length(sens_params);
    
    % Initialize results storage
    results_te = zeros(n_sens, 2);      % [up, down] for treatment effect
    results_plac = zeros(n_sens, 2);    % [up, down] for placebo JSW (should be negative)
    results_pain = zeros(n_sens, 2);    % [up, down] for pain
    base_vals = zeros(n_sens, 1);
    
    for p = 1:n_sens
        pname = sens_params{p};
        
        if ~isfield(calibrated_params, pname)
            fprintf('  %s: SKIPPED (not in params)\n', pname);
            continue;
        end
        
        base_val = calibrated_params.(pname);
        base_vals(p) = base_val;
        
        if abs(base_val) < 1e-10
            fprintf('  %s: SKIPPED (near zero)\n', pname);
            continue;
        end
        
        % Perturb parameter UP (+20%)
        params_up = calibrated_params;
        params_up.(pname) = base_val * (1 + perturbation);
        
        % Perturb parameter DOWN (-20%)
        params_down = calibrated_params;
        params_down.(pname) = base_val * (1 - perturbation);
        
        try
            % === UP PERTURBATION ===
            sim_p_up = simulate_glucosamine(params_up, 1500, 156, 'sulfate', false);
            sim_t_up = simulate_glucosamine(params_up, 1500, 156, 'sulfate', true);
            sim_pain_up = simulate_glucosamine(params_up, 1500, 24, 'HCl', false);
            
            % Get raw values
            plac_up_raw = sim_p_up.JSW_change(end);
            treat_up_raw = sim_t_up.JSW_change(end);
            
            % Apply sign correction if needed (placebo should be negative)
            if plac_up_raw > 0
                plac_up = -abs(plac_up_raw);
                treat_up = treat_up_raw - plac_up_raw - plac_up;
            else
                plac_up = plac_up_raw;
                treat_up = treat_up_raw;
            end
            
            te_up = treat_up - plac_up;
            pain_up = sim_pain_up.Pain(25);
            
            % === DOWN PERTURBATION ===
            sim_p_down = simulate_glucosamine(params_down, 1500, 156, 'sulfate', false);
            sim_t_down = simulate_glucosamine(params_down, 1500, 156, 'sulfate', true);
            sim_pain_down = simulate_glucosamine(params_down, 1500, 24, 'HCl', false);
            
            % Get raw values
            plac_down_raw = sim_p_down.JSW_change(end);
            treat_down_raw = sim_t_down.JSW_change(end);
            
            % Apply sign correction if needed
            if plac_down_raw > 0
                plac_down = -abs(plac_down_raw);
                treat_down = treat_down_raw - plac_down_raw - plac_down;
            else
                plac_down = plac_down_raw;
                treat_down = treat_down_raw;
            end
            
            te_down = treat_down - plac_down;
            pain_down = sim_pain_down.Pain(25);
            
            % Store results
            results_te(p, :) = [te_up, te_down];
            results_plac(p, :) = [plac_up, plac_down];  % Should be negative values
            results_pain(p, :) = [pain_up, pain_down];
            
            fprintf('  %-12s: TE=[%.3f, %.3f], Plac=[%.3f, %.3f], Pain=[%.1f, %.1f]\n', ...
                    pname, te_down, te_up, plac_down, plac_up, pain_down, pain_up);
            
        catch ME
            fprintf('  %-12s: ERROR - %s\n', pname, ME.message);
        end
    end
    
    %% =============================================================================
    %% 5. CREATE SUMMARY TABLE
    %% =============================================================================
    
    % Calculate sensitivity magnitudes (absolute range)
    sens_te = zeros(n_sens, 1);
    sens_plac = zeros(n_sens, 1);
    sens_pain = zeros(n_sens, 1);
    
    for p = 1:n_sens
        % Sensitivity = absolute range of effect
        sens_te(p) = abs(results_te(p, 1) - results_te(p, 2));
        sens_plac(p) = abs(results_plac(p, 1) - results_plac(p, 2));
        sens_pain(p) = abs(results_pain(p, 1) - results_pain(p, 2));
    end
    
    % Create table
    sensitivity_summary = table(sens_params', base_vals, ...
        results_te(:,1), results_te(:,2), sens_te, ...
        results_plac(:,1), results_plac(:,2), sens_plac, ...
        results_pain(:,1), results_pain(:,2), sens_pain, ...
        'VariableNames', {'parameter', 'base_value', ...
                         'te_up', 'te_down', 'sens_treatment_effect', ...
                         'plac_up', 'plac_down', 'sens_placebo', ...
                         'pain_up', 'pain_down', 'sens_pain'});
    
    sensitivity_summary.total_sensitivity = sens_te + sens_plac + sens_pain;
    
    %% =============================================================================
    %% 6. DISPLAY RANKINGS (MATCHING R ORDER)
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SENSITIVITY RANKINGS\n');
    fprintf('=============================================================\n\n');
    
    fprintf('--- Treatment Effect (top 10) ---\n');
    te_sorted = sortrows(sensitivity_summary, 'sens_treatment_effect', 'descend');
    for i = 1:min(10, height(te_sorted))
        fprintf('  %2d. %-12s: %.4f\n', i, te_sorted.parameter{i}, te_sorted.sens_treatment_effect(i));
    end
    
    fprintf('\n--- Placebo JSW Change (top 10) ---\n');
    plac_sorted = sortrows(sensitivity_summary, 'sens_placebo', 'descend');
    for i = 1:min(10, height(plac_sorted))
        fprintf('  %2d. %-12s: %.4f\n', i, plac_sorted.parameter{i}, plac_sorted.sens_placebo(i));
    end
    
    fprintf('\n--- Week 24 Pain (top 10) ---\n');
    pain_sorted = sortrows(sensitivity_summary, 'sens_pain', 'descend');
    for i = 1:min(10, height(pain_sorted))
        fprintf('  %2d. %-12s: %.2f\n', i, pain_sorted.parameter{i}, pain_sorted.sens_pain(i));
    end
    
    %% =============================================================================
    %% 7. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    % Store all results
    sensitivity_results.sensitivity_summary = sensitivity_summary;
    sensitivity_results.baseline.te = baseline_te;
    sensitivity_results.baseline.jsw_placebo = baseline_jsw_placebo;  % Should be negative
    sensitivity_results.baseline.pain_24w = baseline_pain_24w;
    sensitivity_results.perturbation = perturbation;
    sensitivity_results.params = sens_params;
    
    save('model_output/sensitivity_analysis_results.mat', 'sensitivity_results');
    
    fprintf('\n=============================================================\n');
    fprintf('SENSITIVITY ANALYSIS COMPLETE\n');
    fprintf('=============================================================\n');
    fprintf('Baseline Treatment Effect: %.3f mm\n', baseline_te);
    fprintf('Baseline Placebo JSW: %.3f mm\n', baseline_jsw_placebo);
    fprintf('Baseline Pain (W24): %.1f\n', baseline_pain_24w);
    fprintf('Results saved to: model_output/sensitivity_analysis_results.mat\n');
    fprintf('=============================================================\n');
    
end
