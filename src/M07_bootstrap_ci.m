%% =============================================================================
%% M07_bootstrap_ci.m
%% Bootstrap Confidence Intervals for Model Parameters
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB SimBiology Implementation
%%
%% Method: Non-parametric bootstrap with stratified resampling by arm
%%         Re-optimize for each bootstrap sample
%%
%% Output: model_output/bootstrap_results.mat
%%
%% =============================================================================

function bootstrap_results = M07_bootstrap_ci()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('BOOTSTRAP CONFIDENCE INTERVALS\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD DATA AND MODEL
    %% =============================================================================
    
    fprintf('Loading calibrated model and data...\n');
    
    load('model_output/calibration_results.mat', 'calibration_results');
    calibrated_params = calibration_results.params;
    param_info = calibration_results.param_info;
    
    load('data_clean/calibration_data.mat', 'calibration_data');
    jsw_data = calibration_data.jsw_data;
    womac_data = calibration_data.womac_pain_data;
    
    jsw_targets = jsw_data(jsw_data.week > 0, :);
    womac_targets = womac_data(strcmp(womac_data.study, 'GAIT 2006'), :);
    
    %% =============================================================================
    %% 2. SETUP
    %% =============================================================================
    
    opt_idx = param_info.Optimized;
    par_names = param_info.Name(opt_idx);
    par_lower = param_info.Lower(opt_idx)';
    par_upper = param_info.Upper(opt_idx)';
    
    n_params = length(par_names);
    par_init = zeros(n_params, 1);
    for i = 1:n_params
        par_init(i) = calibrated_params.(par_names{i});
    end
    
    %% =============================================================================
    %% 3. OBJECTIVE FUNCTION
    %% =============================================================================
    
    function sse = objective_boot(par_vector, jsw_boot, womac_boot)
        params = calibrated_params;
        for i = 1:n_params
            params.(par_names{i}) = par_vector(i);
        end
        
        if any(par_vector <= 0), sse = 1e12; return; end
        if params.Imax_deg > 1, sse = 1e12; return; end
        if params.k_deg <= params.k_syn, sse = 1e12; return; end
        
        total_sse = 0;
        
        try
            sim_p = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
            sim_t = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
            
            for i = 1:height(jsw_boot)
                wk = jsw_boot.week(i);
                obs = jsw_boot.jsw_change_mm(i);
                arm = jsw_boot.arm{i};
                
                if contains(arm, 'Placebo')
                    idx = find(sim_p.week == wk);
                    if ~isempty(idx), pred = sim_p.JSW_change(idx); else, continue; end
                else
                    idx = find(sim_t.week == wk);
                    if ~isempty(idx), pred = sim_t.JSW_change(idx); else, continue; end
                end
                
                total_sse = total_sse + ((pred - obs) / 0.08)^2;
            end
        catch
            total_sse = 1e10;
        end
        
        try
            sim_p24 = simulate_glucosamine(params, 1500, 24, 'HCl', false);
            sim_t24 = simulate_glucosamine(params, 1500, 24, 'HCl', true);
            
            for i = 1:height(womac_boot)
                wk = womac_boot.week(i);
                obs = womac_boot.womac_pain_normalized(i);
                arm = womac_boot.arm{i};
                
                if contains(arm, 'Placebo')
                    idx = find(sim_p24.week == wk);
                    if ~isempty(idx), pred = sim_p24.Pain(idx); else, continue; end
                else
                    idx = find(sim_t24.week == wk);
                    if ~isempty(idx), pred = sim_t24.Pain(idx); else, continue; end
                end
                
                total_sse = total_sse + 0.1 * ((pred - obs) / 10)^2;
            end
        catch
            total_sse = total_sse + 1e8;
        end
        
        sse = total_sse;
    end
    
    %% =============================================================================
    %% 4. BOOTSTRAP RESAMPLING
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('BOOTSTRAP RESAMPLING\n');
    fprintf('=============================================================\n\n');
    
    n_boot = 1000;
    fprintf('Running %d bootstrap replicates...\n\n', n_boot);
    
    boot_params = zeros(n_boot, n_params);
    boot_sse = zeros(n_boot, 1);
    boot_converged = false(n_boot, 1);
    
    options = optimoptions('fmincon', 'Display', 'off', 'MaxIterations', 200);
    
    pb_interval = ceil(n_boot / 20);
    tic;
    
    % Get unique arms for stratified sampling
    jsw_arms = unique(jsw_targets.arm);
    womac_arms = unique(womac_targets.arm);
    
    for b = 1:n_boot
        if mod(b, pb_interval) == 0
            fprintf('  %d/%d (%.0f%%)\n', b, n_boot, 100*b/n_boot);
        end
        
        % Stratified resampling of JSW data
        jsw_boot = table();
        for a = 1:length(jsw_arms)
            arm_data = jsw_targets(strcmp(jsw_targets.arm, jsw_arms{a}), :);
            n_arm = height(arm_data);
            boot_idx = randi(n_arm, n_arm, 1);
            jsw_boot = [jsw_boot; arm_data(boot_idx, :)];
        end
        
        % Stratified resampling of WOMAC data
        womac_boot = table();
        for a = 1:length(womac_arms)
            arm_data = womac_targets(strcmp(womac_targets.arm, womac_arms{a}), :);
            n_arm = height(arm_data);
            boot_idx = randi(n_arm, n_arm, 1);
            womac_boot = [womac_boot; arm_data(boot_idx, :)];
        end
        
        % Re-calibrate on bootstrap sample
        try
            [x_opt, fval, exitflag] = fmincon(@(x) objective_boot(x, jsw_boot, womac_boot), ...
                                               par_init, [], [], [], [], ...
                                               par_lower, par_upper, [], options);
            
            boot_params(b, :) = x_opt';
            boot_sse(b) = fval;
            boot_converged(b) = (exitflag > 0);
            
        catch
            boot_params(b, :) = NaN;
            boot_sse(b) = NaN;
            boot_converged(b) = false;
        end
    end
    
    elapsed_time = toc;
    fprintf('\nBootstrap completed in %.1f minutes\n', elapsed_time/60);
    
    %% =============================================================================
    %% 5. COMPUTE CONFIDENCE INTERVALS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('BOOTSTRAP RESULTS\n');
    fprintf('=============================================================\n\n');
    
    valid_idx = ~isnan(boot_params(:, 1)) & boot_converged;
    n_valid = sum(valid_idx);
    fprintf('Valid samples: %d/%d (%.1f%%)\n\n', n_valid, n_boot, 100*n_valid/n_boot);
    
    boot_params_valid = boot_params(valid_idx, :);
    
    boot_summary = table('Size', [n_params, 7], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Parameter', 'Estimate', 'Boot_Mean', 'Boot_SD', 'CI_Lower', 'CI_Upper', 'CV_percent'});
    
    for i = 1:n_params
        boot_summary.Parameter(i) = par_names{i};
        boot_summary.Estimate(i) = calibrated_params.(par_names{i});
        boot_summary.Boot_Mean(i) = mean(boot_params_valid(:, i));
        boot_summary.Boot_SD(i) = std(boot_params_valid(:, i));
        boot_summary.CI_Lower(i) = prctile(boot_params_valid(:, i), 2.5);
        boot_summary.CI_Upper(i) = prctile(boot_params_valid(:, i), 97.5);
        boot_summary.CV_percent(i) = 100 * std(boot_params_valid(:, i)) / mean(boot_params_valid(:, i));
    end
    
    fprintf('Parameter Estimates with 95%% Bootstrap CIs:\n\n');
    fprintf('  %-12s  %10s  %20s  %8s\n', 'Parameter', 'Estimate', '95% CI', 'CV%');
    fprintf('%s\n', repmat('-', 1, 55));
    
    for i = 1:n_params
        ci_str = sprintf('[%.4f, %.4f]', boot_summary.CI_Lower(i), boot_summary.CI_Upper(i));
        fprintf('  %-12s  %10.4f  %20s  %8.1f\n', ...
                boot_summary.Parameter{i}, ...
                boot_summary.Estimate(i), ...
                ci_str, ...
                boot_summary.CV_percent(i));
    end
    
    %% =============================================================================
    %% 6. PARAMETER CORRELATION
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('PARAMETER CORRELATION MATRIX\n');
    fprintf('=============================================================\n\n');
    
    cor_matrix = corrcoef(boot_params_valid);
    
    % Display correlation matrix
    fprintf('         ');
    for i = 1:n_params
        fprintf('%8s ', par_names{i}(1:min(8, length(par_names{i}))));
    end
    fprintf('\n');
    
    for i = 1:n_params
        fprintf('%-8s ', par_names{i}(1:min(8, length(par_names{i}))));
        for j = 1:n_params
            fprintf('%8.3f ', cor_matrix(i, j));
        end
        fprintf('\n');
    end
    
    % Find high correlations
    fprintf('\nHigh correlations (|r| > 0.8):\n');
    found_high = false;
    for i = 1:n_params
        for j = (i+1):n_params
            if abs(cor_matrix(i, j)) > 0.8
                fprintf('  %s - %s: r = %.3f\n', par_names{i}, par_names{j}, cor_matrix(i, j));
                found_high = true;
            end
        end
    end
    if ~found_high
        fprintf('  No parameter pairs with |correlation| > 0.8\n');
    end
    
    %% =============================================================================
    %% 7. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    bootstrap_results.boot_params = boot_params_valid;
    bootstrap_results.boot_summary = boot_summary;
    bootstrap_results.correlation_matrix = cor_matrix;
    bootstrap_results.n_boot = n_boot;
    bootstrap_results.n_valid = n_valid;
    bootstrap_results.runtime_minutes = elapsed_time / 60;
    bootstrap_results.par_names = par_names;
    
    save('model_output/bootstrap_results.mat', 'bootstrap_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/bootstrap_results.mat\n');
    fprintf('=============================================================\n');
    
end
