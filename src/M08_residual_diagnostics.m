%% =============================================================================
%% M08_residual_diagnostics.m
%% Residual Diagnostics and Goodness-of-Fit Assessment
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB SimBiology Implementation
%%
%% Outputs:
%%   - Weighted residuals for all calibration points
%%   - Goodness-of-fit metrics (RMSE, MAE, bias)
%%   - Bias test (t-test on weighted residuals)
%%
%% Output: model_output/residual_diagnostics.mat
%%
%% =============================================================================

function diagnostic_results = M08_residual_diagnostics()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('RESIDUAL DIAGNOSTICS\n');
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
    
    fprintf('  JSW observations: %d\n', height(jsw_targets));
    fprintf('  WOMAC observations: %d\n', height(womac_targets));
    
    %% =============================================================================
    %% 2. GENERATE PREDICTIONS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('CALCULATING RESIDUALS\n');
    fprintf('=============================================================\n\n');
    
    sim_placebo_3y = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', false);
    sim_treat_3y = simulate_glucosamine(calibrated_params, 1500, 156, 'sulfate', true);
    sim_placebo_24w = simulate_glucosamine(calibrated_params, 1500, 24, 'HCl', false);
    sim_treat_24w = simulate_glucosamine(calibrated_params, 1500, 24, 'HCl', true);
    
    % Initialize residual tables
    n_jsw = height(jsw_targets);
    n_womac = height(womac_targets);
    
    all_residuals = table('Size', [n_jsw + n_womac, 10], ...
        'VariableTypes', {'string', 'double', 'double', 'string', 'double', 'double', 'double', 'double', 'double', 'string'}, ...
        'VariableNames', {'study', 'week', 'time_years', 'arm', 'observed', 'predicted', 'residual', 'se_observed', 'weighted_residual', 'endpoint'});
    
    % JSW residuals
    for i = 1:n_jsw
        wk = jsw_targets.week(i);
        arm = jsw_targets.arm{i};
        
        if contains(arm, 'Placebo')
            pred = sim_placebo_3y.JSW_change(sim_placebo_3y.week == wk);
        else
            pred = sim_treat_3y.JSW_change(sim_treat_3y.week == wk);
        end
        
        obs = jsw_targets.jsw_change_mm(i);
        se = 0.08;
        
        all_residuals.study(i) = jsw_targets.study{i};
        all_residuals.week(i) = wk;
        all_residuals.time_years(i) = wk / 52;
        all_residuals.arm(i) = arm;
        all_residuals.observed(i) = obs;
        all_residuals.predicted(i) = pred;
        all_residuals.residual(i) = obs - pred;
        all_residuals.se_observed(i) = se;
        all_residuals.weighted_residual(i) = (obs - pred) / se;
        all_residuals.endpoint(i) = 'JSW';
    end
    
    % WOMAC residuals
    for i = 1:n_womac
        wk = womac_targets.week(i);
        arm = womac_targets.arm{i};
        
        if contains(arm, 'Placebo')
            pred = sim_placebo_24w.Pain(sim_placebo_24w.week == wk);
        else
            pred = sim_treat_24w.Pain(sim_treat_24w.week == wk);
        end
        
        obs = womac_targets.womac_pain_normalized(i);
        se = 5;
        
        idx = n_jsw + i;
        all_residuals.study(idx) = womac_targets.study{i};
        all_residuals.week(idx) = wk;
        all_residuals.time_years(idx) = wk / 52;
        all_residuals.arm(idx) = arm;
        all_residuals.observed(idx) = obs;
        all_residuals.predicted(idx) = pred;
        all_residuals.residual(idx) = obs - pred;
        all_residuals.se_observed(idx) = se;
        all_residuals.weighted_residual(idx) = (obs - pred) / se;
        all_residuals.endpoint(idx) = 'WOMAC';
    end
    
    % Summary by endpoint
    fprintf('Residual Summary by Endpoint:\n\n');
    
    endpoints = unique(all_residuals.endpoint);
    resid_summary = table('Size', [length(endpoints), 5], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'endpoint', 'n', 'mean_resid', 'sd_resid', 'mean_wresid'});
    
    for e = 1:length(endpoints)
        ep = endpoints{e};
        ep_data = all_residuals(strcmp(all_residuals.endpoint, ep), :);
        
        resid_summary.endpoint(e) = ep;
        resid_summary.n(e) = height(ep_data);
        resid_summary.mean_resid(e) = mean(ep_data.residual);
        resid_summary.sd_resid(e) = std(ep_data.residual);
        resid_summary.mean_wresid(e) = mean(ep_data.weighted_residual);
    end
    
    disp(resid_summary);
    
    %% =============================================================================
    %% 3. GOODNESS-OF-FIT METRICS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('GOODNESS-OF-FIT METRICS\n');
    fprintf('=============================================================\n\n');
    
    gof_metrics = table('Size', [length(endpoints), 7], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'endpoint', 'n', 'RMSE', 'MAE', 'Mean_Bias', 'Max_Abs_Resid', 'Mean_Weighted_Resid'});
    
    for e = 1:length(endpoints)
        ep = endpoints{e};
        ep_data = all_residuals(strcmp(all_residuals.endpoint, ep), :);
        
        gof_metrics.endpoint(e) = ep;
        gof_metrics.n(e) = height(ep_data);
        gof_metrics.RMSE(e) = sqrt(mean(ep_data.residual.^2));
        gof_metrics.MAE(e) = mean(abs(ep_data.residual));
        gof_metrics.Mean_Bias(e) = mean(ep_data.residual);
        gof_metrics.Max_Abs_Resid(e) = max(abs(ep_data.residual));
        gof_metrics.Mean_Weighted_Resid(e) = mean(ep_data.weighted_residual);
    end
    
    fprintf('By Endpoint:\n');
    disp(gof_metrics);
    
    % Overall weighted residuals
    overall_wresid = all_residuals.weighted_residual;
    
    fprintf('\nOverall Weighted Residuals:\n');
    mean_wresid = mean(overall_wresid);
    se_mean = std(overall_wresid) / sqrt(length(overall_wresid));
    fprintf('  Mean: %.3f (95%% CI: %.3f to %.3f)\n', ...
            mean_wresid, mean_wresid - 1.96*se_mean, mean_wresid + 1.96*se_mean);
    fprintf('  SD: %.3f\n', std(overall_wresid));
    
    % One-sample t-test for bias
    [~, p_value, ~, stats] = ttest(overall_wresid);
    
    if p_value > 0.05
        bias_status = 'No significant bias';
    else
        bias_status = 'SIGNIFICANT BIAS';
    end
    fprintf('  Bias test p-value: %.4f (%s)\n', p_value, bias_status);
    
    %% =============================================================================
    %% 4. OBSERVED VS PREDICTED TABLE
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('OBSERVED VS PREDICTED\n');
    fprintf('=============================================================\n\n');
    
    obs_vs_pred_table = all_residuals(:, {'endpoint', 'study', 'week', 'arm', 'observed', 'predicted', 'residual'});
    obs_vs_pred_table.Properties.VariableNames = {'Endpoint', 'Study', 'Week', 'Arm', 'Observed', 'Predicted', 'Residual'};
    
    disp(obs_vs_pred_table);
    
    %% =============================================================================
    %% 5. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    diagnostic_results.residuals = all_residuals;
    diagnostic_results.gof_metrics = gof_metrics;
    diagnostic_results.bias_test.p_value = p_value;
    diagnostic_results.bias_test.t_stat = stats.tstat;
    diagnostic_results.bias_test.df = stats.df;
    diagnostic_results.obs_vs_pred = obs_vs_pred_table;
    
    save('model_output/residual_diagnostics.mat', 'diagnostic_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/residual_diagnostics.mat\n');
    fprintf('=============================================================\n');
    
end
