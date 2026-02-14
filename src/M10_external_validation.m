%% =============================================================================
%% M10_external_validation.m
%% External Validation Using GUIDE 2007 and MOVES 2015 Trials
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB SimBiology Implementation
%%
%% Validation trials:
%%   - GUIDE 2007: Glucosamine sulfate vs placebo (6 months)
%%   - MOVES 2015: CS+GH combination vs celecoxib (6 months)
%%
%% Output: model_output/validation_results.mat
%%
%% =============================================================================

function validation_results = M10_external_validation()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('EXTERNAL VALIDATION\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD MODEL AND IIV STRUCTURE
    %% =============================================================================
    
    fprintf('Loading calibrated model...\n');
    
    load('model_output/virtual_population_results.mat', 'vpop_results');
    base_params = vpop_results.calibrated_params;
    iiv_structure = vpop_results.iiv_structure;
    
    fprintf('  Parameters loaded\n');
    
    %% =============================================================================
    %% 2. VALIDATION TRIAL DATA
    %% =============================================================================
    
    fprintf('\n--- GUIDE 2007 Trial Data ---\n');
    
    guide_week = [0; 24; 0; 24];
    guide_arm = {'Placebo'; 'Placebo'; 'GS 1500mg'; 'GS 1500mg'};
    guide_pain_likert = [7.9; 6.1; 7.8; 5.1];
    guide_pain_likert_sd = [3.0; 3.2; 3.0; 3.0];
    guide_n = [104; 104; 106; 106];
    
    guide_pain_mean = guide_pain_likert * 5;  % Scale to 0-100
    guide_pain_sd = guide_pain_likert_sd * 5;
    guide_pain_se = guide_pain_sd ./ sqrt(guide_n);
    
    guide_observed = table(repmat({'GUIDE 2007'}, 4, 1), guide_week, guide_arm, ...
                           guide_pain_mean, guide_pain_sd, guide_pain_se, guide_n, ...
                           'VariableNames', {'study', 'week', 'arm', 'pain_mean', 'pain_sd', 'pain_se', 'n'});
    
    disp(guide_observed(:, {'arm', 'week', 'pain_mean', 'n'}));
    
    fprintf('\n--- MOVES 2015 Trial Data ---\n');
    
    moves_day = [0; 30; 60; 120; 180; 0; 30; 60; 120; 180];
    moves_week = [0; 4; 9; 17; 26; 0; 4; 9; 17; 26];
    moves_arm = {'CS+GH'; 'CS+GH'; 'CS+GH'; 'CS+GH'; 'CS+GH'; ...
                 'Celecoxib'; 'Celecoxib'; 'Celecoxib'; 'Celecoxib'; 'Celecoxib'};
    moves_pain_mean = [74.4; 53.5; 46.2; 42.0; 37.2; 74.1; 47.3; 41.2; 36.7; 36.9];
    moves_pain_se = [1.5; 1.8; 2.0; 2.1; 2.2; 1.5; 1.7; 1.9; 2.0; 2.1];
    moves_n = [264; 264; 264; 264; 264; 258; 258; 258; 258; 258];
    
    moves_observed = table(repmat({'MOVES 2015'}, 10, 1), moves_day, moves_week, moves_arm, ...
                           moves_pain_mean, moves_pain_se, moves_n, ...
                           'VariableNames', {'study', 'day', 'week', 'arm', 'pain_mean', 'pain_se', 'n'});
    
    disp(moves_observed(strcmp(moves_observed.arm, 'CS+GH'), {'arm', 'day', 'pain_mean', 'n'}));
    
    %% =============================================================================
    %% 3. SIMULATE GUIDE 2007 TRIAL
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATING GUIDE 2007\n');
    fprintf('=============================================================\n\n');
    
    n_subjects = 500;
    sim_duration = 28;
    
    guide_params = base_params;
    guide_params.Pain_0 = 39.5;
    
    % Placebo arm
    fprintf('Simulating Placebo arm...\n');
    guide_placebo_pain = zeros(sim_duration + 1, n_subjects);
    for i = 1:n_subjects
        ind_params = generate_individual_params(guide_params, iiv_structure);
        ind_params.Pain_0 = guide_params.Pain_0;
        sim = simulate_glucosamine(ind_params, 0, sim_duration, 'sulfate', false);
        guide_placebo_pain(:, i) = sim.Pain;
    end
    
    % GS 1500mg arm
    fprintf('Simulating GS 1500mg arm...\n');
    guide_gs_pain = zeros(sim_duration + 1, n_subjects);
    for i = 1:n_subjects
        ind_params = generate_individual_params(guide_params, iiv_structure);
        ind_params.Pain_0 = guide_params.Pain_0;
        sim = simulate_glucosamine(ind_params, 1500, sim_duration, 'sulfate', true);
        guide_gs_pain(:, i) = sim.Pain;
    end
    
    % Summarize predictions
    guide_pred_placebo = table((0:sim_duration)', repmat({'Placebo'}, sim_duration+1, 1), ...
                               mean(guide_placebo_pain, 2), ...
                               prctile(guide_placebo_pain, 5, 2), ...
                               prctile(guide_placebo_pain, 95, 2), ...
                               'VariableNames', {'week', 'arm', 'pred_mean', 'pred_p5', 'pred_p95'});
    
    guide_pred_gs = table((0:sim_duration)', repmat({'GS 1500mg'}, sim_duration+1, 1), ...
                          mean(guide_gs_pain, 2), ...
                          prctile(guide_gs_pain, 5, 2), ...
                          prctile(guide_gs_pain, 95, 2), ...
                          'VariableNames', {'week', 'arm', 'pred_mean', 'pred_p5', 'pred_p95'});
    
    guide_predictions = [guide_pred_placebo; guide_pred_gs];
    
    %% =============================================================================
    %% 4. SIMULATE MOVES 2015 TRIAL
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATING MOVES 2015\n');
    fprintf('=============================================================\n\n');
    
    moves_params = base_params;
    moves_params.Pain_0 = 74.4;
    
    fprintf('Simulating CS+GH arm...\n');
    moves_pain = zeros(sim_duration + 1, n_subjects);
    for i = 1:n_subjects
        ind_params = generate_individual_params(moves_params, iiv_structure);
        ind_params.Pain_0 = moves_params.Pain_0;
        sim = simulate_glucosamine(ind_params, 1500, sim_duration, 'HCl', true);
        moves_pain(:, i) = sim.Pain;
    end
    
    moves_predictions = table((0:sim_duration)', repmat({'CS+GH'}, sim_duration+1, 1), ...
                              mean(moves_pain, 2), ...
                              prctile(moves_pain, 5, 2), ...
                              prctile(moves_pain, 95, 2), ...
                              'VariableNames', {'week', 'arm', 'pred_mean', 'pred_p5', 'pred_p95'});
    
    %% =============================================================================
    %% 5. VALIDATION METRICS - GUIDE 2007
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('VALIDATION METRICS: GUIDE 2007\n');
    fprintf('=============================================================\n\n');
    
    guide_pred_placebo_w24 = guide_pred_placebo(guide_pred_placebo.week == 24, :);
    guide_pred_gs_w24 = guide_pred_gs(guide_pred_gs.week == 24, :);
    
    guide_validation = table({'Placebo'; 'GS 1500mg'}, [24; 24], ...
                             [guide_observed.pain_mean(strcmp(guide_observed.arm, 'Placebo') & guide_observed.week == 24); ...
                              guide_observed.pain_mean(strcmp(guide_observed.arm, 'GS 1500mg') & guide_observed.week == 24)], ...
                             [guide_pred_placebo_w24.pred_mean; guide_pred_gs_w24.pred_mean], ...
                             [guide_pred_placebo_w24.pred_p5; guide_pred_gs_w24.pred_p5], ...
                             [guide_pred_placebo_w24.pred_p95; guide_pred_gs_w24.pred_p95], ...
                             'VariableNames', {'arm', 'week', 'observed', 'predicted', 'pred_p5', 'pred_p95'});
    
    guide_validation.error = guide_validation.predicted - guide_validation.observed;
    guide_validation.abs_error = abs(guide_validation.error);
    guide_validation.within_PI = guide_validation.observed >= guide_validation.pred_p5 & ...
                                  guide_validation.observed <= guide_validation.pred_p95;
    
    guide_te_predicted = guide_pred_gs_w24.pred_mean - guide_pred_placebo_w24.pred_mean;
    guide_te_observed = guide_validation.observed(2) - guide_validation.observed(1);
    
    guide_mae = mean(guide_validation.abs_error);
    guide_rmse = sqrt(mean(guide_validation.error.^2));
    guide_coverage = mean(guide_validation.within_PI);
    
    fprintf('GUIDE 2007 Results:\n');
    fprintf('  MAE:  %.1f points\n', guide_mae);
    fprintf('  RMSE: %.1f points\n', guide_rmse);
    fprintf('  90%% PI Coverage: %.0f%%\n', guide_coverage * 100);
    fprintf('  Treatment Effect (Pred): %.1f points\n', guide_te_predicted);
    fprintf('  Treatment Effect (Obs):  %.1f points\n', guide_te_observed);
    
    %% =============================================================================
    %% 6. VALIDATION METRICS - MOVES 2015
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('VALIDATION METRICS: MOVES 2015\n');
    fprintf('=============================================================\n\n');
    
    moves_obs_cs_gh = moves_observed(strcmp(moves_observed.arm, 'CS+GH'), :);
    validation_weeks = [0, 4, 9, 17, 26];
    
    moves_validation = table('Size', [length(validation_weeks), 7], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'day', 'week', 'observed', 'predicted', 'pred_p5', 'pred_p95', 'within_PI'});
    
    for i = 1:length(validation_weeks)
        wk = validation_weeks(i);
        obs_idx = find(moves_obs_cs_gh.week == wk);
        pred_idx = find(moves_predictions.week == wk);
        
        moves_validation.day(i) = moves_obs_cs_gh.day(obs_idx);
        moves_validation.week(i) = wk;
        moves_validation.observed(i) = moves_obs_cs_gh.pain_mean(obs_idx);
        moves_validation.predicted(i) = moves_predictions.pred_mean(pred_idx);
        moves_validation.pred_p5(i) = moves_predictions.pred_p5(pred_idx);
        moves_validation.pred_p95(i) = moves_predictions.pred_p95(pred_idx);
        moves_validation.within_PI(i) = moves_validation.observed(i) >= moves_validation.pred_p5(i) & ...
                                        moves_validation.observed(i) <= moves_validation.pred_p95(i);
    end
    
    moves_validation.error = moves_validation.predicted - moves_validation.observed;
    moves_validation.abs_error = abs(moves_validation.error);
    
    moves_mae = mean(moves_validation.abs_error);
    moves_rmse = sqrt(mean(moves_validation.error.^2));
    moves_coverage = mean(moves_validation.within_PI);
    
    fprintf('MOVES 2015 Results:\n');
    fprintf('  MAE:  %.1f points\n', moves_mae);
    fprintf('  RMSE: %.1f points\n', moves_rmse);
    fprintf('  90%% PI Coverage: %.0f%%\n', moves_coverage * 100);
    
    %% =============================================================================
    %% 7. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    validation_results.guide_observed = guide_observed;
    validation_results.guide_predictions = guide_predictions;
    validation_results.guide_validation = guide_validation;
    validation_results.moves_observed = moves_observed;
    validation_results.moves_predictions = moves_predictions;
    validation_results.moves_validation = moves_validation;
    validation_results.metrics = table({'GUIDE 2007'; 'MOVES 2015'}, ...
                                       [guide_mae; moves_mae], ...
                                       [guide_rmse; moves_rmse], ...
                                       [guide_coverage; moves_coverage], ...
                                       'VariableNames', {'Study', 'MAE', 'RMSE', 'Coverage'});
    
    save('model_output/validation_results.mat', 'validation_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/validation_results.mat\n');
    fprintf('=============================================================\n');
    
end

%% =============================================================================
%% HELPER FUNCTION
%% =============================================================================

function ind_params = generate_individual_params(base_params, iiv_structure)
    % Generate individual parameters with IIV
    
    ind_params = base_params;
    
    for j = 1:height(iiv_structure)
        pname = iiv_structure.parameter{j};
        cv = iiv_structure.cv_percent(j) / 100;
        dist = iiv_structure.distribution{j};
        
        if isfield(base_params, pname)
            base_val = base_params.(pname);
            if strcmp(dist, 'lognormal')
                omega = sqrt(log(1 + cv^2));
                eta = randn() * omega;
                ind_params.(pname) = base_val * exp(eta);
            else
                sd_val = base_val * cv;
                ind_params.(pname) = base_val + randn() * sd_val;
            end
        end
    end
    
    % Apply constraints
    ind_params.Imax_deg = min(max(ind_params.Imax_deg, 0.05), 1.00);
    ind_params.GAG_0 = min(max(ind_params.GAG_0, 0.70), 1.05);
    ind_params.Pain_0 = min(max(ind_params.Pain_0, 20), 80);
    
    if ind_params.k_deg <= ind_params.k_syn
        ind_params.k_deg = ind_params.k_syn * 1.2;
    end
end
