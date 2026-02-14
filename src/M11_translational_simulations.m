%% =============================================================================
%% M11_translational_simulations.m
%% Translational Simulations for Clinical Decision Support
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB Implementation
%%
%% CRITICAL: This script calls the external simulate_glucosamine() function
%% (the SAME function used for calibration in M04 and virtual population in
%% M09) to ensure model consistency.
%%
%% IIV APPROACH:
%%   Calibration (M09): Uses trial-specific variability (Table S1 CVs)
%%   Translational (M11): Uses EXPANDED variability reflecting real-world
%%   population heterogeneity beyond controlled trial settings.
%%   This is standard QSP practice — translational predictions should account
%%   for broader interindividual variability than observed in controlled trials.
%%   See manuscript Section 3.5: "expanded population-variability assumptions"
%%
%% Simulations:
%%   A. Dose-Response Analysis (250-3000 mg/day)
%%   B. Enhanced Bioavailability Formulations (22-60%)
%%   C. Treatment Duration Optimization (6 months to 5 years)
%%   D. Patient Stratification (Early vs Late OA)
%%
%% Output: model_output/translational_simulations.mat
%%
%% =============================================================================

function trans_results = M11_translational_simulations()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('TRANSLATIONAL SIMULATIONS\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD MODEL AND DEFINE EXPANDED IIV
    %% =============================================================================
    
    fprintf('Loading calibrated model...\n');
    
    load('model_output/virtual_population_results.mat', 'vpop_results');
    base_params = vpop_results.calibrated_params;
    
    fprintf('  Current bioavailability: F_sulfate = %.0f%%\n', base_params.F_sulfate * 100);
    fprintf('  Imax_deg = %.2f\n', base_params.Imax_deg);
    fprintf('  Using external simulate_glucosamine() for consistency\n');
    
    % -------------------------------------------------------------------------
    % EXPANDED IIV for translational predictions
    % -------------------------------------------------------------------------
    % Rationale: Controlled trials (Reginster, Pavelka) enrolled selected
    % populations with specific inclusion criteria. Real-world patients show
    % greater heterogeneity in drug metabolism, disease severity, and placebo
    % response. CVs are expanded ~1.5x from Table S1 calibration values to
    % reflect this broader population variability.
    %
    %   Parameter       Calibration CV   Translational CV   Rationale
    %   ---------       --------------   ----------------   ---------
    %   k_syn           40%              55%                Wider turnover variability
    %   k_deg           40%              55%                Wider turnover variability
    %   Imax_deg        30%              50%                Drug potency heterogeneity
    %   GAG_0           15%              25%                Broader disease staging
    %   Pain_0          25%              35%                Real-world pain variability
    %   Pmax            50%              65%                Placebo response heterogeneity
    %   kpl             30%              45%                Placebo onset variability
    %   F_sulfate       35%              55%                Absorption/formulation variability
    %   F_HCl           35%              55%                Absorption/formulation variability
    % -------------------------------------------------------------------------
    
    trans_iiv_params = {'k_syn', 'k_deg', 'Imax_deg', 'GAG_0', 'Pain_0', 'Pmax', 'kpl', 'F_sulfate', 'F_HCl'};
    trans_iiv_cv     = [55,      55,      50,         25,      35,       65,     45,    55,         55];
    trans_iiv_dist   = {'lognormal','lognormal','lognormal','normal','normal','lognormal','lognormal','lognormal','lognormal'};
    
    trans_iiv = table(trans_iiv_params', trans_iiv_cv', trans_iiv_dist', ...
                      'VariableNames', {'parameter', 'cv_percent', 'distribution'});
    
    fprintf('\n  Translational IIV (expanded from Table S1):\n');
    for i = 1:height(trans_iiv)
        fprintf('    %-12s: CV = %2d%% (%s)\n', ...
                trans_iiv.parameter{i}, trans_iiv.cv_percent(i), trans_iiv.distribution{i});
    end
    
    n_subjects = 500;  % Larger population for better tail sampling
    fprintf('\n  Population size: N = %d per scenario\n', n_subjects);
    
    %% =============================================================================
    %% 2. SIMULATION A: DOSE-RESPONSE
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATION A: DOSE-RESPONSE ANALYSIS\n');
    fprintf('=============================================================\n\n');
    
    doses = [250, 500, 750, 1000, 1500, 2000, 2500, 3000];
    
    fprintf('Simulating %d dose levels x %d subjects...\n', length(doses), n_subjects);
    
    dose_response_results = table('Size', [length(doses), 6], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'dose', 'jsw_te_mean', 'jsw_te_sd', 'jsw_te_p5', 'jsw_te_p95', 'responder_pct'});
    
    for d = 1:length(doses)
        dose = doses(d);
        jsw_te = zeros(n_subjects, 1);
        
        for i = 1:n_subjects
            ind_params = generate_translational_params(base_params, trans_iiv);
            
            sim_placebo = simulate_glucosamine(ind_params, dose, 156, 'sulfate', false);
            sim_treat   = simulate_glucosamine(ind_params, dose, 156, 'sulfate', true);
            
            jsw_te(i) = sim_treat.JSW_change(end) - sim_placebo.JSW_change(end);
        end
        
        dose_response_results.dose(d) = dose;
        dose_response_results.jsw_te_mean(d) = mean(jsw_te);
        dose_response_results.jsw_te_sd(d) = std(jsw_te);
        dose_response_results.jsw_te_p5(d) = prctile(jsw_te, 5);
        dose_response_results.jsw_te_p95(d) = prctile(jsw_te, 95);
        dose_response_results.responder_pct(d) = mean(jsw_te > 0.1) * 100;
        
        fprintf('  %4d mg/day: TE = %.3f mm (%.0f%% responders)\n', ...
                dose, mean(jsw_te), mean(jsw_te > 0.1) * 100);
    end
    
    fprintf('\nDose-Response Summary:\n');
    fprintf('%-10s %12s %12s\n', 'Dose (mg)', 'Mean TE (mm)', 'Responders');
    fprintf('%s\n', repmat('-', 1, 40));
    for i = 1:height(dose_response_results)
        fprintf('%-10d %12.3f %11.0f%%\n', ...
                dose_response_results.dose(i), ...
                dose_response_results.jsw_te_mean(i), ...
                dose_response_results.responder_pct(i));
    end
    
    %% =============================================================================
    %% 3. SIMULATION B: BIOAVAILABILITY ENHANCEMENT
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATION B: ENHANCED BIOAVAILABILITY\n');
    fprintf('=============================================================\n\n');
    
    F_values = [0.22, 0.26, 0.35, 0.40, 0.50, 0.60];
    
    fprintf('Simulating %d bioavailability scenarios...\n', length(F_values));
    
    bioav_results = table('Size', [length(F_values), 6], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'F_percent', 'jsw_te_mean', 'jsw_te_sd', 'jsw_te_p5', 'jsw_te_p95', 'responder_pct'});
    
    for f = 1:length(F_values)
        F_val = F_values(f);
        jsw_te = zeros(n_subjects, 1);
        
        for i = 1:n_subjects
            ind_params = generate_translational_params(base_params, trans_iiv);
            
            % Override F_sulfate — this is a fixed formulation property
            ind_params.F_sulfate = F_val;
            
            sim_placebo = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', false);
            sim_treat   = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', true);
            
            jsw_te(i) = sim_treat.JSW_change(end) - sim_placebo.JSW_change(end);
        end
        
        bioav_results.F_percent(f) = F_val * 100;
        bioav_results.jsw_te_mean(f) = mean(jsw_te);
        bioav_results.jsw_te_sd(f) = std(jsw_te);
        bioav_results.jsw_te_p5(f) = prctile(jsw_te, 5);
        bioav_results.jsw_te_p95(f) = prctile(jsw_te, 95);
        bioav_results.responder_pct(f) = mean(jsw_te > 0.1) * 100;
        
        fprintf('  F = %2.0f%%: TE = %.3f mm (%.0f%% responders)\n', ...
                F_val * 100, mean(jsw_te), mean(jsw_te > 0.1) * 100);
    end
    
    fprintf('\nBioavailability Enhancement Summary:\n');
    fprintf('%-10s %12s %12s\n', 'F (%)', 'TE (mm)', 'Responders');
    fprintf('%s\n', repmat('-', 1, 40));
    for i = 1:height(bioav_results)
        fprintf('%-10.0f %12.3f %11.0f%%\n', ...
                bioav_results.F_percent(i), ...
                bioav_results.jsw_te_mean(i), ...
                bioav_results.responder_pct(i));
    end
    
    %% =============================================================================
    %% 4. SIMULATION C: TREATMENT DURATION
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATION C: TREATMENT DURATION\n');
    fprintf('=============================================================\n\n');
    
    duration_weeks = 260;  % 5 years
    
    fprintf('Simulating 5-year treatment courses (%d subjects)...\n', n_subjects);
    
    time_results = zeros(duration_weeks + 1, n_subjects);
    
    for i = 1:n_subjects
        ind_params = generate_translational_params(base_params, trans_iiv);
        
        sim_placebo = simulate_glucosamine(ind_params, 1500, duration_weeks, 'sulfate', false);
        sim_treat   = simulate_glucosamine(ind_params, 1500, duration_weeks, 'sulfate', true);
        
        time_results(:, i) = sim_treat.JSW_change - sim_placebo.JSW_change;
        
        if mod(i, 100) == 0, fprintf('.'); end
    end
    fprintf(' Done!\n');
    
    eval_weeks = [26, 52, 104, 156, 208, 260];
    yearly_summary = table('Size', [length(eval_weeks), 5], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'week', 'year', 'jsw_te_mean', 'jsw_te_p5', 'jsw_te_p95'});
    
    for w = 1:length(eval_weeks)
        wk = eval_weeks(w);
        yearly_summary.week(w) = wk;
        yearly_summary.year(w) = wk / 52;
        yearly_summary.jsw_te_mean(w) = mean(time_results(wk + 1, :));
        yearly_summary.jsw_te_p5(w) = prctile(time_results(wk + 1, :), 5);
        yearly_summary.jsw_te_p95(w) = prctile(time_results(wk + 1, :), 95);
    end
    
    fprintf('\nTreatment Duration Summary:\n');
    fprintf('%-12s %12s\n', 'Duration', 'JSW TE (mm)');
    fprintf('%s\n', repmat('-', 1, 30));
    for i = 1:height(yearly_summary)
        if yearly_summary.week(i) == 26
            dur_label = '6 months';
        else
            dur_label = sprintf('Year %.0f', yearly_summary.year(i));
        end
        fprintf('%-12s %12.3f\n', dur_label, yearly_summary.jsw_te_mean(i));
    end
    
    %% =============================================================================
    %% 5. SIMULATION D: PATIENT STRATIFICATION
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATION D: PATIENT STRATIFICATION\n');
    fprintf('=============================================================\n\n');
    
    subgroups = {'Very Early OA', 'Early OA', 'Moderate OA', 'Late OA', 'Severe OA'};
    GAG_0_values = [0.98, 0.95, 0.90, 0.85, 0.80];
    
    fprintf('Simulating %d patient subgroups...\n', length(subgroups));
    
    strat_results = table('Size', [length(subgroups), 8], ...
        'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'subgroup', 'GAG_0', 'jsw_te_mean', 'jsw_te_sd', 'jsw_te_p5', 'jsw_te_p95', 'responder_pct', 'NNT'});
    
    for s = 1:length(subgroups)
        jsw_te = zeros(n_subjects, 1);
        
        for i = 1:n_subjects
            ind_params = generate_translational_params(base_params, trans_iiv);
            
            % Fix GAG_0 for this subgroup
            ind_params.GAG_0 = GAG_0_values(s);
            
            sim_placebo = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', false);
            sim_treat   = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', true);
            
            jsw_te(i) = sim_treat.JSW_change(end) - sim_placebo.JSW_change(end);
        end
        
        responder_rate = mean(jsw_te > 0.1);
        
        strat_results.subgroup(s) = subgroups{s};
        strat_results.GAG_0(s) = GAG_0_values(s) * 100;
        strat_results.jsw_te_mean(s) = mean(jsw_te);
        strat_results.jsw_te_sd(s) = std(jsw_te);
        strat_results.jsw_te_p5(s) = prctile(jsw_te, 5);
        strat_results.jsw_te_p95(s) = prctile(jsw_te, 95);
        strat_results.responder_pct(s) = responder_rate * 100;
        if responder_rate > 0
            strat_results.NNT(s) = ceil(1 / responder_rate);
        else
            strat_results.NNT(s) = NaN;
        end
        
        fprintf('  %-15s (GAG_0=%2.0f%%): TE = %.3f mm, %.0f%% resp, NNT = %d\n', ...
                subgroups{s}, GAG_0_values(s)*100, mean(jsw_te), responder_rate*100, strat_results.NNT(s));
    end
    
    fprintf('\nPatient Stratification Summary:\n');
    fprintf('%-15s %8s %12s %10s %6s\n', 'Subgroup', 'GAG_0', 'TE (mm)', 'Responders', 'NNT');
    fprintf('%s\n', repmat('-', 1, 55));
    for i = 1:height(strat_results)
        fprintf('%-15s %7.0f%% %12.3f %9.0f%% %6d\n', ...
                strat_results.subgroup{i}, ...
                strat_results.GAG_0(i), ...
                strat_results.jsw_te_mean(i), ...
                strat_results.responder_pct(i), ...
                strat_results.NNT(i));
    end
    
    %% =============================================================================
    %% 6. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    trans_results.dose_response = dose_response_results;
    trans_results.bioavailability = bioav_results;
    trans_results.duration = yearly_summary;
    trans_results.stratification = strat_results;
    trans_results.time_courses = time_results;
    trans_results.translational_iiv = trans_iiv;
    
    save('model_output/translational_simulations.mat', 'trans_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/translational_simulations.mat\n');
    fprintf('=============================================================\n');
    
end

%% =============================================================================
%% HELPER: Generate individual parameters with EXPANDED translational IIV
%% =============================================================================
%%
%% Uses the same sampling methodology as M09 (lognormal/normal with constraints)
%% but with WIDER CVs reflecting real-world population heterogeneity.
%% Constraints are aligned with M09 generate_virtual_population().
%% =============================================================================

function ind_params = generate_translational_params(base_params, trans_iiv)
    
    ind_params = base_params;
    
    % Apply expanded IIV
    for j = 1:height(trans_iiv)
        pname = trans_iiv.parameter{j};
        cv = trans_iiv.cv_percent(j) / 100;
        dist = trans_iiv.distribution{j};
        
        if isfield(base_params, pname)
            base_val = base_params.(pname);
            if strcmp(dist, 'lognormal')
                omega = sqrt(log(1 + cv^2));
                eta = randn() * omega;
                ind_params.(pname) = base_val * exp(eta);
            else  % normal
                sd_val = base_val * cv;
                ind_params.(pname) = base_val + randn() * sd_val;
            end
        end
    end
    
    % Apply constraints — ALIGNED WITH M09
    ind_params.Imax_deg = min(max(ind_params.Imax_deg, 0.05), 1.00);
    ind_params.GAG_0 = min(max(ind_params.GAG_0, 0.70), 1.05);
    ind_params.Pain_0 = min(max(ind_params.Pain_0, 20), 80);
    ind_params.F_sulfate = min(max(ind_params.F_sulfate, 0.05), 0.50);
    ind_params.F_HCl = min(max(ind_params.F_HCl, 0.03), 0.30);
    
    % Ensure disease progression
    if ind_params.k_deg <= ind_params.k_syn
        ind_params.k_deg = ind_params.k_syn * 1.2;
    end
end
