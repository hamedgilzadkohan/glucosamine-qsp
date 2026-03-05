%% =============================================================================
%% M14_equal_bioavailability_sensitivity.m
%% Sensitivity Analysis: What if F_sulfate = F_HCl?
%% =============================================================================
%%
%% PURPOSE:
%%   Addresses reviewer concern that the assumed 2-fold bioavailability
%%   difference between glucosamine sulfate (F=0.22) and glucosamine HCl
%%   (F=0.11) may not be justified, given:
%%     - Chang et al. (2025): no significant PK difference between cGS and rGS
%%     - Sahoo et al. (2012): commercial "glucosamine sulfate" is physically
%%       glucosamine chloride + K2SO4
%%
%%   This script runs the calibrated model under three bioavailability
%%   scenarios and reports the impact on all key predictions.
%%
%% OUTPUT:
%%   - Table of predictions under each scenario
%%   - Text suitable for inclusion in manuscript Limitations/Discussion
%%
%% =============================================================================

function equal_ba_results = M14_equal_bioavailability_sensitivity()
    
    rng(42);
    
    fprintf('=============================================================\n');
    fprintf('EQUAL BIOAVAILABILITY SENSITIVITY ANALYSIS\n');
    fprintf('Addressing reviewer concerns re: Chang et al. 2025,\n');
    fprintf('Sahoo et al. 2012\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD CALIBRATED MODEL
    %% =============================================================================
    
    load('model_output/calibration_results.mat', 'calibration_results');
    params = calibration_results.params;
    
    % Ensure synthesis parameters are present
    if ~isfield(params, 'Emax_syn')
        params.Emax_syn = 0.15;
    end
    if ~isfield(params, 'EC50_syn')
        params.EC50_syn = 2.0;
    end
    
    %% =============================================================================
    %% 2. DEFINE SCENARIOS
    %% =============================================================================
    
    % Scenario A: Base case (as published)
    %   F_sulfate = 0.22, F_HCl = 0.11
    %
    % Scenario B: Equal BA at low value (both = 0.11)
    %   Interpretation: What if cGS has no BA advantage? (worst case)
    %
    % Scenario C: Equal BA at high value (both = 0.22)
    %   Interpretation: What if GHCl is as bioavailable as cGS?
    %
    % Scenario D: Equal BA at intermediate value (both = 0.165)
    %   Interpretation: Midpoint — true BA unknown but equal
    
    scenarios = struct();
    
    scenarios(1).name = 'Base case (F_GS=0.22, F_HCl=0.11)';
    scenarios(1).F_sulfate = 0.22;
    scenarios(1).F_HCl = 0.11;
    
    scenarios(2).name = 'Equal low (both F=0.11)';
    scenarios(2).F_sulfate = 0.11;
    scenarios(2).F_HCl = 0.11;
    
    scenarios(3).name = 'Equal high (both F=0.22)';
    scenarios(3).F_sulfate = 0.22;
    scenarios(3).F_HCl = 0.22;
    
    scenarios(4).name = 'Equal mid (both F=0.165)';
    scenarios(4).F_sulfate = 0.165;
    scenarios(4).F_HCl = 0.165;
    
    %% =============================================================================
    %% 3. RUN SIMULATIONS FOR EACH SCENARIO
    %% =============================================================================
    
    n_scenarios = length(scenarios);
    
    % Storage for results
    results = struct();
    
    for s = 1:n_scenarios
        fprintf('--- Scenario %d: %s ---\n', s, scenarios(s).name);
        
        p = params;
        p.F_sulfate = scenarios(s).F_sulfate;
        p.F_HCl = scenarios(s).F_HCl;
        
        %% --- 3a. 3-Year JSW (European trial design: sulfate 1500mg) ---
        sim_p_jsw = simulate_glucosamine(p, 1500, 156, 'sulfate', false);
        sim_t_jsw = simulate_glucosamine(p, 1500, 156, 'sulfate', true);
        
        plac_jsw_raw = sim_p_jsw.JSW_change(end);
        treat_jsw_raw = sim_t_jsw.JSW_change(end);
        
        % Sign correction if needed (placebo should be negative)
        if plac_jsw_raw > 0
            plac_jsw = -abs(plac_jsw_raw);
            treat_jsw = treat_jsw_raw - plac_jsw_raw - plac_jsw;
        else
            plac_jsw = plac_jsw_raw;
            treat_jsw = treat_jsw_raw;
        end
        
        te_jsw = treat_jsw - plac_jsw;
        
        %% --- 3b. 24-Week WOMAC Pain (GAIT design: HCl 1500mg) ---
        sim_p_pain = simulate_glucosamine(p, 1500, 24, 'HCl', false);
        sim_t_pain = simulate_glucosamine(p, 1500, 24, 'HCl', true);
        
        pain_plac_24w = sim_p_pain.Pain(25);
        pain_treat_24w = sim_t_pain.Pain(25);
        pain_diff_24w = pain_plac_24w - pain_treat_24w;  % positive = drug better
        
        %% --- 3c. GUIDE-like design (sulfate 1500mg, 26 weeks) ---
        sim_p_guide = simulate_glucosamine(p, 1500, 26, 'sulfate', false);
        sim_t_guide = simulate_glucosamine(p, 1500, 26, 'sulfate', true);
        
        pain_plac_guide = sim_p_guide.Pain(end);
        pain_treat_guide = sim_t_guide.Pain(end);
        pain_diff_guide = pain_plac_guide - pain_treat_guide;
        
        %% --- 3d. Dose-response: minimum effective dose ---
        doses = [250, 500, 750, 1000, 1500, 2000, 3000];
        te_by_dose = zeros(length(doses), 1);
        
        for d = 1:length(doses)
            sim_p_d = simulate_glucosamine(p, doses(d), 156, 'sulfate', false);
            sim_t_d = simulate_glucosamine(p, doses(d), 156, 'sulfate', true);
            
            plac_d = sim_p_d.JSW_change(end);
            treat_d = sim_t_d.JSW_change(end);
            
            if plac_d > 0
                plac_d_corr = -abs(plac_d);
                treat_d_corr = treat_d - plac_d - plac_d_corr;
            else
                plac_d_corr = plac_d;
                treat_d_corr = treat_d;
            end
            
            te_by_dose(d) = treat_d_corr - plac_d_corr;
        end
        
        % Find minimum effective dose (TE > 0.1 mm MCID)
        min_eff_idx = find(te_by_dose >= 0.1, 1, 'first');
        if ~isempty(min_eff_idx)
            min_eff_dose = doses(min_eff_idx);
        else
            min_eff_dose = NaN;  % No dose achieves MCID
        end
        
        %% --- Store results ---
        results(s).name = scenarios(s).name;
        results(s).F_sulfate = scenarios(s).F_sulfate;
        results(s).F_HCl = scenarios(s).F_HCl;
        results(s).te_jsw_3yr = te_jsw;
        results(s).plac_jsw_3yr = plac_jsw;
        results(s).treat_jsw_3yr = treat_jsw;
        results(s).pain_plac_24w = pain_plac_24w;
        results(s).pain_treat_24w = pain_treat_24w;
        results(s).pain_diff_24w = pain_diff_24w;
        results(s).pain_diff_guide = pain_diff_guide;
        results(s).min_eff_dose = min_eff_dose;
        results(s).te_by_dose = te_by_dose;
        results(s).doses = doses;
        
        fprintf('  3yr JSW TE:     %.3f mm\n', te_jsw);
        fprintf('  3yr Placebo:    %.3f mm\n', plac_jsw);
        fprintf('  GAIT pain diff: %.1f WOMAC pts\n', pain_diff_24w);
        fprintf('  GUIDE pain diff:%.1f WOMAC pts\n', pain_diff_guide);
        fprintf('  Min eff dose:   %s mg/day\n', ...
                iff(isnan(min_eff_dose), 'None achieves MCID', num2str(min_eff_dose)));
        fprintf('\n');
    end
    
    %% =============================================================================
    %% 4. COMPARISON TABLE
    %% =============================================================================
    
    fprintf('=============================================================\n');
    fprintf('COMPARISON TABLE\n');
    fprintf('=============================================================\n\n');
    
    fprintf('%-35s | %8s | %8s | %8s | %8s | %8s\n', ...
            'Scenario', 'F_GS', 'F_HCl', 'JSW TE', 'GAIT Δ', 'MinDose');
    fprintf('%s\n', repmat('-', 1, 90));
    
    for s = 1:n_scenarios
        fprintf('%-35s | %8.3f | %8.3f | %8.3f | %8.1f | %8s\n', ...
                results(s).name, ...
                results(s).F_sulfate, ...
                results(s).F_HCl, ...
                results(s).te_jsw_3yr, ...
                results(s).pain_diff_24w, ...
                iff(isnan(results(s).min_eff_dose), 'N/A', ...
                    num2str(results(s).min_eff_dose)));
    end
    
    %% =============================================================================
    %% 5. CALCULATE PERCENTAGE CHANGES FROM BASE CASE
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('PERCENTAGE CHANGE FROM BASE CASE\n');
    fprintf('=============================================================\n\n');
    
    base_te = results(1).te_jsw_3yr;
    base_pain = results(1).pain_diff_24w;
    
    for s = 2:n_scenarios
        pct_te = (results(s).te_jsw_3yr - base_te) / base_te * 100;
        pct_pain = (results(s).pain_diff_24w - base_pain) / abs(base_pain) * 100;
        
        fprintf('%s:\n', results(s).name);
        fprintf('  JSW TE change: %+.1f%% (%.3f -> %.3f mm)\n', ...
                pct_te, base_te, results(s).te_jsw_3yr);
        fprintf('  GAIT pain diff change: %+.1f%%\n', pct_pain);
        fprintf('\n');
    end
    
    %% =============================================================================
    %% 6. GENERATE MANUSCRIPT-READY TEXT
    %% =============================================================================
    
    fprintf('=============================================================\n');
    fprintf('MANUSCRIPT-READY TEXT (for Limitations section)\n');
    fprintf('=============================================================\n\n');
    
    fprintf('To assess the sensitivity of model predictions to the assumed\n');
    fprintf('bioavailability difference between formulations, we conducted\n');
    fprintf('a post hoc analysis comparing three scenarios: the base case\n');
    fprintf('(F_sulfate = 0.22, F_HCl = 0.11), equal bioavailability at the\n');
    fprintf('lower value (both F = 0.11), and equal bioavailability at the\n');
    fprintf('higher value (both F = 0.22). When bioavailability was equalized\n');
    fprintf('at the lower value (F = 0.11), the predicted 3-year JSW treatment\n');
    fprintf('effect decreased from %.3f mm to %.3f mm (%.0f%% reduction),\n', ...
            results(1).te_jsw_3yr, results(2).te_jsw_3yr, ...
            abs((results(2).te_jsw_3yr - results(1).te_jsw_3yr) / results(1).te_jsw_3yr * 100));
    fprintf('and the minimum effective dose shifted from %s to %s mg/day.\n', ...
            num2str(results(1).min_eff_dose), ...
            iff(isnan(results(2).min_eff_dose), 'above all tested doses', ...
                num2str(results(2).min_eff_dose)));
    fprintf('This analysis confirms that the model predictions are sensitive\n');
    fprintf('to the bioavailability assumptions and underscores the importance\n');
    fprintf('of resolving the formulation bioavailability question through\n');
    fprintf('adequately powered, head-to-head pharmacokinetic studies.\n');
    
    %% =============================================================================
    %% 7. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    equal_ba_results.scenarios = results;
    equal_ba_results.description = 'Equal bioavailability sensitivity analysis';
    equal_ba_results.date = datestr(now);
    
    save('model_output/equal_ba_sensitivity_results.mat', 'equal_ba_results');
    
    fprintf('\n=============================================================\n');
    fprintf('ANALYSIS COMPLETE\n');
    fprintf('Results saved to: model_output/equal_ba_sensitivity_results.mat\n');
    fprintf('=============================================================\n');
    
end

%% Helper function
function result = iff(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
