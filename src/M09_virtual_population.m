%% =============================================================================
%% M09_virtual_population.m
%% Virtual Population Simulation with Inter-Individual Variability
%% =============================================================================
%%
%%
%% Output: model_output/virtual_population_results.mat
%%
%% =============================================================================

function vpop_results = M09_virtual_population()
    
    rng(42);  % For reproducibility
    
    fprintf('=============================================================\n');
    fprintf('VIRTUAL POPULATION SIMULATION\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. LOAD CALIBRATED MODEL
    %% =============================================================================
    
    fprintf('Loading calibrated model...\n');
    
    load('model_output/calibration_results.mat', 'calibration_results');
    calibrated_params = calibration_results.params;
    
    fprintf('  Calibrated parameters loaded\n');
    fprintf('    k_syn = %.4f /year\n', calibrated_params.k_syn);
    fprintf('    k_deg = %.4f /year\n', calibrated_params.k_deg);
    fprintf('    Imax_deg = %.4f\n', calibrated_params.Imax_deg);
    
    %% =============================================================================
    %% 2. IIV STRUCTURE
    %% =============================================================================
    
    fprintf('\n--- IIV Structure ---\n\n');
    
    iiv_params = {'k_syn', 'k_deg', 'Imax_deg', 'GAG_0', 'Pain_0', 'Pmax', 'F_sulfate', 'F_HCl'};
    iiv_cv = [40, 40, 30, 15, 25, 50, 35, 35];  % CV in percent
    iiv_dist = {'lognormal', 'lognormal', 'lognormal', 'normal', 'normal', 'lognormal', 'lognormal', 'lognormal'};
    
    iiv_structure = table(iiv_params', iiv_cv', iiv_dist', ...
                          'VariableNames', {'parameter', 'cv_percent', 'distribution'});
    
    for i = 1:height(iiv_structure)
        fprintf('  %-12s: CV = %2d%% (%s)\n', ...
                iiv_structure.parameter{i}, ...
                iiv_structure.cv_percent(i), ...
                iiv_structure.distribution{i});
    end
    
    %% =============================================================================
    %% 3. GENERATE VIRTUAL POPULATION
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('GENERATING VIRTUAL POPULATION\n');
    fprintf('=============================================================\n\n');
    
    n_subjects = 500;
    fprintf('Generating %d virtual subjects...\n', n_subjects);
    
    % Generate individual parameters
    vpop_params = generate_virtual_population(calibrated_params, iiv_structure, n_subjects);
    
    fprintf('\nVirtual population parameter summary:\n');
    fprintf('%-12s %10s %10s %10s %10s\n', 'Parameter', 'Mean', 'SD', 'Min', 'Max');
    fprintf('%s\n', repmat('-', 1, 55));
    
    summary_params = {'k_syn', 'k_deg', 'Imax_deg', 'GAG_0', 'Pain_0', 'Pmax'};
    for p = 1:length(summary_params)
        pname = summary_params{p};
        vals = [vpop_params.(pname)];
        fprintf('%-12s %10.4f %10.4f %10.4f %10.4f\n', ...
                pname, mean(vals), std(vals), min(vals), max(vals));
    end
    
    %% =============================================================================
    %% 4. SIMULATE VIRTUAL POPULATION - 3-YEAR JSW (SULFATE)
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('SIMULATING VIRTUAL POPULATION\n');
    fprintf('=============================================================\n\n');
    
    fprintf('Simulating %d subjects (3-year JSW)...\n', n_subjects);
    
    all_sims_placebo_3yr = cell(n_subjects, 1);
    all_sims_treat_3yr = cell(n_subjects, 1);
    
    fprintf('Progress: ');
    pb_interval = ceil(n_subjects / 20);
    
    for i = 1:n_subjects
        ind_params = vpop_params(i);
        
        all_sims_placebo_3yr{i} = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', false);
        all_sims_treat_3yr{i} = simulate_glucosamine(ind_params, 1500, 156, 'sulfate', true);
        
        if mod(i, pb_interval) == 0
            fprintf('.');
        end
    end
    fprintf(' Done!\n');
    
    %% =============================================================================
    %% 4b. SIMULATE VIRTUAL POPULATION - 24-WEEK PAIN (HCl)
    %%     NEW: Required for proper Figure 3C generation
    %% =============================================================================
    
    fprintf('Simulating %d subjects (24-week pain, HCl)...\n', n_subjects);
    
    all_sims_placebo_pain = cell(n_subjects, 1);
    all_sims_treat_pain = cell(n_subjects, 1);
    
    fprintf('Progress: ');
    
    for i = 1:n_subjects
        ind_params = vpop_params(i);
        
        all_sims_placebo_pain{i} = simulate_glucosamine(ind_params, 1500, 24, 'HCl', false);
        all_sims_treat_pain{i} = simulate_glucosamine(ind_params, 1500, 24, 'HCl', true);
        
        if mod(i, pb_interval) == 0
            fprintf('.');
        end
    end
    fprintf(' Done!\n');
    
    %% =============================================================================
    %% 5. CALCULATE SUMMARY STATISTICS
    %% =============================================================================
    
    fprintf('\n--- Calculating Summary Statistics ---\n');
    
    jsw_summary_placebo = summarize_sims(all_sims_placebo_3yr, 'JSW_change');
    jsw_summary_treatment = summarize_sims(all_sims_treat_3yr, 'JSW_change');
    gag_summary_placebo = summarize_sims(all_sims_placebo_3yr, 'GAG');
    gag_summary_treatment = summarize_sims(all_sims_treat_3yr, 'GAG');
    
    % NEW: Pain trajectory summaries for Figure 3C
    pain_summary_placebo = summarize_sims(all_sims_placebo_pain, 'Pain');
    pain_summary_treatment = summarize_sims(all_sims_treat_pain, 'Pain');
    
    %% =============================================================================
    %% 6. KEY RESULTS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('KEY RESULTS WITH 90%% PREDICTION INTERVALS\n');
    fprintf('=============================================================\n\n');
    
    y3_idx = 157;  % Week 156
    
    fprintf('3-Year JSW Change (mm):\n');
    fprintf('  Placebo:   %.3f [%.3f, %.3f]\n', ...
            jsw_summary_placebo.median(y3_idx), ...
            jsw_summary_placebo.p5(y3_idx), ...
            jsw_summary_placebo.p95(y3_idx));
    fprintf('  Treatment: %.3f [%.3f, %.3f]\n', ...
            jsw_summary_treatment.median(y3_idx), ...
            jsw_summary_treatment.p5(y3_idx), ...
            jsw_summary_treatment.p95(y3_idx));
    
    te_median = jsw_summary_treatment.median(y3_idx) - jsw_summary_placebo.median(y3_idx);
    fprintf('  Treatment Effect: %.3f mm\n', te_median);
    
    fprintf('\n3-Year GAG Content (%% of healthy):\n');
    fprintf('  Placebo:   %.1f%% [%.1f%%, %.1f%%]\n', ...
            gag_summary_placebo.median(y3_idx) * 100, ...
            gag_summary_placebo.p5(y3_idx) * 100, ...
            gag_summary_placebo.p95(y3_idx) * 100);
    fprintf('  Treatment: %.1f%% [%.1f%%, %.1f%%]\n', ...
            gag_summary_treatment.median(y3_idx) * 100, ...
            gag_summary_treatment.p5(y3_idx) * 100, ...
            gag_summary_treatment.p95(y3_idx) * 100);
    
    % NEW: Pain summary at week 24
    w24_idx = 25;  % Week 24
    fprintf('\n24-Week Pain Score (WOMAC 0-100):\n');
    fprintf('  Placebo:   %.1f [%.1f, %.1f]\n', ...
            pain_summary_placebo.median(w24_idx), ...
            pain_summary_placebo.p5(w24_idx), ...
            pain_summary_placebo.p95(w24_idx));
    fprintf('  Treatment: %.1f [%.1f, %.1f]\n', ...
            pain_summary_treatment.median(w24_idx), ...
            pain_summary_treatment.p5(w24_idx), ...
            pain_summary_treatment.p95(w24_idx));
    
    %% =============================================================================
    %% 7. RESPONDER ANALYSIS
    %% =============================================================================
    
    fprintf('\n--- Responder Analysis ---\n');
    
    individual_te = zeros(n_subjects, 1);
    for i = 1:n_subjects
        individual_te(i) = all_sims_treat_3yr{i}.JSW_change(y3_idx) - ...
                           all_sims_placebo_3yr{i}.JSW_change(y3_idx);
    end
    
    responders_any = sum(individual_te > 0);
    responders_clinically_meaningful = sum(individual_te > 0.1);
    
    fprintf('  Any benefit (TE > 0):        %d/%d (%.1f%%)\n', ...
            responders_any, n_subjects, responders_any/n_subjects*100);
    fprintf('  Clinically meaningful (>0.1mm): %d/%d (%.1f%%)\n', ...
            responders_clinically_meaningful, n_subjects, ...
            responders_clinically_meaningful/n_subjects*100);
    
    %% =============================================================================
    %% 8. SAVE RESULTS
    %% =============================================================================
    
    if ~exist('model_output', 'dir')
        mkdir('model_output');
    end
    
    vpop_results.calibrated_params = calibrated_params;
    vpop_results.vpop_params = vpop_params;
    vpop_results.iiv_structure = iiv_structure;
    vpop_results.n_subjects = n_subjects;
    vpop_results.jsw_summary_placebo = jsw_summary_placebo;
    vpop_results.jsw_summary_treatment = jsw_summary_treatment;
    vpop_results.gag_summary_placebo = gag_summary_placebo;
    vpop_results.gag_summary_treatment = gag_summary_treatment;
    vpop_results.individual_te = individual_te;
    % NEW: Pain trajectory summaries
    vpop_results.pain_summary_placebo = pain_summary_placebo;
    vpop_results.pain_summary_treatment = pain_summary_treatment;
    vpop_results.responder_summary = table({'Any benefit (TE > 0)'; 'Clinically meaningful (TE > 0.1mm)'}, ...
                                           [responders_any; responders_clinically_meaningful], ...
                                           [responders_any/n_subjects*100; responders_clinically_meaningful/n_subjects*100], ...
                                           'VariableNames', {'definition', 'n_responders', 'pct_responders'});
    
    save('model_output/virtual_population_results.mat', 'vpop_results');
    
    fprintf('\n=============================================================\n');
    fprintf('Results saved to: model_output/virtual_population_results.mat\n');
    fprintf('=============================================================\n');
    
end

%% =============================================================================
%% HELPER FUNCTIONS
%% =============================================================================

function vpop_params = generate_virtual_population(base_params, iiv_structure, n)
    % Generate individual parameter sets for virtual population
    
    % Get field names from base_params
    param_names = fieldnames(base_params);
    
    % Initialize struct array
    vpop_params(n) = base_params;  % Pre-allocate
    
    for i = 1:n
        % Start with base params
        vpop_params(i) = base_params;
        
        % Apply IIV
        for j = 1:height(iiv_structure)
            pname = iiv_structure.parameter{j};
            cv = iiv_structure.cv_percent(j) / 100;
            dist = iiv_structure.distribution{j};
            
            if isfield(base_params, pname)
                base_val = base_params.(pname);
                
                if strcmp(dist, 'lognormal')
                    omega = sqrt(log(1 + cv^2));
                    eta = randn() * omega;
                    vpop_params(i).(pname) = base_val * exp(eta);
                else  % normal
                    sd_val = base_val * cv;
                    vpop_params(i).(pname) = base_val + randn() * sd_val;
                end
            end
        end
        
        % Apply constraints
        vpop_params(i).Imax_deg = min(max(vpop_params(i).Imax_deg, 0.05), 1.00);
        vpop_params(i).GAG_0 = min(max(vpop_params(i).GAG_0, 0.70), 1.05);
        vpop_params(i).Pain_0 = min(max(vpop_params(i).Pain_0, 20), 80);
        vpop_params(i).F_sulfate = min(max(vpop_params(i).F_sulfate, 0.05), 0.50);
        vpop_params(i).F_HCl = min(max(vpop_params(i).F_HCl, 0.03), 0.30);
        
        % k_deg must be > k_syn
        if vpop_params(i).k_deg <= vpop_params(i).k_syn
            vpop_params(i).k_deg = vpop_params(i).k_syn * 1.2;
        end
    end
end

function summary = summarize_sims(sim_list, endpoint)
    % Summarize simulations across virtual population
    
    n_weeks = height(sim_list{1});
    n_sims = length(sim_list);
    
    endpoint_matrix = zeros(n_weeks, n_sims);
    for i = 1:n_sims
        endpoint_matrix(:, i) = sim_list{i}.(endpoint);
    end
    
    summary = table();
    summary.week = sim_list{1}.week;
    summary.mean = mean(endpoint_matrix, 2, 'omitnan');
    summary.median = median(endpoint_matrix, 2, 'omitnan');
    summary.p5 = prctile(endpoint_matrix, 5, 2);
    summary.p95 = prctile(endpoint_matrix, 95, 2);
end
