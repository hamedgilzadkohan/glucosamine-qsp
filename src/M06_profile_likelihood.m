%% =============================================================================
%% M06_profile_likelihood.m
%% Profile Likelihood Analysis for Practical Identifiability (CORRECTED)
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB Implementation
%%
%% Method:
%%   For each estimated parameter, fix it at each value on a grid and
%%   RE-OPTIMIZE the remaining four parameters (fmincon, interior-point),
%%   using the same weighted objective, parameter bounds (M01) and
%%   k_deg > k_syn constraint as calibration and bootstrap (M07).
%%   Delta objective is measured from the global minimum.
%%
%%   A parameter is classed as practically identifiable only if its profile
%%   exceeds the threshold (1.92 = chi-square(1 df, 0.95)/2) on BOTH sides of
%%   the optimum within the evaluated range. Confidence bounds are obtained by
%%   linear interpolation at the threshold crossings.
%%
%%   For comparison, the previous one-at-a-time "slice" (other parameters held
%%   at their calibrated values, no re-optimization; delta measured from the
%%   slice minimum) is computed on the same grid with the same objective.
%%
%% Profiling is done outward from the optimum in each direction, warm-started
%% from the neighbouring solution, with a second start from the optimum; the
%% better of the two fits is kept.
%%
%% Requirements: MATLAB R2023b+, Optimization Toolbox
%% Inputs:  data_clean/calibration_data.mat     (from M03)
%%          model_output/calibration_results.mat (from M04)
%% Outputs: model_output/profile_likelihood_results.mat
%%          model_output/profile_likelihood_summary.csv
%%
%% =============================================================================

function profile_results = M06_profile_likelihood()

    rng(42);

    fprintf('=============================================================\n');
    fprintf('PROFILE LIKELIHOOD ANALYSIS (re-optimized profiles)\n');
    fprintf('=============================================================\n\n');

    %% =========================================================================
    %% 1. LOAD MODEL AND DATA
    %% =========================================================================

    % Explicit outputs from load() are required in functions with nested functions
    S_cal = load('model_output/calibration_results.mat', 'calibration_results');
    calibrated_params = S_cal.calibration_results.params;
    param_info = S_cal.calibration_results.param_info;

    S_dat = load('data_clean/calibration_data.mat', 'calibration_data');
    jsw_data = S_dat.calibration_data.jsw_data;
    womac_data = S_dat.calibration_data.womac_pain_data;

    jsw_targets = jsw_data(jsw_data.week > 0, :);
    womac_targets = womac_data(strcmp(womac_data.study, 'GAIT 2006'), :);

    %% =========================================================================
    %% 2. SETUP
    %% =========================================================================

    par_names = {'k_syn', 'k_deg', 'Imax_deg', 'Pain_0', 'Pmax'};
    n_par = numel(par_names);

    lb = zeros(1, n_par);
    ub = zeros(1, n_par);
    x_cal = zeros(1, n_par);
    for i = 1:n_par
        row = strcmp(param_info.Name, par_names{i});
        lb(i) = param_info.Lower(row);
        ub(i) = param_info.Upper(row);
        x_cal(i) = calibrated_params.(par_names{i});
    end

    i_syn = 1;          % index of k_syn
    i_deg = 2;          % index of k_deg
    eps_k = 1e-4;       % enforces k_deg >= k_syn + eps_k

    chi2_threshold = 1.92;
    n_grid = 21;

    % Grids (Imax_deg range matches the manuscript: 0.5-1.0)
    grids.k_syn    = linspace(0.08, 0.20, n_grid);
    grids.k_deg    = linspace(0.12, 0.25, n_grid);
    grids.Imax_deg = linspace(0.50, 1.00, n_grid);
    grids.Pain_0   = linspace(35, 55, n_grid);
    grids.Pmax     = linspace(10, 30, n_grid);

    opts = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'off', 'MaxIterations', 500, ...
        'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10);

    %% =========================================================================
    %% 3. GLOBAL REFERENCE MINIMUM
    %% =========================================================================

    obj_published = objective_full(x_cal);
    fprintf('Objective at published parameter values: %.4f\n', obj_published);

    [x_best, f_best] = refit(1:n_par, x_cal);
    if isnan(f_best) || f_best >= obj_published
        x_best = x_cal;
        f_best = obj_published;
        fprintf('Re-fit of all 5 parameters did not improve on published values.\n');
    else
        fprintf('Re-fit of all 5 parameters: objective %.4f\n', f_best);
        for i = 1:n_par
            fprintf('  %-9s published %.6f   re-fit %.6f\n', ...
                par_names{i}, x_cal(i), x_best(i));
        end
    end
    fprintf('\n');

    %% =========================================================================
    %% 4. PROFILES (re-optimized) AND SLICES (previous method)
    %% =========================================================================

    profile_data = struct();
    tic;

    for p = 1:n_par
        pname = par_names{p};
        fprintf('Profiling %-9s ', pname);

        grid_vals = unique([grids.(pname), x_best(p)]);           % sorted
        grid_vals = grid_vals(grid_vals >= lb(p) & grid_vals <= ub(p));
        n_g = numel(grid_vals);
        [~, i_opt] = min(abs(grid_vals - x_best(p)));
        free_idx = setdiff(1:n_par, p);

        prof_obj = nan(n_g, 1);
        prof_x = nan(n_g, n_par);
        slice_obj = nan(n_g, 1);

        for dirn = [1, -1]
            if dirn == 1
                idx_seq = i_opt:n_g;
            else
                idx_seq = (i_opt - 1):-1:1;
            end
            x_prev = x_best;

            for j = idx_seq
                % Start 1: warm start from neighbouring grid solution
                xs1 = x_prev;  xs1(p) = grid_vals(j);
                [x1, f1] = refit(free_idx, xs1);

                % Start 2: from the global optimum
                xs2 = x_best;  xs2(p) = grid_vals(j);
                [x2, f2] = refit(free_idx, xs2);

                if isnan(f1) || (~isnan(f2) && f2 < f1)
                    x1 = x2;  f1 = f2;
                end
                prof_obj(j) = f1;
                prof_x(j, :) = x1;
                if ~isnan(f1), x_prev = x1; end

                % Previous method: one-at-a-time slice, no re-optimization
                xsl = x_cal;  xsl(p) = grid_vals(j);
                slice_obj(j) = objective_full(xsl);
            end
        end

        profile_data.(pname).param_value = grid_vals(:);
        profile_data.(pname).profile_obj = prof_obj;
        profile_data.(pname).refit_params = prof_x;
        profile_data.(pname).slice_obj = slice_obj;
        profile_data.(pname).i_opt = i_opt;
        profile_data.(pname).estimate = x_best(p);
        profile_data.(pname).estimate_published = x_cal(p);

        fprintf('done (%.0f s elapsed)\n', toc);
    end

    % Guard: if any profile point beat the reference, use it as the minimum
    all_prof = [];
    for p = 1:n_par
        all_prof = [all_prof; profile_data.(par_names{p}).profile_obj]; %#ok<AGROW>
    end
    f_ref = min([f_best; all_prof], [], 'omitnan');
    if f_ref < f_best - 1e-6
        fprintf(['\nNOTE: a profile point reached a lower objective (%.4f) than the\n' ...
                 'global re-fit (%.4f); deltas are measured from the lower value.\n'], ...
                 f_ref, f_best);
    end

    %% =========================================================================
    %% 5. CLASSIFY AND COMPUTE CONFIDENCE BOUNDS
    %% =========================================================================

    summary = table('Size', [n_par, 10], ...
        'VariableTypes', {'string', 'double', 'double', ...
                          'string', 'double', 'double', ...
                          'string', 'double', 'double', 'double'}, ...
        'VariableNames', {'parameter', 'estimate_published', 'estimate_refit', ...
                          'slice_status', 'slice_ci_lower', 'slice_ci_upper', ...
                          'profile_status', 'profile_ci_lower', 'profile_ci_upper', ...
                          'profile_max_delta'});

    for p = 1:n_par
        pname = par_names{p};
        d = profile_data.(pname);
        vals = d.param_value;

        % Re-optimized profile
        d.profile_delta = d.profile_obj - f_ref;
        [st_p, ci_p, flo_p, fhi_p] = classify_profile(vals, d.profile_delta, ...
            d.i_opt, chi2_threshold, lb(p), ub(p));

        % Previous slice (delta from slice minimum, as in the earlier code)
        d.slice_delta = d.slice_obj - min(d.slice_obj, [], 'omitnan');
        [~, i_cal] = min(abs(vals - x_cal(p)));
        [st_s, ci_s] = classify_profile(vals, d.slice_delta, ...
            i_cal, chi2_threshold, lb(p), ub(p));

        d.profile_status = st_p;
        d.profile_ci = ci_p;
        d.profile_found_lo = flo_p;
        d.profile_found_hi = fhi_p;
        d.slice_status = st_s;
        d.slice_ci = ci_s;
        profile_data.(pname) = d;

        summary.parameter(p) = pname;
        summary.estimate_published(p) = x_cal(p);
        summary.estimate_refit(p) = x_best(p);
        summary.slice_status(p) = st_s;
        summary.slice_ci_lower(p) = ci_s(1);
        summary.slice_ci_upper(p) = ci_s(2);
        summary.profile_status(p) = st_p;
        summary.profile_ci_lower(p) = ci_p(1);
        summary.profile_ci_upper(p) = ci_p(2);
        summary.profile_max_delta(p) = max(d.profile_delta, [], 'omitnan');
    end

    fprintf('\n=============================================================\n');
    fprintf('IDENTIFIABILITY: previous slice vs re-optimized profile\n');
    fprintf('=============================================================\n');
    fprintf('%-9s | %-26s | %-26s\n', 'Param', 'Slice (previous)', 'Profile (corrected)');
    fprintf('%s\n', repmat('-', 1, 68));
    for p = 1:n_par
        fprintf('%-9s | %-16s [%6.3g,%6.3g] | %-16s [%6.3g,%6.3g]\n', ...
            summary.parameter(p), ...
            summary.slice_status(p), summary.slice_ci_lower(p), summary.slice_ci_upper(p), ...
            summary.profile_status(p), summary.profile_ci_lower(p), summary.profile_ci_upper(p));
    end
    fprintf(['\nCI bounds equal to a grid edge mean the profile did not cross the\n' ...
             'threshold on that side within the evaluated range.\n']);
    fprintf(['Caveat: pain residuals are down-weighted (0.1 x ((pred-obs)/10)^2),\n' ...
             'as in calibration, so the 1.92 threshold is approximate.\n']);

    %% =========================================================================
    %% 6. SAVE RESULTS
    %% =========================================================================

    if ~exist('model_output', 'dir'), mkdir('model_output'); end

    profile_results.par_names = par_names;
    profile_results.profiles = profile_data;
    profile_results.summary = summary;
    profile_results.objective_published = obj_published;
    profile_results.objective_reference = f_ref;
    profile_results.params_refit = x_best;
    profile_results.chi2_threshold = chi2_threshold;
    profile_results.n_grid = n_grid;
    profile_results.method = ['Profile likelihood with re-optimization of the ' ...
        'remaining parameters (fmincon interior-point); slice shown for comparison.'];

    save('model_output/profile_likelihood_results.mat', 'profile_results');
    writetable(summary, 'model_output/profile_likelihood_summary.csv');

    fprintf('\nSaved: model_output/profile_likelihood_results.mat\n');
    fprintf('Saved: model_output/profile_likelihood_summary.csv\n');
    fprintf('Total time: %.1f min\n', toc / 60);

    %% =========================================================================
    %% NESTED HELPERS
    %% =========================================================================

    function [x_full, f] = refit(free_idx, x_start)
        % Optimize parameters in free_idx with the others fixed at x_start.
        x_full = x_start;
        f = NaN;

        lbf = lb(free_idx);
        ubf = ub(free_idx);
        A = [];
        b = [];

        syn_free = any(free_idx == i_syn);
        deg_free = any(free_idx == i_deg);

        if syn_free && deg_free
            A = zeros(1, numel(free_idx));
            A(free_idx == i_syn) = 1;
            A(free_idx == i_deg) = -1;
            b = -eps_k;                                   % k_syn - k_deg <= -eps
        elseif deg_free                                   % k_syn fixed
            k = find(free_idx == i_deg);
            lbf(k) = max(lbf(k), x_start(i_syn) + eps_k);
        elseif syn_free                                   % k_deg fixed
            k = find(free_idx == i_syn);
            ubf(k) = min(ubf(k), x_start(i_deg) - eps_k);
        end

        if any(lbf > ubf), return; end                    % infeasible grid point

        x0 = min(max(x_start(free_idx), lbf), ubf);
        if syn_free && deg_free
            ks = find(free_idx == i_syn);
            kd = find(free_idx == i_deg);
            if x0(kd) <= x0(ks) + eps_k
                x0(kd) = min(ubf(kd), x0(ks) + 2 * eps_k);
                if x0(kd) <= x0(ks) + eps_k
                    x0(ks) = max(lbf(ks), x0(kd) - 2 * eps_k);
                end
            end
        end

        fun = @(xf) objective_full(assemble(xf, free_idx, x_start));
        try
            [xf_opt, f_opt] = fmincon(fun, x0, A, b, [], [], lbf, ubf, [], opts);
            x_full(free_idx) = xf_opt;
            f = f_opt;
        catch
            % leave f = NaN
        end
    end

    function x_full = assemble(xf, free_idx, x_template)
        x_full = x_template;
        x_full(free_idx) = xf;
    end

    function sse = objective_full(x)
        % Same weighted objective as calibration (M04/M07); bounds and the
        % k_deg > k_syn constraint are enforced by the optimizer.
        params = calibrated_params;
        for ii = 1:n_par
            params.(par_names{ii}) = x(ii);
        end
        sse = compute_objective(params, jsw_targets, womac_targets);
    end

end

%% =============================================================================
%% LOCAL FUNCTIONS
%% =============================================================================

function [status, ci, found_lo, found_hi] = classify_profile(vals, delta, i_opt, thr, lb_p, ub_p)
    [lo, found_lo] = interp_crossing(vals, delta, i_opt, thr, -1);
    [hi, found_hi] = interp_crossing(vals, delta, i_opt, thr, +1);
    ci = [lo, hi];

    at_bound = abs(vals(i_opt) - lb_p) < 1e-9 || abs(vals(i_opt) - ub_p) < 1e-9;

    if found_lo && found_hi
        status = "identifiable";
    elseif found_lo || found_hi
        if at_bound
            status = "one-sided (bound)";
        else
            status = "one-sided";
        end
    else
        status = "not identifiable";
    end
end

function [bound, found] = interp_crossing(vals, delta, i_opt, thr, dirn)
    found = false;
    if dirn > 0
        seq = (i_opt + 1):numel(vals);
        bound = vals(end);
    else
        seq = (i_opt - 1):-1:1;
        bound = vals(1);
    end
    j_prev = i_opt;
    for j = seq
        if isnan(delta(j)), continue; end
        if delta(j) > thr
            d0 = delta(j_prev);
            d1 = delta(j);
            t = (thr - d0) / (d1 - d0);
            bound = vals(j_prev) + t * (vals(j) - vals(j_prev));
            found = true;
            return;
        end
        j_prev = j;
    end
end

function sse = compute_objective(params, jsw_targets, womac_targets)
    % Weighted sum of squared errors (identical to the calibration objective)

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
                if ~isempty(idx), pred = sim_p.JSW_change(idx); else, continue; end
            else
                idx = find(sim_t.week == wk);
                if ~isempty(idx), pred = sim_t.JSW_change(idx); else, continue; end
            end

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
                if ~isempty(idx), pred = sim_p24.Pain(idx); else, continue; end
            else
                idx = find(sim_t24.week == wk);
                if ~isempty(idx), pred = sim_t24.Pain(idx); else, continue; end
            end

            total_sse = total_sse + 0.1 * ((pred - obs) / 10)^2;
        end
    catch
        total_sse = total_sse + 50;
    end

    sse = total_sse;
end
