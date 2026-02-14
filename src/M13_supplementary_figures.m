%% =============================================================================
%% M13_supplementary_figures.m
%% Generate Supplementary Figures Matching R Manuscript
%% =============================================================================
%%
%% Figures:
%%   S1. Profile Likelihood Analysis for Parameter Identifiability
%%   S2. Residual Diagnostics for Model Calibration
%%   S3. Bootstrap Parameter Distributions
%%   S4. Bioavailability Impact on Efficacy
%%   S5. kpl Parameter Sensitivity
%%
%%
%% =============================================================================

fprintf('=============================================================\n');
fprintf('GENERATING SUPPLEMENTARY FIGURES\n');
fprintf('=============================================================\n');

% Resolve repository root (works whether called from run_all.m or directly)
this_script = mfilename('fullpath');
if ~isempty(this_script)
    src_dir = fileparts(this_script);
    repo_root_dir = fileparts(src_dir);
    cd(repo_root_dir);
end

if ~exist('figures', 'dir'), mkdir('figures'); end

% Load calibration data
load('model_output/calibration_results.mat', 'calibration_results');
params = calibration_results.params;
load('data_clean/calibration_data.mat', 'calibration_data');
jsw_data = calibration_data.jsw_data;
womac_data = calibration_data.womac_pain_data;

% Colors matching R ggplot2
col_blue = [0.00, 0.45, 0.70];
col_orange = [0.90, 0.45, 0.00];
col_green = [0.00, 0.60, 0.30];
col_gray = [0.5, 0.5, 0.5];
col_lightgray = [0.95, 0.95, 0.95];

%% =============================================================================
%% FIGURE S1: PROFILE LIKELIHOOD ANALYSIS
%% =============================================================================

fprintf('Generating Figure S1: Profile Likelihood Analysis...\n');

fig_s1 = figure('Position', [50, 50, 900, 1000], 'Color', 'w');

% Parameters to profile
param_names = {'k_syn', 'k_deg', 'Imax_deg', 'Pain_0', 'Pmax'};
param_labels = {'k\_syn (synthesis rate)', 'k\_deg (degradation rate)', ...
                'Imax\_deg (max inhibition)', 'Pain\_0 (baseline pain)', ...
                'Pmax (placebo effect)'};
param_values = [params.k_syn, params.k_deg, params.Imax_deg, params.Pain_0, params.Pmax];
param_ranges = {[0.08, 0.22], [0.10, 0.25], [0.5, 1.0], [35, 55], [10, 30]};

% 95% CI threshold (chi-square with 1 df)
chi2_threshold = 1.92;

for p = 1:5
    if p <= 4
        subplot(3, 2, p);
    else
        subplot(3, 2, 5);
    end
    hold on;
    set(gca, 'Color', col_lightgray);
    
    % Generate profile likelihood
    pname = param_names{p};
    pval_opt = param_values(p);
    prange = param_ranges{p};
    
    n_points = 20;
    pvals = linspace(prange(1), prange(2), n_points);
    delta_sse = zeros(size(pvals));
    
    % Compute SSE at each parameter value
    for i = 1:n_points
        test_params = params;
        test_params.(pname) = pvals(i);
        
        % Simulate and compute SSE
        try
            sim_p = simulate_glucosamine(test_params, 1500, 156, 'sulfate', false);
            sim_t = simulate_glucosamine(test_params, 1500, 156, 'sulfate', true);
            
            % JSW SSE
            sse_jsw = 0;
            for j = 1:height(jsw_data)
                week_idx = find(sim_p.week == jsw_data.week(j), 1);
                if isempty(week_idx), week_idx = round(jsw_data.week(j)) + 1; end
                week_idx = min(week_idx, length(sim_p.JSW_change));
                
                if contains(jsw_data.arm{j}, 'Placebo')
                    pred = sim_p.JSW_change(week_idx);
                else
                    pred = sim_t.JSW_change(week_idx);
                end
                sse_jsw = sse_jsw + (jsw_data.jsw_change_mm(j) - pred)^2;
            end
            
            % Pain SSE (simplified)
            sim_p_pain = simulate_glucosamine(test_params, 1500, 24, 'HCl', false);
            sim_t_pain = simulate_glucosamine(test_params, 1500, 24, 'HCl', true);
            
            sse_pain = 0;
            for j = 1:height(womac_data)
                week_idx = womac_data.week(j) + 1;
                week_idx = min(week_idx, length(sim_p_pain.Pain));
                
                if contains(womac_data.arm{j}, 'Placebo')
                    pred = sim_p_pain.Pain(week_idx);
                else
                    pred = sim_t_pain.Pain(week_idx);
                end
                sse_pain = sse_pain + (womac_data.womac_pain_normalized(j) - pred)^2;
            end
            
            delta_sse(i) = sse_jsw + sse_pain * 0.01;  % Weight pain less
        catch
            delta_sse(i) = NaN;
        end
    end
    
    % Normalize to minimum
    delta_sse = delta_sse - min(delta_sse);
    
    % Plot
    plot(pvals, delta_sse, '-o', 'Color', col_blue, 'LineWidth', 1.5, 'MarkerFaceColor', col_blue, 'MarkerSize', 4);
    
    % 95% CI threshold
    yline(chi2_threshold, '--', 'Color', col_orange, 'LineWidth', 1.5);
    text(prange(2)*0.98, chi2_threshold*1.1, '95% CI threshold', 'Color', col_orange, ...
         'HorizontalAlignment', 'right', 'FontSize', 8);
    
    % Optimal value
    xline(pval_opt, '-', 'Color', col_green, 'LineWidth', 2);
    
    % CI bounds (where delta_sse crosses threshold)
    ci_mask = delta_sse <= chi2_threshold;
    if any(ci_mask)
        ci_vals = pvals(ci_mask);
        xline(min(ci_vals), ':', 'Color', col_gray, 'LineWidth', 1);
        xline(max(ci_vals), ':', 'Color', col_gray, 'LineWidth', 1);
    end
    
    xlabel('Parameter Value', 'FontSize', 10);
    ylabel('\Delta SSE', 'FontSize', 10);
    title(param_labels{p}, 'FontSize', 11, 'FontWeight', 'bold');
    xlim(prange);
    box on;
    hold off;
end

sgtitle('Figure S1. Profile Likelihood Analysis for Parameter Identifiability', ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(fig_s1, 'figures/FigS1_profile_likelihood.png');
saveas(fig_s1, 'figures/FigS1_profile_likelihood.fig');
fprintf('  Saved: figures/FigS1_profile_likelihood.png\n');

%% =============================================================================
%% FIGURE S2: RESIDUAL DIAGNOSTICS
%% =============================================================================

fprintf('Generating Figure S2: Residual Diagnostics...\n');

fig_s2 = figure('Position', [50, 50, 1000, 800], 'Color', 'w');

% Simulate at optimal parameters
sim_p = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
sim_t = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
sim_p_pain = simulate_glucosamine(params, 1500, 24, 'HCl', false);
sim_t_pain = simulate_glucosamine(params, 1500, 24, 'HCl', true);

% Compute JSW residuals
jsw_pred = zeros(height(jsw_data), 1);
jsw_resid = zeros(height(jsw_data), 1);
jsw_time = zeros(height(jsw_data), 1);
jsw_arm = cell(height(jsw_data), 1);
jsw_study = cell(height(jsw_data), 1);

for j = 1:height(jsw_data)
    week_idx = round(jsw_data.week(j)) + 1;
    week_idx = min(max(week_idx, 1), length(sim_p.JSW_change));
    
    if contains(jsw_data.arm{j}, 'Placebo')
        jsw_pred(j) = sim_p.JSW_change(week_idx);
        jsw_arm{j} = 'Placebo';
    else
        jsw_pred(j) = sim_t.JSW_change(week_idx);
        jsw_arm{j} = 'Treatment';
    end
    jsw_resid(j) = jsw_data.jsw_change_mm(j) - jsw_pred(j);
    jsw_time(j) = jsw_data.week(j) / 52;
    jsw_study{j} = jsw_data.study{j};
end

% Compute Pain residuals
pain_pred = zeros(height(womac_data), 1);
pain_resid = zeros(height(womac_data), 1);
pain_arm = cell(height(womac_data), 1);

for j = 1:height(womac_data)
    week_idx = womac_data.week(j) + 1;
    week_idx = min(max(week_idx, 1), length(sim_p_pain.Pain));
    
    if contains(womac_data.arm{j}, 'Placebo')
        pain_pred(j) = sim_p_pain.Pain(week_idx);
        pain_arm{j} = 'Placebo';
    else
        pain_pred(j) = sim_t_pain.Pain(week_idx);
        pain_arm{j} = 'Treatment';
    end
    pain_resid(j) = womac_data.womac_pain_normalized(j) - pain_pred(j);
end

% Panel A: JSW Residuals vs Predicted
subplot(2, 2, 1);
hold on;
set(gca, 'Color', col_lightgray);

for j = 1:length(jsw_pred)
    if strcmp(jsw_arm{j}, 'Placebo')
        col = col_gray;
    else
        col = col_green;
    end
    if contains(jsw_study{j}, 'Pavelka')
        marker = 'o';
    else
        marker = '^';
    end
    scatter(jsw_pred(j), jsw_resid(j), 80, col, 'filled', marker);
end

yline(0, '--', 'Color', col_gray, 'LineWidth', 1);

% Add LOESS-like curve (polynomial fit)
[sorted_pred, sort_idx] = sort(jsw_pred);
sorted_resid = jsw_resid(sort_idx);
p_fit = polyfit(sorted_pred, sorted_resid, 3);
x_fit = linspace(min(jsw_pred), max(jsw_pred), 50);
y_fit = polyval(p_fit, x_fit);
plot(x_fit, y_fit, ':', 'Color', col_gray, 'LineWidth', 1.5);

xlabel('Predicted JSW Change (mm)', 'FontSize', 10);
ylabel('Residual (mm)', 'FontSize', 10);
title('A. Residuals vs Predicted (JSW)', 'FontSize', 11, 'FontWeight', 'bold');
legend({'Placebo', 'Treatment'}, 'Location', 'best', 'FontSize', 8);
box on; hold off;

% Panel B: JSW Residuals vs Time
subplot(2, 2, 2);
hold on;
set(gca, 'Color', col_lightgray);

for j = 1:length(jsw_time)
    if strcmp(jsw_arm{j}, 'Placebo')
        col = col_gray;
    else
        col = col_green;
    end
    if contains(jsw_study{j}, 'Pavelka')
        marker = 'o';
    else
        marker = '^';
    end
    scatter(jsw_time(j), jsw_resid(j), 80, col, 'filled', marker);
end

yline(0, '--', 'Color', col_gray, 'LineWidth', 1);

xlabel('Time (years)', 'FontSize', 10);
ylabel('Residual (mm)', 'FontSize', 10);
title('B. Residuals vs Time (JSW)', 'FontSize', 11, 'FontWeight', 'bold');
box on; hold off;

% Panel C: Pain Residuals vs Predicted
subplot(2, 2, 3);
hold on;
set(gca, 'Color', col_lightgray);

for j = 1:length(pain_pred)
    if strcmp(pain_arm{j}, 'Placebo')
        col = col_gray;
    else
        col = col_green;
    end
    scatter(pain_pred(j), pain_resid(j), 80, col, 'filled', 'o');
end

yline(0, '--', 'Color', col_gray, 'LineWidth', 1);

xlabel('Predicted Pain Score', 'FontSize', 10);
ylabel('Residual (points)', 'FontSize', 10);
title('C. Residuals vs Predicted (WOMAC)', 'FontSize', 11, 'FontWeight', 'bold');
box on; hold off;

% Panel D: Weighted Residuals Distribution
subplot(2, 2, 4);
hold on;
set(gca, 'Color', col_lightgray);

% Weighted residuals (standardized)
jsw_se = 0.05;  % Approximate SE for JSW
pain_se = 3;    % Approximate SE for pain
weighted_jsw = jsw_resid / jsw_se;
weighted_pain = pain_resid / pain_se;

% Histograms
edges = -3:0.5:3;
histogram(weighted_jsw, edges, 'FaceColor', col_blue, 'FaceAlpha', 0.6, 'EdgeColor', 'w', 'Normalization', 'pdf');
histogram(weighted_pain, edges, 'FaceColor', col_orange, 'FaceAlpha', 0.6, 'EdgeColor', 'w', 'Normalization', 'pdf');

% Standard normal overlay
x_norm = linspace(-3, 3, 100);
y_norm = normpdf(x_norm, 0, 1);
plot(x_norm, y_norm, 'k:', 'LineWidth', 2);

xline(0, '--', 'Color', col_gray, 'LineWidth', 1);

xlabel('Weighted Residual', 'FontSize', 10);
ylabel('Density', 'FontSize', 10);
title('D. Weighted Residuals Distribution', 'FontSize', 11, 'FontWeight', 'bold');
text(0.02, 0.98, 'Dotted line: Standard normal N(0,1)', 'Units', 'normalized', ...
     'VerticalAlignment', 'top', 'FontSize', 8, 'Color', col_gray);
legend({'JSW', 'WOMAC'}, 'Location', 'northeast', 'FontSize', 8);
box on; hold off;

sgtitle('Figure S2. Residual Diagnostics for Model Calibration', ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(fig_s2, 'figures/FigS2_residual_diagnostics.png');
saveas(fig_s2, 'figures/FigS2_residual_diagnostics.fig');
fprintf('  Saved: figures/FigS2_residual_diagnostics.png\n');

%% =============================================================================
%% FIGURE S3: BOOTSTRAP PARAMETER DISTRIBUTIONS
%% =============================================================================
%% CORRECTED: Now loads ACTUAL bootstrap results from M07 output
%% instead of generating fake lognormal approximations.
%%
%% The real bootstrap distributions will show:
%%   - Imax_deg: boundary pileup near 1.0 (NOT a smooth bell curve)
%%   - k_syn, k_deg: possible bimodality from structural correlation
%%   - Pain_0, Pmax: may be tight but should come from real resampling
%%
%% If M07 has not been run, this section will fall back to a warning.
%% =============================================================================

fprintf('Generating Figure S3: Bootstrap Parameter Distributions...\n');

fig_s3 = figure('Position', [50, 50, 900, 1000], 'Color', 'w');

param_names = {'k_syn', 'k_deg', 'Imax_deg', 'Pain_0', 'Pmax'};
param_labels = {'k\_syn', 'k\_deg', 'Imax\_deg', 'Pain\_0', 'Pmax'};

% --- LOAD ACTUAL BOOTSTRAP RESULTS FROM M07 ---
try
    load('model_output/bootstrap_results.mat', 'bootstrap_results');
    
    boot_data = bootstrap_results.boot_params;  % N x 5 matrix from M07
    n_boot = size(boot_data, 1);
    par_names_m07 = bootstrap_results.par_names;
    
    fprintf('  Loaded %d actual bootstrap replicates from M07\n', n_boot);
    
    % Map M07 column order to our parameter order
    boot_params = cell(5, 1);
    for p = 1:5
        col_idx = find(strcmp(par_names_m07, param_names{p}));
        if ~isempty(col_idx)
            boot_params{p} = boot_data(:, col_idx);
        else
            error('Parameter %s not found in M07 output', param_names{p});
        end
    end
    
    boot_estimates = [params.k_syn, params.k_deg, params.Imax_deg, params.Pain_0, params.Pmax];
    
    for p = 1:5
        if p <= 4
            subplot(3, 2, p);
        else
            subplot(3, 2, 5);
        end
        hold on;
        set(gca, 'Color', col_lightgray);
        
        data = boot_params{p};
        est = boot_estimates(p);
        
        % Histogram
        histogram(data, 30, 'FaceColor', col_blue, 'EdgeColor', 'w', 'FaceAlpha', 0.7, 'Normalization', 'pdf');
        
        % Kernel density
        [f, xi] = ksdensity(data);
        plot(xi, f, '-', 'Color', col_orange, 'LineWidth', 2);
        
        % Point estimate
        xline(est, '-', 'Color', col_green, 'LineWidth', 2);
        
        % 95% CI bounds
        ci_lo = prctile(data, 2.5);
        ci_hi = prctile(data, 97.5);
        xline(ci_lo, '--', 'Color', col_gray, 'LineWidth', 1.5);
        xline(ci_hi, '--', 'Color', col_gray, 'LineWidth', 1.5);
        
        % Statistics annotation
        cv_pct = std(data) / mean(data) * 100;
        pct_at_bound = 0;
        if p == 3  % Imax_deg
            pct_at_bound = sum(data >= 0.95) / length(data) * 100;
            text(0.02, 0.92, sprintf('%.1f%% at bound', pct_at_bound), ...
                 'Units', 'normalized', 'FontSize', 8, 'Color', 'r');
        end
        text(0.02, 0.82, sprintf('CV = %.1f%%', cv_pct), ...
             'Units', 'normalized', 'FontSize', 8, 'Color', col_gray);
        text(0.02, 0.72, sprintf('95%% CI: [%.4f, %.4f]', ci_lo, ci_hi), ...
             'Units', 'normalized', 'FontSize', 7, 'Color', col_gray);
        
        xlabel('Parameter Value', 'FontSize', 10);
        ylabel('Density', 'FontSize', 10);
        title(param_labels{p}, 'FontSize', 11, 'FontWeight', 'bold');
        box on;
        hold off;
    end
    
    % Add correlation info in panel 6
    subplot(3, 2, 6);
    set(gca, 'Color', col_lightgray);
    
    cor_matrix = bootstrap_results.correlation_matrix;
    imagesc(cor_matrix);
    colorbar;
    caxis([-1, 1]);
    colormap(gca, parula);
    set(gca, 'XTick', 1:5, 'XTickLabel', param_labels, 'YTick', 1:5, 'YTickLabel', param_labels);
    title('Parameter Correlations', 'FontSize', 11, 'FontWeight', 'bold');
    
    sgtitle(sprintf('Figure S3. Bootstrap Parameter Distributions (N = %d actual replicates)', n_boot), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
catch ME
    % If M07 results not available, display warning
    warning('Could not load M07 bootstrap results: %s', ME.message);
    
    subplot(1, 1, 1);
    text(0.5, 0.5, {'\bfFigure S3 cannot be generated'; ''; ...
         'M07 bootstrap results not found.'; ...
         'You must run M07_bootstrap_ci.m first.'; ''; ...
         'DO NOT use fake/synthetic bootstrap distributions.'}, ...
         'Units', 'normalized', 'HorizontalAlignment', 'center', ...
         'FontSize', 14, 'Color', 'r');
    axis off;
    
    sgtitle('Figure S3. Bootstrap Parameter Distributions (NOT YET GENERATED)', ...
            'FontSize', 14, 'FontWeight', 'bold');
end

saveas(fig_s3, 'figures/FigS3_bootstrap_distributions.png');
saveas(fig_s3, 'figures/FigS3_bootstrap_distributions.fig');
fprintf('  Saved: figures/FigS3_bootstrap_distributions.png\n');

%% =============================================================================
%% FIGURE S4: BIOAVAILABILITY IMPACT ON EFFICACY - FIXED TITLE OVERLAP
%% =============================================================================

fprintf('Generating Figure S4: Bioavailability Impact...\n');

fig_s4 = figure('Position', [100, 100, 700, 500], 'Color', 'w');
hold on;
set(gca, 'Color', col_lightgray);

% Bioavailability range
F_values = [10, 15, 20, 25, 30, 35, 40];
te_values = zeros(size(F_values));

for i = 1:length(F_values)
    test_params = params;
    test_params.F_sulfate = F_values(i) / 100;
    
    sim_p = simulate_glucosamine(test_params, 1500, 156, 'sulfate', false);
    sim_t = simulate_glucosamine(test_params, 1500, 156, 'sulfate', true);
    
    te_values(i) = sim_t.JSW_change(end) - sim_p.JSW_change(end);
end

% Fill area under curve above MCID
fill([F_values, fliplr(F_values)], [te_values, 0.1*ones(size(te_values))], ...
     col_blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot line with markers
plot(F_values, te_values, '-o', 'Color', col_blue, 'LineWidth', 2, ...
     'MarkerFaceColor', col_orange, 'MarkerEdgeColor', col_orange, 'MarkerSize', 8);

% MCID line
yline(0.1, ':', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.5);
text(38, 0.1, 'MCID = 0.1 mm', 'Color', [0, 0.6, 0.5], 'FontSize', 10, ...
     'VerticalAlignment', 'bottom');

% Current F line
xline(22, '--', 'Color', col_gray, 'LineWidth', 1.5);
text(22.5, 0.25, 'Current', 'Color', col_gray, 'FontSize', 10);
text(22.5, 0.23, 'F = 22%', 'Color', col_gray, 'FontSize', 10);

xlabel('Bioavailability (%)', 'FontSize', 12);
ylabel('3-Year JSW Treatment Effect (mm)', 'FontSize', 12);
% FIX: Two-line title to avoid overlap with subtitle
title({'Bioavailability Impact on Efficacy'; ...
       'Potential for formulation improvement to enhance structure modification'}, ...
      'FontSize', 14, 'FontWeight', 'bold');
% FIX: Moved note inside plot area to avoid overlap with xlabel
text(0.98, 0.06, 'MCID = Minimum Clinically Important Difference', ...
     'Units', 'normalized', 'FontSize', 8, 'Color', col_gray, ...
     'HorizontalAlignment', 'right');

xlim([8, 42]);
ylim([0, 0.4]);
box on;
hold off;

saveas(fig_s4, 'figures/FigS4_bioavailability_impact.png');
saveas(fig_s4, 'figures/FigS4_bioavailability_impact.fig');
fprintf('  Saved: figures/FigS4_bioavailability_impact.png\n');

%% =============================================================================
%% FIGURE S5: KPL PARAMETER SENSITIVITY - FIXED TITLE OVERLAP
%% =============================================================================

fprintf('Generating Figure S5: kpl Parameter Sensitivity...\n');

fig_s5 = figure('Position', [100, 100, 700, 500], 'Color', 'w');
hold on;

% kpl values on log scale
kpl_values = [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.5, 1, 10];
pain_values = zeros(size(kpl_values));

for i = 1:length(kpl_values)
    test_params = params;
    test_params.kpl = kpl_values(i);
    
    sim_p = simulate_glucosamine(test_params, 1500, 24, 'HCl', false);
    pain_values(i) = sim_p.Pain(end);
end

% Background shading
fill([0.001, 0.05, 0.05, 0.001], [25, 25, 40, 40], ...
     [0.8, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
fill([0.05, 10, 10, 0.05], [25, 25, 40, 40], ...
     [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot
plot(kpl_values, pain_values, '-o', 'Color', col_blue, 'LineWidth', 2, ...
     'MarkerFaceColor', col_orange, 'MarkerEdgeColor', col_orange, 'MarkerSize', 8);

% Observed line
yline(30.2, '--', 'Color', col_gray, 'LineWidth', 1.5);
text(0.3, 30.5, 'Observed = 30.2', 'Color', col_gray, 'FontSize', 10);

% Labels
text(0.008, 26, 'Physiological', 'Color', [0, 0.5, 0], 'FontSize', 9, 'FontAngle', 'italic');
text(0.008, 25.3, 'range', 'Color', [0, 0.5, 0], 'FontSize', 9, 'FontAngle', 'italic');

set(gca, 'XScale', 'log', 'Color', 'none');
xlabel('k_{pl} (day^{-1}, log scale)', 'FontSize', 12);
ylabel('Week 24 Pain Score (0-100)', 'FontSize', 12);
% FIX: Two-line title to avoid overlap with subtitle
title({'kpl Parameter Sensitivity'; ...
       'Effect of placebo rate constant on Week 24 pain score'}, ...
      'FontSize', 14, 'FontWeight', 'bold');
% FIX: Moved note inside plot area (top-right) to avoid overlap with xlabel
text(0.98, 0.95, 'Note: Pain plateaus at high k_{pl} (''sloppy'' beyond ~0.05)', ...
     'Units', 'normalized', 'FontSize', 8, 'Color', col_gray, 'FontAngle', 'italic', ...
     'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

xlim([0.001, 10]);
ylim([28, 40]);
box on;
hold off;

saveas(fig_s5, 'figures/FigS5_kpl_sensitivity.png');
saveas(fig_s5, 'figures/FigS5_kpl_sensitivity.fig');
fprintf('  Saved: figures/FigS5_kpl_sensitivity.png\n');

%% =============================================================================
%% COMPLETION
%% =============================================================================

fprintf('=============================================================\n');
fprintf('ALL SUPPLEMENTARY FIGURES COMPLETE\n');
fprintf('=============================================================\n');
fprintf('Saved to: figures/\n');
fprintf('  - FigS1_profile_likelihood.png\n');
fprintf('  - FigS2_residual_diagnostics.png\n');
fprintf('  - FigS3_bootstrap_distributions.png\n');
fprintf('  - FigS4_bioavailability_impact.png\n');
fprintf('  - FigS5_kpl_sensitivity.png\n');
fprintf('=============================================================\n');
