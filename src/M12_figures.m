%% =============================================================================
%% M12_figures.m
%% Generate Publication-Ready Figures Matching R Manuscript
%% =============================================================================
%%
%%
%% =============================================================================

fprintf('=============================================================\n');
fprintf('GENERATING PUBLICATION FIGURES\n');
fprintf('=============================================================\n');

% Resolve repository root (works whether called from run_all.m or directly)
this_script = mfilename('fullpath');
if ~isempty(this_script)
    src_dir = fileparts(this_script);
    repo_root_dir = fileparts(src_dir);
    cd(repo_root_dir);
end

if ~exist('figures', 'dir'), mkdir('figures'); end

% Load data
load('model_output/calibration_results.mat', 'calibration_results');
params = calibration_results.params;
load('data_clean/calibration_data.mat', 'calibration_data');
jsw_data = calibration_data.jsw_data;
womac_data = calibration_data.womac_pain_data;

% R-style colors
col_treat = [0.00, 0.45, 0.70];   % Blue
col_plac = [0.90, 0.45, 0.00];    % Orange
col_gray = [0.5, 0.5, 0.5];

% Font sizes
fs_title = 14;
fs_axis = 11;
fs_legend = 9;

% Parameter display labels (proper TeX subscripts for figures)
param_display_map = containers.Map({...
    'k_syn', 'k_deg', 'GAG_0', 'JSW_0', 'F_sulfate', 'F_HCl', ...
    'Imax_deg', 'IC50_deg', 'Emax_syn', 'EC50_syn', 'beta_drug', ...
    'alpha_struct', 'kpl', 'Pmax', 'Pain_0', 'EC50_pain'}, {...
    'k_{syn}', 'k_{deg}', 'GAG_0', 'JSW_0', 'F_{sulfate}', 'F_{HCl}', ...
    'I_{max,deg}', 'IC_{50,deg}', 'E_{max,syn}', 'EC_{50,syn}', '\beta_{drug}', ...
    '\alpha_{struct}', 'k_{pl}', 'P_{max}', 'Pain_0', 'EC_{50,pain}'});

%% =============================================================================
%% FIGURE 2: MODEL CALIBRATION (3 vertical panels)
%% =============================================================================

fprintf('Generating Figure 2: Model Calibration...\n');

sim_p = simulate_glucosamine(params, 1500, 156, 'sulfate', false);
sim_t = simulate_glucosamine(params, 1500, 156, 'sulfate', true);
sim_p_pain = simulate_glucosamine(params, 1500, 24, 'HCl', false);
sim_t_pain = simulate_glucosamine(params, 1500, 24, 'HCl', true);

fig1 = figure('Position', [100, 50, 700, 900], 'Color', 'w');

% --- Panel A: JSW ---
subplot(3, 1, 1);
hold on;
set(gca, 'Color', [0.95, 0.95, 0.95]);

plot(sim_t.time_years, sim_t.JSW_change, 'Color', col_treat, 'LineWidth', 2);
plot(sim_p.time_years, sim_p.JSW_change, 'Color', col_plac, 'LineWidth', 2);

% Reginster 2001 (triangles)
reg_p = jsw_data(strcmp(jsw_data.study, 'Reginster 2001') & contains(jsw_data.arm, 'Placebo'), :);
reg_t = jsw_data(strcmp(jsw_data.study, 'Reginster 2001') & contains(jsw_data.arm, 'GS'), :);
scatter(reg_t.week/52, reg_t.jsw_change_mm, 70, col_treat, 'filled', '^');
scatter(reg_p.week/52, reg_p.jsw_change_mm, 70, col_plac, 'filled', '^');

% Pavelka 2002 (circles)
pav_p = jsw_data(strcmp(jsw_data.study, 'Pavelka 2002') & contains(jsw_data.arm, 'Placebo'), :);
pav_t = jsw_data(strcmp(jsw_data.study, 'Pavelka 2002') & contains(jsw_data.arm, 'GS'), :);
scatter(pav_t.week/52, pav_t.jsw_change_mm, 70, col_treat, 'filled', 'o');
scatter(pav_p.week/52, pav_p.jsw_change_mm, 70, col_plac, 'filled', 'o');

yline(0, '--', 'Color', col_gray);
xlabel('Time (years)', 'FontSize', fs_axis);
ylabel('JSW Change from Baseline (mm)', 'FontSize', fs_axis);
title('A. Joint Space Width Change', 'FontSize', fs_title, 'FontWeight', 'bold');
text(0.02, 0.98, 'Calibration: Reginster 2001 & Pavelka 2002 (crystalline glucosamine sulfate)', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', col_gray);
xlim([0, 3.2]); ylim([-0.4, 0.2]);
legend({'GS 1500mg', 'Placebo'}, 'Location', 'southwest', 'FontSize', fs_legend);
box on; hold off;

% --- Panel B: GAG Content ---
subplot(3, 1, 2);
hold on;
set(gca, 'Color', [0.95, 0.95, 0.95]);

plot(sim_t.time_years, sim_t.GAG, 'Color', col_treat, 'LineWidth', 2);
plot(sim_p.time_years, sim_p.GAG, 'Color', col_plac, 'LineWidth', 2);
yline(1.0, ':', 'Color', col_gray, 'LineWidth', 1.5);
text(3.1, 1.0, 'Healthy', 'FontSize', 8, 'Color', col_gray);

xlabel('Time (years)', 'FontSize', fs_axis);
ylabel('GAG Content (normalized)', 'FontSize', fs_axis);
title('B. Cartilage GAG Content', 'FontSize', fs_title, 'FontWeight', 'bold');
text(0.02, 0.98, 'Model prediction (normalized to healthy cartilage = 1.0)', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', col_gray);
xlim([0, 3.2]); ylim([0.88, 1.02]);
legend({'GS 1500mg', 'Placebo'}, 'Location', 'southwest', 'FontSize', fs_legend);
box on; hold off;

% --- Panel C: Pain ---
subplot(3, 1, 3);
hold on;
set(gca, 'Color', [0.95, 0.95, 0.95]);

plot(sim_t_pain.week, sim_t_pain.Pain, 'Color', col_treat, 'LineWidth', 2);
plot(sim_p_pain.week, sim_p_pain.Pain, 'Color', col_plac, 'LineWidth', 2);

gait_p = womac_data(contains(womac_data.arm, 'Placebo'), :);
gait_t = womac_data(contains(womac_data.arm, 'GH'), :);
scatter(gait_t.week, gait_t.womac_pain_normalized, 70, col_treat, 'filled', 'o');
scatter(gait_p.week, gait_p.womac_pain_normalized, 70, col_plac, 'filled', 'o');

xlabel('Time (weeks)', 'FontSize', fs_axis);
ylabel('WOMAC Pain (0-100)', 'FontSize', fs_axis);
title('C. WOMAC Pain Score', 'FontSize', fs_title, 'FontWeight', 'bold');
text(0.02, 0.98, 'Calibration: GAIT Trial 2006 (glucosamine hydrochloride formulation)', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', col_gray);
xlim([0, 26]); ylim([25, 55]);
legend({'GH 1500mg (HCl)', 'Placebo'}, 'Location', 'northeast', 'FontSize', fs_legend);
box on; hold off;

saveas(fig1, 'figures/Fig2_calibration.png');
saveas(fig1, 'figures/Fig2_calibration.fig');
fprintf('  Saved: figures/Fig2_calibration.png\n');
fprintf('  Treatment Effect: %.3f mm\n', sim_t.JSW_change(end) - sim_p.JSW_change(end));

%% =============================================================================
%% FIGURE 3: SENSITIVITY ANALYSIS (3 tornado plots)
%% =============================================================================

fprintf('Generating Figure 3: Sensitivity Analysis...\n');

try
    load('model_output/sensitivity_analysis_results.mat', 'sensitivity_results');
    sens = sensitivity_results.sensitivity_summary;
    baseline = sensitivity_results.baseline;
    
    % FIX: Larger figure, show top 10 params (not 15) to avoid cramping
    fig2 = figure('Position', [100, 50, 900, 1050], 'Color', 'w');
    
    n_show = min(10, height(sens));
    
    % --- Panel A: Treatment Effect ---
    subplot(3, 1, 1);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    [~, idx] = sort(sens.sens_treatment_effect, 'descend');
    top_idx = idx(1:n_show);
    sens_a = sens(top_idx, :);
    sens_a = flipud(sens_a);
    
    for i = 1:height(sens_a)
        lo_val = sens_a.te_down(i);
        hi_val = sens_a.te_up(i);
        lo = min(lo_val, hi_val);
        hi = max(lo_val, hi_val);
        barh(i, hi - lo, 'BaseValue', lo, 'FaceColor', col_treat, 'EdgeColor', 'none');
        plot(lo_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_plac, 'MarkerEdgeColor', 'k');
        plot(hi_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_treat, 'MarkerEdgeColor', 'k');
    end
    xline(baseline.te, '--', 'Color', col_gray, 'LineWidth', 1.5);
    labels_a = cellfun(@(x) param_display_map(x), sens_a.parameter, 'UniformOutput', false);
    set(gca, 'YTick', 1:height(sens_a), 'YTickLabel', labels_a, 'FontSize', 9, 'TickLabelInterpreter', 'tex');
    ylim([0.5, n_show + 0.5]);
    xlabel('JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    title({'A. Sensitivity: 3-Year JSW Treatment Effect'; ...
           sprintf('Baseline: %.3f mm | \\pm20%% parameter perturbation', baseline.te)}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    box on; hold off;
    
    % --- Panel B: Placebo JSW (NEGATIVE VALUES) ---
    subplot(3, 1, 2);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    [~, idx] = sort(sens.sens_placebo, 'descend');
    top_idx = idx(1:n_show);
    sens_b = sens(top_idx, :);
    sens_b = flipud(sens_b);
    
    for i = 1:height(sens_b)
        lo_val = sens_b.plac_down(i);
        hi_val = sens_b.plac_up(i);
        lo = min(lo_val, hi_val);
        hi = max(lo_val, hi_val);
        barh(i, hi - lo, 'BaseValue', lo, 'FaceColor', col_treat, 'EdgeColor', 'none');
        plot(lo_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_plac, 'MarkerEdgeColor', 'k');
        plot(hi_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_treat, 'MarkerEdgeColor', 'k');
    end
    xline(baseline.jsw_placebo, '--', 'Color', col_gray, 'LineWidth', 1.5);
    labels_b = cellfun(@(x) param_display_map(x), sens_b.parameter, 'UniformOutput', false);
    set(gca, 'YTick', 1:height(sens_b), 'YTickLabel', labels_b, 'FontSize', 9, 'TickLabelInterpreter', 'tex');
    ylim([0.5, n_show + 0.5]);
    xlabel('JSW Y3 Placebo (mm)', 'FontSize', fs_axis);
    title({'B. Sensitivity: 3-Year Placebo JSW Change'; ...
           sprintf('Baseline: %.3f mm | \\pm20%% parameter perturbation', baseline.jsw_placebo)}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    box on; hold off;
    
    % --- Panel C: Week 24 Pain ---
    subplot(3, 1, 3);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    [~, idx] = sort(sens.sens_pain, 'descend');
    top_idx = idx(1:n_show);
    sens_c = sens(top_idx, :);
    sens_c = flipud(sens_c);
    
    for i = 1:height(sens_c)
        lo_val = sens_c.pain_down(i);
        hi_val = sens_c.pain_up(i);
        lo = min(lo_val, hi_val);
        hi = max(lo_val, hi_val);
        barh(i, hi - lo, 'BaseValue', lo, 'FaceColor', col_treat, 'EdgeColor', 'none');
        plot(lo_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_plac, 'MarkerEdgeColor', 'k');
        plot(hi_val, i, 'o', 'MarkerSize', 6, 'MarkerFaceColor', col_treat, 'MarkerEdgeColor', 'k');
    end
    xline(baseline.pain_24w, '--', 'Color', col_gray, 'LineWidth', 1.5);
    labels_c = cellfun(@(x) param_display_map(x), sens_c.parameter, 'UniformOutput', false);
    set(gca, 'YTick', 1:height(sens_c), 'YTickLabel', labels_c, 'FontSize', 9, 'TickLabelInterpreter', 'tex');
    ylim([0.5, n_show + 0.5]);
    xlabel('Pain W24 Placebo (points)', 'FontSize', fs_axis);
    title({'C. Sensitivity: Week 24 Placebo Pain'; ...
           sprintf('Baseline: %.1f | \\pm20%% parameter perturbation', baseline.pain_24w)}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    box on; hold off;
    
    saveas(fig2, 'figures/Fig3_sensitivity.png');
    saveas(fig2, 'figures/Fig3_sensitivity.fig');
    fprintf('  Saved: figures/Fig3_sensitivity.png\n');
    fprintf('  Baseline TE: %.3f mm, Placebo: %.3f mm, Pain: %.1f\n', ...
            baseline.te, baseline.jsw_placebo, baseline.pain_24w);
catch ME
    warning('Sensitivity figure error: %s', ME.message);
end

%% =============================================================================
%% FIGURE 4: VIRTUAL POPULATION (4 panels)
%% =============================================================================

fprintf('Generating Figure 4: Virtual Population...\n');

try
    load('model_output/virtual_population_results.mat', 'vpop_results');
    
    % FIX: Taller figure for more spacing between panels
    fig3 = figure('Position', [100, 50, 1050, 920], 'Color', 'w');
    
    % Get the data
    jsw_p = vpop_results.jsw_summary_placebo;
    jsw_t = vpop_results.jsw_summary_treatment;
    gag_p = vpop_results.gag_summary_placebo;
    gag_t = vpop_results.gag_summary_treatment;
    
    time_y = jsw_p.week / 52;
    
    % Panel A: JSW with 90% PI
    subplot(2, 2, 1);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    fill([time_y; flipud(time_y)], [jsw_t.p5; flipud(jsw_t.p95)], col_treat, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([time_y; flipud(time_y)], [jsw_p.p5; flipud(jsw_p.p95)], col_plac, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(time_y, jsw_t.median, 'Color', col_treat, 'LineWidth', 2);
    plot(time_y, jsw_p.median, 'Color', col_plac, 'LineWidth', 2);
    
    scatter(reg_t.week/52, reg_t.jsw_change_mm, 50, col_treat, 'filled', '^');
    scatter(reg_p.week/52, reg_p.jsw_change_mm, 50, col_plac, 'filled', '^');
    scatter(pav_t.week/52, pav_t.jsw_change_mm, 50, col_treat, 'filled', 'o');
    scatter(pav_p.week/52, pav_p.jsw_change_mm, 50, col_plac, 'filled', 'o');
    
    yline(0, '--', 'Color', col_gray);
    xlabel('Time (years)', 'FontSize', fs_axis);
    ylabel('JSW Change from Baseline (mm)', 'FontSize', fs_axis);
    title({'A. Joint Space Width Change Over 3 Years'; ...
           'Shaded: 90% PI | Lines: median'}, 'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GS 1500mg', 'Placebo'}, 'Location', 'southwest', 'FontSize', fs_legend);
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_min_a = min([jsw_p.p5; jsw_t.p5]);
    y_max_a = max([jsw_p.p95; jsw_t.p95; reg_t.jsw_change_mm; pav_t.jsw_change_mm]);
    y_pad_a = 0.10 * (y_max_a - y_min_a);
    xlim([0, 3.2]); ylim([y_min_a - y_pad_a, y_max_a + y_pad_a]); box on; hold off;
    
    % Panel B: GAG with 90% PI
    subplot(2, 2, 2);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    fill([time_y; flipud(time_y)], [gag_t.p5*100; flipud(gag_t.p95*100)], col_treat, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    fill([time_y; flipud(time_y)], [gag_p.p5*100; flipud(gag_p.p95*100)], col_plac, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(time_y, gag_t.median*100, 'Color', col_treat, 'LineWidth', 2);
    plot(time_y, gag_p.median*100, 'Color', col_plac, 'LineWidth', 2);
    yline(100, ':', 'Color', col_gray, 'LineWidth', 1.5);
    text(3.1, 100, 'Healthy', 'FontSize', 8, 'Color', col_gray);
    
    xlabel('Time (years)', 'FontSize', fs_axis);
    ylabel('GAG Content (% of healthy)', 'FontSize', fs_axis);
    title({'B. Cartilage GAG Content'; ...
           'Shaded: 90% PI | Dotted: healthy reference'}, 'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GS 1500mg', 'Placebo'}, 'Location', 'southwest', 'FontSize', fs_legend);
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_min_b = min([gag_p.p5*100; gag_t.p5*100]);
    y_max_b = max([gag_p.p95*100; gag_t.p95*100]);
    y_pad_b = 0.10 * (y_max_b - y_min_b);
    xlim([0, 3.2]); ylim([y_min_b - y_pad_b, y_max_b + y_pad_b]); box on; hold off;
    
    % Panel C: Pain with 90% PI from ACTUAL VPop simulations
    % CORRECTED: Uses real simulation output from M09, not hardcoded linspace
    subplot(2, 2, 3);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    % Load actual pain trajectories from VPop results
    pain_p_summary = vpop_results.pain_summary_placebo;
    pain_t_summary = vpop_results.pain_summary_treatment;
    weeks_pain = pain_p_summary.week;
    
    pain_t_mean = pain_t_summary.median;
    pain_p_mean = pain_p_summary.median;
    pain_t_lo = pain_t_summary.p5;
    pain_t_hi = pain_t_summary.p95;
    pain_p_lo = pain_p_summary.p5;
    pain_p_hi = pain_p_summary.p95;
    
    fill([weeks_pain; flipud(weeks_pain)], [pain_t_lo; flipud(pain_t_hi)], col_treat, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([weeks_pain; flipud(weeks_pain)], [pain_p_lo; flipud(pain_p_hi)], col_plac, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(weeks_pain, pain_t_mean, 'Color', col_treat, 'LineWidth', 2);
    plot(weeks_pain, pain_p_mean, 'Color', col_plac, 'LineWidth', 2);
    
    scatter(gait_t.week, gait_t.womac_pain_normalized, 50, col_treat, 'filled', 'o');
    scatter(gait_p.week, gait_p.womac_pain_normalized, 50, col_plac, 'filled', 'o');
    
    xlabel('Time (weeks)', 'FontSize', fs_axis);
    ylabel('WOMAC Pain (0-100)', 'FontSize', fs_axis);
    title({'C. WOMAC Pain Score'; ...
           'GAIT trial (GH formulation) | 90% PI'}, 'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GH 1500mg (HCl)', 'Placebo'}, 'Location', 'northeast', 'FontSize', fs_legend);
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_min_c = min([pain_t_lo; pain_p_lo]);
    y_max_c = max([pain_t_hi; pain_p_hi]);
    y_pad_c = 0.10 * (y_max_c - y_min_c);
    xlim([0, 26]); ylim([y_min_c - y_pad_c, y_max_c + y_pad_c]); box on; hold off;
    
    % Panel D: Treatment Effect Distribution
    subplot(2, 2, 4);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    te_dist = vpop_results.individual_te;
    
    histogram(te_dist, 30, 'FaceColor', col_treat, 'EdgeColor', 'w', 'FaceAlpha', 0.7, 'Normalization', 'pdf');
    
    te_mean = mean(te_dist);
    te_std = std(te_dist);
    x_fit = linspace(min(te_dist)-0.1, max(te_dist)+0.1, 100);
    y_fit = normpdf(x_fit, te_mean, te_std);
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
    
    te_median = median(te_dist);
    pct_benefit = 100 * sum(te_dist > 0) / length(te_dist);
    
    xline(te_median, '--', 'Color', col_plac, 'LineWidth', 2);
    xline(0, ':', 'Color', col_gray, 'LineWidth', 1.5);
    
    % FIX: Use normalized coordinates for labels so they don't get hidden
    text(0.65, 0.92, sprintf('Median = %.3f mm', te_median), ...
         'Units', 'normalized', 'FontSize', 9, 'Color', col_plac);
    text(0.03, 0.85, 'No effect \rightarrow', 'Units', 'normalized', ...
         'FontSize', 8, 'Color', col_gray);
    
    xlabel('3-Year JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    ylabel('Density', 'FontSize', fs_axis);
    title(sprintf('D. Treatment Effect Distribution (N=500, %.1f%% benefit)', pct_benefit), ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    % Add mean/median as annotation instead of cramming into title
    text(0.98, 0.85, sprintf('Mean = %.3f mm', te_mean), ...
         'Units', 'normalized', 'FontSize', 9, 'Color', 'k', ...
         'HorizontalAlignment', 'right');
    xlim([-0.2, 1]); box on; hold off;
    
    saveas(fig3, 'figures/Fig4_virtual_population.png');
    saveas(fig3, 'figures/Fig4_virtual_population.fig');
    fprintf('  Saved: figures/Fig4_virtual_population.png\n');
    fprintf('  Median TE: %.3f mm, %.1f%% show benefit\n', te_median, pct_benefit);
catch ME
    warning('Virtual population figure error: %s', ME.message);
end

%% =============================================================================
%% FIGURE 5: EXTERNAL VALIDATION (2x2) - FIXED TITLE/SUBTITLE OVERLAP
%% =============================================================================
%%
%% FIX SUMMARY:
%%   Panels A-D: Merged title() and subtitle text() into two-line title
%%               using cell array {'Main Title'; 'Subtitle'} to prevent overlap.
%%   Panel D:    Moved threshold labels left (x=2.3) to prevent right-edge
%%               clipping, and increased ylim to accommodate label text.
%%

fprintf('Generating Figure 5: External Validation...\n');

try
    load('model_output/validation_results.mat', 'validation_results');
    
    fig4 = figure('Position', [100, 50, 1000, 800], 'Color', 'w');
    
    col_guide_t = [0.3, 0.6, 0.8];
    col_guide_p = [0.6, 0.6, 0.6];
    col_moves = [0.7, 0.3, 0.5];
    col_celecoxib = [0.9, 0.7, 0.2];
    
    % Get data
    guide_obs = validation_results.guide_observed;
    guide_pred = validation_results.guide_predictions;
    moves_obs = validation_results.moves_observed;
    moves_pred = validation_results.moves_predictions;
    metrics_table = validation_results.metrics;
    
    % Separate GUIDE observed by arm
    guide_obs_treat = guide_obs(contains(guide_obs.arm, 'GS'), :);
    guide_obs_plac = guide_obs(contains(guide_obs.arm, 'Placebo'), :);
    
    % Separate GUIDE predictions by arm
    guide_pred_treat = guide_pred(contains(guide_pred.arm, 'GS'), :);
    guide_pred_plac = guide_pred(contains(guide_pred.arm, 'Placebo'), :);
    
    % Separate MOVES observed by arm
    moves_obs_comb = moves_obs(contains(moves_obs.arm, 'CS+GH'), :);
    moves_obs_celec = moves_obs(contains(moves_obs.arm, 'Celecoxib'), :);
    
    % Panel A: GUIDE 2007
    subplot(2, 2, 1);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    % Model predictions with PI - Treatment
    fill([guide_pred_treat.week; flipud(guide_pred_treat.week)], ...
         [guide_pred_treat.pred_p5; flipud(guide_pred_treat.pred_p95)], ...
         col_guide_t, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(guide_pred_treat.week, guide_pred_treat.pred_mean, 'Color', col_guide_t, 'LineWidth', 2);
    
    % Model predictions with PI - Placebo
    fill([guide_pred_plac.week; flipud(guide_pred_plac.week)], ...
         [guide_pred_plac.pred_p5; flipud(guide_pred_plac.pred_p95)], ...
         col_guide_p, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(guide_pred_plac.week, guide_pred_plac.pred_mean, 'Color', col_guide_p, 'LineWidth', 2);
    
    % Observed data with error bars
    errorbar(guide_obs_treat.week, guide_obs_treat.pain_mean, guide_obs_treat.pain_se, 'o', ...
             'Color', col_guide_t, 'MarkerFaceColor', col_guide_t, 'LineWidth', 1.5);
    errorbar(guide_obs_plac.week, guide_obs_plac.pain_mean, guide_obs_plac.pain_se, 's', ...
             'Color', col_guide_p, 'MarkerFaceColor', col_guide_p, 'LineWidth', 1.5);
    
    xlabel('Time (weeks)', 'FontSize', fs_axis);
    ylabel('WOMAC Pain Score (0-100)', 'FontSize', fs_axis);
    % FIX: Two-line title to avoid overlap with subtitle
    title({'GUIDE 2007 Trial'; 'True Placebo-Controlled External Validation'}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GS 1500mg', 'Placebo'}, 'Location', 'northeast', 'FontSize', fs_legend);
    text(0, -0.12, 'A', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
    xlim([0, 26]); ylim([10, 50]); box on; hold off;
    
    % Panel B: MOVES 2015
    subplot(2, 2, 2);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    % Model prediction with PI
    fill([moves_pred.week; flipud(moves_pred.week)], ...
         [moves_pred.pred_p5; flipud(moves_pred.pred_p95)], ...
         col_moves, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(moves_pred.week, moves_pred.pred_mean, '--', 'Color', col_moves, 'LineWidth', 2);
    
    % Observed CS+GH
    plot(moves_obs_comb.week, moves_obs_comb.pain_mean, '-o', 'Color', col_moves, 'MarkerFaceColor', col_moves, 'LineWidth', 2);
    
    % Celecoxib reference
    plot(moves_obs_celec.week, moves_obs_celec.pain_mean, ':^', 'Color', col_celecoxib, 'MarkerFaceColor', col_celecoxib, 'LineWidth', 1.5);
    
    xlabel('Time (weeks)', 'FontSize', fs_axis);
    ylabel('WOMAC Pain Score (0-100)', 'FontSize', fs_axis);
    % FIX: Two-line title to avoid overlap with subtitle
    title({'MOVES 2015 Trial'; 'Active Comparator (No Placebo Arm)'}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'Model (GH only)', 'CS+GH Observed', 'Celecoxib (ref)'}, 'Location', 'northeast', 'FontSize', fs_legend-1);
    text(0, -0.12, 'B', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
    xlim([0, 28]); ylim([30, 80]); box on; hold off;
    
    % Panel C: Observed vs Predicted - FIXED: Match weeks properly
    subplot(2, 2, 3);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    % For GUIDE Treatment: find predictions at observed weeks
    guide_treat_obs_weeks = guide_obs_treat.week;
    guide_treat_pred_at_obs = zeros(size(guide_treat_obs_weeks));
    for i = 1:length(guide_treat_obs_weeks)
        idx = find(guide_pred_treat.week == guide_treat_obs_weeks(i));
        if ~isempty(idx)
            guide_treat_pred_at_obs(i) = guide_pred_treat.pred_mean(idx);
        end
    end
    
    % For GUIDE Placebo
    guide_plac_obs_weeks = guide_obs_plac.week;
    guide_plac_pred_at_obs = zeros(size(guide_plac_obs_weeks));
    for i = 1:length(guide_plac_obs_weeks)
        idx = find(guide_pred_plac.week == guide_plac_obs_weeks(i));
        if ~isempty(idx)
            guide_plac_pred_at_obs(i) = guide_pred_plac.pred_mean(idx);
        end
    end
    
    % For MOVES
    moves_obs_weeks = moves_obs_comb.week;
    moves_pred_at_obs = zeros(size(moves_obs_weeks));
    for i = 1:length(moves_obs_weeks)
        idx = find(moves_pred.week == moves_obs_weeks(i));
        if ~isempty(idx)
            moves_pred_at_obs(i) = moves_pred.pred_mean(idx);
        end
    end
    
    % Plot scatter
    scatter(guide_obs_treat.pain_mean, guide_treat_pred_at_obs, 80, col_guide_t, 'filled', 'o');
    scatter(guide_obs_plac.pain_mean, guide_plac_pred_at_obs, 80, col_guide_p, 'filled', 's');
    scatter(moves_obs_comb.pain_mean, moves_pred_at_obs, 80, col_moves, 'filled', '^');
    
    % Unity line
    plot([0, 80], [0, 80], '--', 'Color', col_gray, 'LineWidth', 1.5);
    
    xlabel('Observed Pain Score', 'FontSize', fs_axis);
    ylabel('Predicted Pain Score', 'FontSize', fs_axis);
    % FIX: Two-line title to avoid overlap with subtitle
    title({'Model Prediction Accuracy'; 'Combined External Validation'}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GS 1500mg (GUIDE)', 'Placebo (GUIDE)', 'CS+GH (MOVES)'}, 'Location', 'southeast', 'FontSize', fs_legend-1);
    text(0, -0.12, 'C', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
    xlim([0, 80]); ylim([0, 80]); axis square; box on; hold off;
    
    % Panel D: Validation Metrics
    subplot(2, 2, 4);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    % Read metrics from validation_results.metrics (2-row table from M10)
    % Row 1 = GUIDE 2007, Row 2 = MOVES 2015
    % Use direct row indexing for robustness
    mae_g  = metrics_table.MAE(1);
    rmse_g = metrics_table.RMSE(1);
    mae_m  = metrics_table.MAE(2);
    rmse_m = metrics_table.RMSE(2);
    
    % 2x2 matrix: rows = [MAE; RMSE], cols = [GUIDE, MOVES]
    bar_data = [mae_g, mae_m; rmse_g, rmse_m];
    
    b = bar(1:2, bar_data, 'grouped');
    b(1).FaceColor = col_guide_t;
    b(2).FaceColor = col_moves;
    
    yline(5, '--', 'Color', [0, 0.6, 0], 'LineWidth', 1.5);
    text(0.4, 5.5, 'Excellent (<5 pts)', 'FontSize', 8, 'Color', [0, 0.6, 0]);
    yline(10, '--', 'Color', col_celecoxib, 'LineWidth', 1.5);
    text(0.4, 10.5, 'Acceptable (<10 pts)', 'FontSize', 8, 'Color', col_celecoxib);
    
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'MAE', 'RMSE'});
    xlim([0.4, 2.6]);
    ylabel('Error (WOMAC points)', 'FontSize', fs_axis);
    title({'Validation Metrics'; 'Prediction Error (points)'}, ...
          'FontSize', fs_title, 'FontWeight', 'bold');
    legend({'GUIDE 2007', 'MOVES 2015'}, 'Location', 'northwest', 'FontSize', fs_legend);
    text(0, -0.12, 'D', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold');
    ylim([0, 18]); box on; hold off;
    
    saveas(fig4, 'figures/Fig5_validation.png');
    saveas(fig4, 'figures/Fig5_validation.fig');
    fprintf('  Saved: figures/Fig5_validation.png\n');
catch ME
    warning('Validation figure error: %s', ME.message);
end

%% =============================================================================
%% FIGURE 6: TRANSLATIONAL SIMULATIONS (2x2)
%% =============================================================================

fprintf('Generating Figure 6: Translational Simulations...\n');

try
    load('model_output/translational_simulations.mat', 'trans_results');
    
    fig5 = figure('Position', [100, 50, 1000, 800], 'Color', 'w');
    
    obs_te = sim_t.JSW_change(end) - sim_p.JSW_change(end);
    
    % Panel A: Dose-Response
    subplot(2, 2, 1);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    dr = trans_results.dose_response;
    
    fill([dr.dose; flipud(dr.dose)], [dr.jsw_te_p5; flipud(dr.jsw_te_p95)], col_treat, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(dr.dose, dr.jsw_te_mean, '-o', 'Color', col_treat, 'LineWidth', 2, 'MarkerFaceColor', col_treat);
    
    yline(0.1, '--', 'Color', col_gray, 'LineWidth', 1.5);
    text(200, 0.105, 'MCID = 0.1 mm', 'FontSize', 8, 'Color', col_gray);
    
    yline(obs_te, ':', 'Color', col_plac, 'LineWidth', 1.5);
    text(2500, obs_te+0.01, 'Observed TE', 'FontSize', 8, 'Color', col_plac);
    
    xlabel('Glucosamine Sulfate Dose (mg/day)', 'FontSize', fs_axis);
    ylabel('JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    title('A. Dose-Response Relationship', 'FontSize', fs_title, 'FontWeight', 'bold');
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_max_dr = max(dr.jsw_te_p95) * 1.12;
    xlim([0, 3200]); ylim([0, y_max_dr]); box on; hold off;
    
    % Panel B: Bioavailability
    subplot(2, 2, 2);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    bio = trans_results.bioavailability;
    col_bio = [0.2, 0.6, 0.4];
    
    fill([bio.F_percent; flipud(bio.F_percent)], [bio.jsw_te_p5; flipud(bio.jsw_te_p95)], col_bio, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(bio.F_percent, bio.jsw_te_mean, '-o', 'Color', col_bio, 'LineWidth', 2, 'MarkerFaceColor', col_bio);
    
    yline(obs_te, '--', 'Color', col_gray, 'LineWidth', 1.5);
    text(55, obs_te+0.01, 'Current clinical TE', 'FontSize', 8, 'Color', col_gray);
    
    text(22, bio.jsw_te_mean(1)-0.02, 'Current (22%)', 'FontSize', 7);
    
    xlabel('Bioavailability (%)', 'FontSize', fs_axis);
    ylabel('JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    title('B. Bioavailability Enhancement', 'FontSize', fs_title, 'FontWeight', 'bold');
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_max_bio = max(bio.jsw_te_p95) * 1.12;
    xlim([15, 65]); ylim([0, y_max_bio]); box on; hold off;
    
    % Panel C: Duration
    subplot(2, 2, 3);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    dur = trans_results.duration;
    col_dur = [0.5, 0.3, 0.7];
    
    fill([dur.year; flipud(dur.year)], [dur.jsw_te_p5; flipud(dur.jsw_te_p95)], col_dur, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(dur.year, dur.jsw_te_mean, '-o', 'Color', col_dur, 'LineWidth', 2, 'MarkerFaceColor', col_dur);
    
    yline(0.1, '--', 'Color', col_gray, 'LineWidth', 1.5);
    text(0.2, 0.105, 'MCID', 'FontSize', 8, 'Color', col_gray);
    
    yline(obs_te, ':', 'Color', col_plac, 'LineWidth', 1.5);
    text(4, obs_te+0.01, 'Observed TE (3y)', 'FontSize', 8, 'Color', col_plac);
    
    for i = [1, 3, 5]
        if i <= height(dur)
            text(dur.year(i)+0.1, dur.jsw_te_mean(i)+0.015, sprintf('%.2f', dur.jsw_te_mean(i)), 'FontSize', 7);
        end
    end
    
    xlabel('Treatment Duration (years)', 'FontSize', fs_axis);
    ylabel('Cumulative JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    title('C. Duration Optimization', 'FontSize', fs_title, 'FontWeight', 'bold');
    % FIX v3: Dynamic y-axis from actual PI bounds with 10% padding
    y_max_dur = max(dur.jsw_te_p95) * 1.12;
    xlim([0, 5.5]); ylim([0, y_max_dur]); box on; hold off;
    
    % Panel D: Patient Stratification
    subplot(2, 2, 4);
    hold on;
    set(gca, 'Color', [0.95, 0.95, 0.95]);
    
    strat = trans_results.stratification;
    colors_strat = [0.3, 0.6, 0.9; 0.2, 0.7, 0.4; 0.95, 0.85, 0.2; 0.9, 0.5, 0.2; 0.85, 0.4, 0.3];
    
    x = 1:height(strat);
    for i = 1:height(strat)
        color_idx = min(i, size(colors_strat, 1));
        bar(i, strat.jsw_te_mean(i), 'FaceColor', colors_strat(color_idx,:), 'EdgeColor', 'none');
        % Use actual 90% PI from simulation (not fake ±30%)
        err_lo = strat.jsw_te_mean(i) - strat.jsw_te_p5(i);
        err_hi = strat.jsw_te_p95(i) - strat.jsw_te_mean(i);
        errorbar(i, strat.jsw_te_mean(i), err_lo, err_hi, 'k', 'LineWidth', 1);
        text(i, strat.jsw_te_p95(i) + 0.025, ...
             sprintf('%.0f%%', strat.responder_pct(i)), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
    
    yline(0.1, '--', 'Color', col_gray, 'LineWidth', 1.5);
    text(height(strat)+0.3, 0.105, 'MCID', 'FontSize', 8, 'Color', col_gray);
    
    set(gca, 'XTick', x, 'XTickLabel', strat.subgroup);
    xtickangle(20);
    xlabel('Patient Subgroup (OA Severity)', 'FontSize', fs_axis);
    ylabel('JSW Treatment Effect (mm)', 'FontSize', fs_axis);
    title('D. Patient Stratification', 'FontSize', fs_title, 'FontWeight', 'bold');
    % FIX v3: Dynamic y-axis to accommodate error bars + responder % labels
    y_max_strat = max(strat.jsw_te_p95) * 1.20;  % extra space for % labels above error bars
    ylim([0, y_max_strat]); box on; hold off;
    
    saveas(fig5, 'figures/Fig6_translational.png');
    saveas(fig5, 'figures/Fig6_translational.fig');
    fprintf('  Saved: figures/Fig6_translational.png\n');
catch ME
    warning('Translational figure error: %s', ME.message);
end

%% =============================================================================
%% COMPLETION
%% =============================================================================

fprintf('=============================================================\n');
fprintf('ALL 5 FIGURES COMPLETE\n');
fprintf('=============================================================\n');
fprintf('Saved to: figures/\n');
fprintf('  - Fig2_calibration.png\n');
fprintf('  - Fig3_sensitivity.png\n');
fprintf('  - Fig4_virtual_population.png\n');
fprintf('  - Fig5_validation.png\n');
fprintf('  - Fig6_translational.png\n');
fprintf('=============================================================\n');
