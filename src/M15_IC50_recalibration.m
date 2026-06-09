%% =============================================================================
%% M15_IC50_recalibration.m
%% Wide-range IC50,deg RE-CALIBRATION analysis (Reviewer 2, Points 7-8)
%% =============================================================================
%% At each IC50_deg value across a biologically realistic range, the model is
%% FULLY RE-CALIBRATED (k_syn, k_deg, Imax_deg, Pain_0, Pmax re-fit by fmincon
%% against the same calibration targets and objective used in M06/M07).
%%
%% Reports, per IC50: re-fit parameters, SSE (fit quality), whether Imax_deg
%% saturates at its upper bound, predicted 3-yr JSW treatment effect, % of
%% virtual responders, and achieved degradation inhibition at TWO synovial
%% concentrations:
%%   (a) model reference 0.30 ug/mL
%%   (b) Persiani 2007 MEASURED synovial 0.78 ug/mL (4.34 uM, free-base MW)
%%
%% Run AFTER M03 and M04 (needs calibration data + a calibrated baseline).
%% =============================================================================

function ic50_recal = M15_IC50_recalibration()

    rng(42);
    fprintf('=== WIDE-RANGE IC50 RE-CALIBRATION (Reviewer 2) ===\n\n');

    % --- Load model structure (for bounds) and calibration data ---
    [param_info, ~, ~] = M01_model_structure();
    load('model_output/calibration_results.mat', 'calibration_results');
    base_params = calibration_results.params;
    load('data_clean/calibration_data.mat', 'calibration_data');
    jsw_targets   = calibration_data.jsw_data(calibration_data.jsw_data.week > 0, :);
    womac_targets = calibration_data.womac_pain_data( ...
                      strcmp(calibration_data.womac_pain_data.study,'GAIT 2006'), :);

    % --- Parameters that get re-fit at each IC50 (same 5 as the paper) ---
    fit_names = {'k_syn','k_deg','Imax_deg','Pain_0','Pmax'};
    lb = zeros(1,5); ub = zeros(1,5); x0 = zeros(1,5);
    for i = 1:5
        row = strcmp(param_info.Name, fit_names{i});
        lb(i) = param_info.Lower(row);
        ub(i) = param_info.Upper(row);
        x0(i) = base_params.(fit_names{i});
    end
    Imax_ub = ub(3);   % upper bound on Imax_deg (1.0)

    % --- IC50 sweep (log-spaced, 3 -> 1000 ug/mL) ---
    IC50_grid = [3, 10, 30, 100, 300, 1000];
    nG = numel(IC50_grid);

    refit           = zeros(nG,5);
    sse             = zeros(nG,1);
    imax_at_bnd     = false(nG,1);
    te_jsw          = zeros(nG,1);
    pct_resp        = zeros(nG,1);
    I_drug_ref      = zeros(nG,1);
    I_drug_measured = zeros(nG,1);

    opts = optimoptions('fmincon','Display','off','MaxIterations',300);

    for g = 1:nG
        p = base_params;
        p.IC50_deg = IC50_grid(g);

        obj = @(x) ic50_objective(x, fit_names, p, jsw_targets, womac_targets);

        % Multi-start (3 starts) for robustness
        best = inf; bx = x0;
        starts = [x0; x0.*[1.2 0.8 1 1 1]; x0.*[0.8 1.2 1 0.9 1.1]];
        for s = 1:size(starts,1)
            s0 = min(max(starts(s,:), lb), ub);
            try
                [xx,ff] = fmincon(obj, s0, [],[],[],[], lb, ub, ...
                          @(x) deal( x(1)-x(2), [] ), opts); % constraint: k_syn < k_deg
                if ff < best, best = ff; bx = xx; end
            catch
            end
        end

        refit(g,:)     = bx;
        sse(g)         = best;
        imax_at_bnd(g) = bx(3) >= 0.999*Imax_ub;

        % Predicted 3-yr JSW treatment effect with re-fit params
        pf = setp(base_params, fit_names, bx); pf.IC50_deg = IC50_grid(g);
        sp = simulate_glucosamine(pf,1500,156,'sulfate',false);
        st = simulate_glucosamine(pf,1500,156,'sulfate',true);
        te_jsw(g) = st.JSW_change(end) - sp.JSW_change(end);

        % Achieved inhibition at (a) model reference 0.30 ug/mL and
        % (b) Persiani's MEASURED steady-state synovial conc 0.78 ug/mL (4.34 uM)
        C_model    = 0.30;   % ug/mL  (model's internal reference)
        C_measured = 0.78;   % ug/mL  (Persiani 2007 synovial: 4.34 uM x 179.17/1000 = 0.778)
        I_drug_ref(g)      = bx(3) * C_model    / (IC50_grid(g) + C_model);
        I_drug_measured(g) = bx(3) * C_measured / (IC50_grid(g) + C_measured);

        % % responders in a quick virtual population (N=200)
        pct_resp(g) = quick_responders(pf, 200);

        fprintf(['IC50=%6.0f | SSE=%7.2f | Imax=%.3f%s | k_syn=%.3f k_deg=%.3f ' ...
                 '| TE=%.3f mm | resp=%.0f%% | Iinh@0.30=%.1f%% | Iinh@0.78=%.1f%%\n'], ...
                 IC50_grid(g), best, bx(3), ternary(imax_at_bnd(g),'*',' '), ...
                 bx(1), bx(2), te_jsw(g), pct_resp(g), ...
                 100*I_drug_ref(g), 100*I_drug_measured(g));
    end

    fprintf('\nNote: SSE at IC50=3 is the published fit. Rising SSE = worse fit.\n');
    fprintf('"*" next to Imax indicates it is pinned at its upper bound (cannot compensate).\n');

    % --- Figure ---
    fig = figure('Position',[100 100 1100 400],'Color','w');
    subplot(1,3,1);
    semilogx(IC50_grid, sse, '-o','LineWidth',2); grid on;
    xlabel('IC_{50,deg} (\mug/mL)'); ylabel('SSE (fit quality)');
    title('A. Goodness of fit'); xlim([2 1200]);

    subplot(1,3,2);
    semilogx(IC50_grid, refit(:,3), '-o','LineWidth',2); grid on; hold on;
    yline(Imax_ub,'--','I_{max} bound');
    xlabel('IC_{50,deg} (\mug/mL)'); ylabel('Re-fit I_{max,deg}');
    title('B. I_{max} compensation'); xlim([2 1200]); ylim([0 1.05]);

    subplot(1,3,3);
    semilogx(IC50_grid, te_jsw, '-o','LineWidth',2); grid on; hold on;
    yline(0.1,'--','MCID');
    xlabel('IC_{50,deg} (\mug/mL)'); ylabel('3-yr JSW TE (mm)');
    title('C. Predicted benefit'); xlim([2 1200]);

    sgtitle('Figure S7. Re-calibration across realistic IC_{50,deg} range','FontWeight','bold');
    if ~exist('model_output','dir'), mkdir('model_output'); end
    exportgraphics(fig,'model_output/FigS7_IC50_recalibration.tiff','Resolution',300);

    % --- Save ---
    ic50_recal = struct('IC50_grid',IC50_grid,'refit',refit,'fit_names',{fit_names}, ...
        'sse',sse,'imax_at_bound',imax_at_bnd,'te_jsw',te_jsw,'pct_resp',pct_resp, ...
        'I_drug_ref',I_drug_ref,'I_drug_measured',I_drug_measured);
    save('model_output/ic50_recalibration_results.mat','ic50_recal');
    fprintf('\nSaved: model_output/ic50_recalibration_results.mat\n');
    fprintf('Figure: model_output/FigS7_IC50_recalibration.tiff\n');
end

%% ---- objective: same structure as M06_profile_likelihood compute_objective ----
function s = ic50_objective(x, fit_names, p, jsw_t, womac_t)
    for i = 1:numel(fit_names), p.(fit_names{i}) = x(i); end
    if p.k_deg <= p.k_syn, s = 1e9; return; end
    s = 0;
    try
        sp = simulate_glucosamine(p,1500,156,'sulfate',false);
        st = simulate_glucosamine(p,1500,156,'sulfate',true);
        for i = 1:height(jsw_t)
            wk = jsw_t.week(i); obs = jsw_t.jsw_change_mm(i);
            if contains(jsw_t.arm{i},'Placebo')
                idx = find(sp.week==wk); if ~isempty(idx), pr = sp.JSW_change(idx); else, continue; end
            else
                idx = find(st.week==wk); if ~isempty(idx), pr = st.JSW_change(idx); else, continue; end
            end
            s = s + ((pr-obs)/0.08)^2;
        end
    catch, s = s + 100; end
    try
        sp24 = simulate_glucosamine(p,1500,24,'HCl',false);
        st24 = simulate_glucosamine(p,1500,24,'HCl',true);
        for i = 1:height(womac_t)
            wk = womac_t.week(i); obs = womac_t.womac_pain_normalized(i);
            if contains(womac_t.arm{i},'Placebo')
                idx = find(sp24.week==wk); if ~isempty(idx), pr = sp24.Pain(idx); else, continue; end
            else
                idx = find(st24.week==wk); if ~isempty(idx), pr = st24.Pain(idx); else, continue; end
            end
            s = s + 0.1*((pr-obs)/10)^2;
        end
    catch, s = s + 50; end
end

%% ---- helper: set named parameters from a vector ----
function p = setp(p, names, x)
    for i = 1:numel(names), p.(names{i}) = x(i); end
end

%% ---- helper: quick virtual-population responder fraction ----
function pct = quick_responders(p, N)
    cnt = 0;
    for i = 1:N
        q = p;
        q.k_syn    = p.k_syn  * exp(0.30*randn());
        q.k_deg    = p.k_deg  * exp(0.30*randn());
        q.Imax_deg = min(max(p.Imax_deg*exp(0.25*randn()),0.05),1);
        if q.k_deg <= q.k_syn, q.k_deg = q.k_syn*1.2; end
        sp = simulate_glucosamine(q,1500,156,'sulfate',false);
        st = simulate_glucosamine(q,1500,156,'sulfate',true);
        if (st.JSW_change(end)-sp.JSW_change(end)) > 0.1, cnt = cnt+1; end
    end
    pct = 100*cnt/N;
end

%% ---- helper: inline ternary ----
function out = ternary(c,a,b)
    if c, out = a; else, out = b; end
end