%% =============================================================================
%% M16_IC50_forward_sweep.m
%% Forward IC50,deg sensitivity (calibrated parameters FIXED) - Reviewer 2
%% =============================================================================
%% Companion to M15. Holds the published calibrated parameters fixed and sweeps
%% IC50_deg from 3 -> 1000 ug/mL, showing how (A) achieved degradation
%% inhibition at the measured synovial concentration and (B) the predicted
%% 3-year JSW treatment effect decline as potency becomes biologically
%% realistic. This is the intuitive illustration; M15 is the full re-fit.
%%
%% Run AFTER M04 (needs model_output/calibration_results.mat).
%% =============================================================================

function fwd = M16_IC50_forward_sweep()

    fprintf('=== FORWARD IC50 SWEEP (calibrated params fixed) ===\n\n');

    load('model_output/calibration_results.mat', 'calibration_results');
    p0 = calibration_results.params;
    if ~isfield(p0,'Emax_syn'), p0.Emax_syn = 0.15; end
    if ~isfield(p0,'EC50_syn'), p0.EC50_syn = 2.0;  end

    IC50_grid = [3, 10, 30, 100, 300, 1000];
    nG = numel(IC50_grid);

    C_model    = 0.30;   % ug/mL model reference
    C_measured = 0.78;   % ug/mL Persiani measured synovial (4.34 uM)

    I_model    = zeros(nG,1);
    I_measured = zeros(nG,1);
    te_jsw     = zeros(nG,1);

    fprintf('%-8s | %-12s | %-12s | %-10s\n','IC50','Iinh@0.30','Iinh@0.78','JSW TE');
    fprintf('%s\n', repmat('-',1,48));
    for g = 1:nG
        p = p0; p.IC50_deg = IC50_grid(g);
        I_model(g)    = p.Imax_deg * C_model    / (IC50_grid(g) + C_model);
        I_measured(g) = p.Imax_deg * C_measured / (IC50_grid(g) + C_measured);
        sp = simulate_glucosamine(p,1500,156,'sulfate',false);
        st = simulate_glucosamine(p,1500,156,'sulfate',true);
        te_jsw(g) = st.JSW_change(end) - sp.JSW_change(end);
        fprintf('%-8.0f | %-12.2f | %-12.2f | %-10.3f\n', ...
                IC50_grid(g), 100*I_model(g), 100*I_measured(g), te_jsw(g));
    end

    fig = figure('Position',[100 100 1000 420],'Color','w');
    subplot(1,2,1);
    semilogx(IC50_grid, 100*I_measured, '-o','LineWidth',2,'Color',[0.2 0.4 0.7], ...
             'MarkerFaceColor',[0.2 0.4 0.7]); hold on;
    semilogx(IC50_grid, 100*I_model, '--s','LineWidth',1.5,'Color',[0.6 0.6 0.6], ...
             'MarkerFaceColor',[0.6 0.6 0.6]);
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)');
    ylabel('Degradation inhibition (%)');
    title('A. Achieved anti-catabolic effect');
    legend({'At measured synovial (0.78 \mug/mL)','At model reference (0.30 \mug/mL)'}, ...
           'Location','northeast','FontSize',8);
    xline(100,':','Largo effective range','Color',[0.7 0.3 0.3], ...
          'LabelVerticalAlignment','bottom','FontSize',8);

    subplot(1,2,2);
    semilogx(IC50_grid, te_jsw, '-o','LineWidth',2,'Color',[0.2 0.4 0.7], ...
             'MarkerFaceColor',[0.2 0.4 0.7]); hold on;
    yline(0.1,'--','MCID = 0.1 mm','Color',[0.5 0.5 0.5]);
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)');
    ylabel('3-year JSW treatment effect (mm)');
    title('B. Predicted structural benefit');

    sgtitle('Figure S6. Forward IC_{50,deg} sensitivity (calibrated parameters fixed)', ...
            'FontWeight','bold');

    if ~exist('model_output','dir'), mkdir('model_output'); end
    exportgraphics(fig,'model_output/FigS6_IC50_forward_sweep.tiff','Resolution',300);

    fwd = struct('IC50_grid',IC50_grid,'I_model',I_model, ...
                 'I_measured',I_measured,'te_jsw',te_jsw);
    save('model_output/ic50_forward_sweep_results.mat','fwd');
    fprintf('\nSaved: model_output/ic50_forward_sweep_results.mat\n');
    fprintf('Figure: model_output/FigS6_IC50_forward_sweep.tiff\n');
end