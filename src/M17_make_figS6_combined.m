%% =============================================================================
%% M17_make_figS6_combined.m
%% Combined 4-panel Figure S6 for Reviewer 2 IC50 analysis (CLEAN export)
%% =============================================================================
%% Top row (A,B): forward sweep (calibrated params fixed)
%% Bottom row (C,D): re-calibration (params re-fit at each IC50)
%% Reads results saved by M15 and M16; no re-computation needed.
%% Toolbar suppressed; legend fixed; axes labeled.
%% =============================================================================

function M17_make_figS6_combined()

    load('model_output/ic50_forward_sweep_results.mat','fwd');
    load('model_output/ic50_recalibration_results.mat','ic50_recal');

    IC = fwd.IC50_grid;
    blue = [0.2 0.4 0.7]; grey = [0.6 0.6 0.6];

    fig = figure('Position',[100 100 1100 850],'Color','w');

    % ---- Panel A: achieved inhibition (forward) ----
    axA = subplot(2,2,1);
    semilogx(IC, 100*fwd.I_measured,'-o','LineWidth',2,'Color',blue,'MarkerFaceColor',blue); hold on;
    semilogx(IC, 100*fwd.I_model,'--s','LineWidth',1.5,'Color',grey,'MarkerFaceColor',grey);
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)'); ylabel('Degradation inhibition (%)');
    title('A. Achieved anti-catabolic effect');
    legend({'At measured synovial (0.78 \mug/mL)','At model reference (0.30 \mug/mL)'}, ...
           'Location','northeast','FontSize',8);
    % shade Largo effective range (100-1000) instead of a legend-polluting line
    yl = ylim; patch([100 1200 1200 100],[yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.3 0.3],'FaceAlpha',0.06,'EdgeColor','none');
    text(110, yl(2)*0.9,'Largo effective range','Color',[0.7 0.3 0.3],'FontSize',8);
    uistack(axA,'top');

    % ---- Panel B: predicted benefit (forward) ----
    subplot(2,2,2);
    semilogx(IC, fwd.te_jsw,'-o','LineWidth',2,'Color',blue,'MarkerFaceColor',blue); hold on;
    yline(0.1,'--','MCID = 0.1 mm','Color',grey,'LabelHorizontalAlignment','left');
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)'); ylabel('3-year JSW treatment effect (mm)');
    title('B. Predicted structural benefit (parameters fixed)');

    % ---- Panel C: goodness of fit (re-calibration) ----
    subplot(2,2,3);
    semilogx(IC, ic50_recal.sse,'-o','LineWidth',2,'Color',blue,'MarkerFaceColor',blue);
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)'); ylabel('SSE (lower = better fit)');
    title('C. Fit degrades when re-calibrated');

    % ---- Panel D: Imax pinned + benefit (re-calibration) ----
    subplot(2,2,4);
    yyaxis left;
    semilogx(IC, ic50_recal.refit(:,3),'-o','LineWidth',2,'Color',blue,'MarkerFaceColor',blue);
    ylabel('Re-fit I_{max,deg}'); ylim([0 1.08]);
    yline(1.0,':','I_{max} upper bound','Color',[0.3 0.3 0.3]);
    yyaxis right;
    semilogx(IC, ic50_recal.te_jsw,'-s','LineWidth',1.8);
    ylabel('3-yr JSW TE (mm)');
    grid on; box on; xlim([2 1200]);
    xlabel('IC_{50,deg} (\mug/mL, log scale)');
    title('D. I_{max} pinned at bound; benefit collapses');

    sgtitle('Figure S6. Sensitivity of predicted efficacy to anti-catabolic potency (IC_{50,deg})', ...
            'FontWeight','bold','FontSize',13);

    % ---- CLEAN export: remove toolbars from every axis, then export ----
    set(findall(fig,'Type','axes'),'Toolbar',[]);
    drawnow;
    if ~exist('model_output','dir'), mkdir('model_output'); end
    exportgraphics(fig,'model_output/FigS6_IC50_sensitivity_combined.tiff','Resolution',300);
    fprintf('Saved: model_output/FigS6_IC50_sensitivity_combined.tiff\n');
end