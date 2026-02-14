%% =============================================================================
%% M03_calibration_data.m
%% Calibration Data Preparation
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%% MATLAB SimBiology Implementation
%%
%% Data Sources:
%%   - Reginster 2001: 3-year JSW trial (glucosamine sulfate)
%%   - Pavelka 2002: 3-year JSW trial (glucosamine sulfate)
%%   - GAIT 2006: 24-week pain trial (glucosamine HCl)
%%
%%   - GAIT WOMAC normalization corrected and documented
%%   - Reginster GS Year 1 JSW change verified
%%   - All values cross-referenced against source publications
%%
%% Output: data_clean/calibration_data.mat
%%
%% =============================================================================

function calibration_data = M03_calibration_data()
    
    fprintf('=============================================================\n');
    fprintf('CALIBRATION DATA PREPARATION (CORRECTED)\n');
    fprintf('=============================================================\n\n');
    
    %% =============================================================================
    %% 1. JSW DATA - REGINSTER 2001
    %% =============================================================================
    %% Source: Reginster et al., Lancet 2001; 357:251-256
    %% Table 2, ITT population
    %%
    %% Values verified against source:
    %%   Placebo: BL=5.45, Y1=5.39, Y3=5.14 mm (Table 2)
    %%   GS:      BL=5.51, Y1=5.55, Y3=5.45 mm (Table 2)
    %%
        %% =============================================================================
    
    fprintf('Loading Reginster 2001 data...\n');
    
    reginster_study = repmat({'Reginster 2001'}, 6, 1);
    reginster_week = [0; 52; 156; 0; 52; 156];
    reginster_arm = {'Placebo'; 'Placebo'; 'Placebo'; 'GS 1500mg'; 'GS 1500mg'; 'GS 1500mg'};
    reginster_jsw = [5.45; 5.39; 5.14; 5.51; 5.55; 5.45];
    reginster_jsw_sd = [1.35; 1.35; 1.45; 1.40; 1.40; 1.45];
    reginster_n = [106; 106; 106; 106; 106; 106];
    
    % Calculate change from baseline
    reginster_baseline_placebo = reginster_jsw(1);
    reginster_baseline_treat = reginster_jsw(4);
    reginster_change = zeros(6, 1);
    for i = 1:3
        reginster_change(i) = reginster_jsw(i) - reginster_baseline_placebo;
    end
    for i = 4:6
        reginster_change(i) = reginster_jsw(i) - reginster_baseline_treat;
    end
    
    % VERIFICATION OUTPUT
    fprintf('  Reginster GS Y1 change: %.2f mm (source paper: +0.04 mm)\n', reginster_change(5));
    fprintf('  Reginster GS Y3 change: %.2f mm (source paper: -0.06 mm)\n', reginster_change(6));
    fprintf('  Reginster Placebo Y3 change: %.2f mm (source paper: -0.31 mm)\n', reginster_change(3));
    
    reginster_jsw_data = table(reginster_study, reginster_week, reginster_arm, ...
                               reginster_jsw, reginster_jsw_sd, reginster_n, reginster_change, ...
                               'VariableNames', {'study', 'week', 'arm', 'jsw_mm', 'jsw_sd', 'n', 'jsw_change_mm'});
    
    %% =============================================================================
    %% 2. JSW DATA - PAVELKA 2002
    %% =============================================================================
    %% Source: Pavelka et al., Arch Intern Med 2002; 162:2113-2123
    %% Table 3, ITT population
    %%
    %% Values verified against source:
    %%   Placebo: BL=3.95, Y1=3.92, Y3=3.76 mm
    %%   GS:      BL=3.87, Y1=4.00, Y3=3.91 mm
    %% =============================================================================
    
    fprintf('Loading Pavelka 2002 data...\n');
    
    pavelka_study = repmat({'Pavelka 2002'}, 6, 1);
    pavelka_week = [0; 52; 156; 0; 52; 156];
    pavelka_arm = {'Placebo'; 'Placebo'; 'Placebo'; 'GS 1500mg'; 'GS 1500mg'; 'GS 1500mg'};
    pavelka_jsw = [3.95; 3.92; 3.76; 3.87; 4.00; 3.91];
    pavelka_jsw_sd = [1.51; 1.48; 1.50; 1.61; 1.57; 1.60];
    pavelka_n = [101; 101; 101; 101; 101; 101];
    
    % Calculate change from baseline
    pavelka_baseline_placebo = pavelka_jsw(1);
    pavelka_baseline_treat = pavelka_jsw(4);
    pavelka_change = zeros(6, 1);
    for i = 1:3
        pavelka_change(i) = pavelka_jsw(i) - pavelka_baseline_placebo;
    end
    for i = 4:6
        pavelka_change(i) = pavelka_jsw(i) - pavelka_baseline_treat;
    end
    
    pavelka_jsw_data = table(pavelka_study, pavelka_week, pavelka_arm, ...
                             pavelka_jsw, pavelka_jsw_sd, pavelka_n, pavelka_change, ...
                             'VariableNames', {'study', 'week', 'arm', 'jsw_mm', 'jsw_sd', 'n', 'jsw_change_mm'});
    
    %% =============================================================================
    %% 3. COMBINE JSW DATA
    %% =============================================================================
    
    jsw_data = [reginster_jsw_data; pavelka_jsw_data];
    
    fprintf('\nJSW Data Summary:\n');
    disp(jsw_data(:, {'study', 'week', 'arm', 'jsw_change_mm'}));
    
    %% =============================================================================
    %% 4. WOMAC PAIN DATA - GAIT 2006
    %% =============================================================================
    %%
    %% Source: Clegg et al., NEJM 2006; 354:795-808
    %% Table 2 and Figure 2
    %%
    %% CRITICAL NORMALIZATION NOTE:
    %%   GAIT reports WOMAC pain on 0-500 VAS scale (5 items x 0-100 mm VAS each)
    %%   To normalize to 0-100: divide by 5
    %%
    %%   GAIT Table 2 reports:
    %%     Placebo baseline: 237.1 (SD 90.8), N=313
    %%     GH baseline:      233.2 (SD 87.5), N=317
    %%
    %%   Normalized to 0-100:
    %%     Placebo baseline: 237.1/5 = 47.42
    %%     GH baseline:      233.2/5 = 46.64
    %%
        %%   These values (256, 252) do NOT match the GAIT publication.
    %%
        %%   Timepoints other than baseline are read from Figure 2 of GAIT paper.
    %%   
        %% =============================================================================
    
    fprintf('\nLoading GAIT 2006 pain data (CORRECTED normalization)...\n');
    
    gait_study = repmat({'GAIT 2006'}, 10, 1);
    gait_week = [0; 4; 8; 16; 24; 0; 4; 8; 16; 24];
    gait_arm = {'Placebo'; 'Placebo'; 'Placebo'; 'Placebo'; 'Placebo'; ...
                'GH 1500mg (HCl)'; 'GH 1500mg (HCl)'; 'GH 1500mg (HCl)'; 'GH 1500mg (HCl)'; 'GH 1500mg (HCl)'};
    
    % CORRECTED: Use GAIT Table 2 values on WOMAC 0-500 VAS scale
    % Baseline values from Table 2; follow-up from Figure 2 (approximate)
    % Placebo: BL=237.1, approximate follow-ups from Fig 2
    % GH:      BL=233.2, approximate follow-ups from Fig 2
    gait_womac = [237.1; 176; 160; 149; 143;    % Placebo trajectory
                  233.2; 172; 154; 140; 137];    % GH trajectory
    gait_womac_sd = [90.8; 105; 110; 113; 116;  % Placebo SDs (BL from Table 2, rest approximate)
                     87.5; 103; 108; 111; 114];  % GH SDs
    gait_n = [313; 313; 313; 313; 313; 317; 317; 317; 317; 317];
    
    % Normalize to 0-100 scale (WOMAC pain 0-500 VAS -> 0-100)
    gait_womac_norm = gait_womac / 5;
    gait_womac_sd_norm = gait_womac_sd / 5;
    
    fprintf('\n  CORRECTED GAIT normalized values:\n');
    fprintf('    Placebo baseline: %.1f (was 51.2 in old code, source paper: 47.4)\n', gait_womac_norm(1));
    fprintf('    GH baseline:      %.1f (was 50.4 in old code, source paper: 46.6)\n', gait_womac_norm(6));
    fprintf('    Placebo W24:      %.1f\n', gait_womac_norm(5));
    fprintf('    GH W24:           %.1f\n', gait_womac_norm(10));
    
    womac_pain_data = table(gait_study, gait_week, gait_arm, ...
                            gait_womac, gait_womac_sd, gait_n, ...
                            gait_womac_norm, gait_womac_sd_norm, ...
                            'VariableNames', {'study', 'week', 'arm', ...
                                             'womac_pain', 'womac_pain_sd', 'n', ...
                                             'womac_pain_normalized', 'womac_pain_sd_normalized'});
    
    fprintf('\nWOMAC Pain Data Summary:\n');
    disp(womac_pain_data(:, {'study', 'week', 'arm', 'womac_pain_normalized'}));
    
    %% =============================================================================
    %% 5. CREATE CALIBRATION TARGETS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('CALIBRATION TARGETS\n');
    fprintf('=============================================================\n\n');
    
    % JSW targets (exclude baseline)
    jsw_targets = jsw_data(jsw_data.week > 0, :);
    
    fprintf('JSW Targets (%d observations):\n', height(jsw_targets));
    disp(jsw_targets(:, {'study', 'week', 'arm', 'jsw_change_mm'}));
    
    % WOMAC targets
    womac_targets = womac_pain_data;
    
    fprintf('\nWOMAC Pain Targets (%d observations):\n', height(womac_targets));
    disp(womac_targets(:, {'study', 'week', 'arm', 'womac_pain_normalized'}));
    
    %% =============================================================================
    %% 6. CALCULATE OBSERVED TREATMENT EFFECTS
    %% =============================================================================
    
    fprintf('\n=============================================================\n');
    fprintf('OBSERVED TREATMENT EFFECTS\n');
    fprintf('=============================================================\n\n');
    
    % JSW treatment effects
    fprintf('JSW Treatment Effects:\n');
    
    studies = unique(jsw_data.study);
    for s = 1:length(studies)
        study_name = studies{s};
        study_data = jsw_data(strcmp(jsw_data.study, study_name) & jsw_data.week == 156, :);
        placebo_idx = contains(study_data.arm, 'Placebo');
        treat_idx = contains(study_data.arm, 'GS');
        
        if sum(placebo_idx) > 0 && sum(treat_idx) > 0
            placebo_change = study_data.jsw_change_mm(placebo_idx);
            treat_change = study_data.jsw_change_mm(treat_idx);
            te = treat_change - placebo_change;
            fprintf('  %s (3-year): Placebo = %.2f mm, Treatment = %.2f mm, TE = %.2f mm\n', ...
                    study_name, placebo_change, treat_change, te);
        end
    end
    
    % Pain treatment effect
    fprintf('\nWOMAC Pain Treatment Effect:\n');
    pain_24 = womac_pain_data(womac_pain_data.week == 24, :);
    placebo_pain = pain_24.womac_pain_normalized(contains(pain_24.arm, 'Placebo'));
    treat_pain = pain_24.womac_pain_normalized(contains(pain_24.arm, 'GH'));
    fprintf('  GAIT 2006 (24-week): Placebo = %.1f, Treatment = %.1f, Diff = %.1f\n', ...
            placebo_pain, treat_pain, placebo_pain - treat_pain);
    
    %% =============================================================================
    %% 7. SAVE CALIBRATION DATA
    %% =============================================================================
    
    if ~exist('data_clean', 'dir')
        mkdir('data_clean');
    end
    
    calibration_data.jsw_data = jsw_data;
    calibration_data.womac_pain_data = womac_pain_data;
    calibration_data.jsw_targets = jsw_targets;
    calibration_data.womac_targets = womac_targets;
    
    save('data_clean/calibration_data.mat', 'calibration_data');
    
    fprintf('\n=============================================================\n');
    fprintf('Data saved to: data_clean/calibration_data.mat\n');
    fprintf('=============================================================\n');
    
    
end

%% =============================================================================
%% HELPER FUNCTIONS
%% =============================================================================

function [jsw_targets, womac_targets] = get_calibration_targets()
    % Load calibration data and return targets
    load('data_clean/calibration_data.mat', 'calibration_data');
    jsw_targets = calibration_data.jsw_targets;
    womac_targets = calibration_data.womac_targets;
end

function targets_array = create_targets_array(jsw_targets, womac_targets)
    % Create a combined array of calibration targets for optimization
    % Format: [week, arm_code, observed_value, SE, endpoint_type]
    % arm_code: 0 = placebo, 1 = treatment
    % endpoint_type: 1 = JSW, 2 = WOMAC
    
    n_jsw = height(jsw_targets);
    n_womac = height(womac_targets);
    
    targets_array = zeros(n_jsw + n_womac, 5);
    
    % JSW targets
    for i = 1:n_jsw
        targets_array(i, 1) = jsw_targets.week(i);
        targets_array(i, 2) = contains(jsw_targets.arm{i}, 'GS');  % 1 for treatment
        targets_array(i, 3) = jsw_targets.jsw_change_mm(i);
        targets_array(i, 4) = 0.08;  % SE
        targets_array(i, 5) = 1;     % JSW endpoint
    end
    
    % WOMAC targets
    for i = 1:n_womac
        targets_array(n_jsw + i, 1) = womac_targets.week(i);
        targets_array(n_jsw + i, 2) = contains(womac_targets.arm{i}, 'GH');
        targets_array(n_jsw + i, 3) = womac_targets.womac_pain_normalized(i);
        targets_array(n_jsw + i, 4) = 5.0;  % SE
        targets_array(n_jsw + i, 5) = 2;    % WOMAC endpoint
    end
end
