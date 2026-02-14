%% =============================================================================
%% M01_model_structure.m
%% Model Structure and Parameter Definitions
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%%
%% Defines all model parameters, default values, bounds, and optimization flags.
%%
%% Output:
%%   param_info    - Table with parameter metadata
%%   default_params - Structure with default parameter values
%%   param_bounds   - Structure with lower/upper bounds
%%
%% =============================================================================

function [param_info, default_params, param_bounds] = M01_model_structure()

    %% =============================================================================
    %% 1. PARAMETER DEFINITIONS
    %% =============================================================================

    param_names = {
        'F_sulfate', 'F_HCl', 'Vd_apparent', 'R_syn', ...
        'k_syn', 'k_deg', 'GAG_0', 'JSW_0', 'gamma_JSW', ...
        'Imax_deg', 'IC50_deg', ...
        'Emax_syn', 'EC50_syn', ...
        'Pain_0', 'alpha_struct', 'beta_drug', 'EC50_pain', 'Pmax', 'kpl'
    };

    descriptions = {
        'Oral bioavailability (sulfate formulation)', ...
        'Oral bioavailability (HCl formulation)', ...
        'Apparent volume of distribution', ...
        'Synovial fluid to plasma ratio', ...
        'GAG synthesis rate', ...
        'GAG degradation rate', ...
        'Initial GAG content', ...
        'Baseline joint space width', ...
        'JSW-GAG power relationship', ...
        'Max inhibition of degradation', ...
        'IC50 for anti-catabolic effect', ...
        'Max synthesis stimulation', ...
        'EC50 for synthesis stimulation', ...
        'Baseline WOMAC pain', ...
        'Pain from cartilage loss', ...
        'Direct analgesic effect', ...
        'EC50 for pain effect', ...
        'Maximum placebo response', ...
        'Placebo effect rate'
    };

    units = {
        'fraction', 'fraction', 'L', 'ratio', ...
        '/year', '/year', 'fraction', 'mm', 'dimensionless', ...
        'fraction', 'ug/mL', ...
        'fraction', 'ug/mL', ...
        '0-100', 'points/%GAG', 'points', 'uM', 'points', '/day'
    };

    % Default (calibrated) values
    defaults = [0.22, 0.11, 200, 0.25, ...               % F_sulfate, F_HCl, Vd, R_syn
                0.175258, 0.207263, 0.95, 4.2, 1.0, ...   % k_syn, k_deg, GAG_0, JSW_0, gamma
                1.0, 3.0, ...                              % Imax_deg, IC50_deg
                0.15, 2.0, ...                             % Emax_syn, EC50_syn
                46.651892, 150, 5, 3, 17.624509, 0.02];   % Pain params

    % Lower bounds
    lower_bounds = [0.15, 0.08, 100, 0.10, ...
                    0.05, 0.06, 0.85, 3.5, 0.5, ...
                    0.10, 0.5, ...
                    0.0, 0.5, ...
                    35, 50, 0, 0.5, 10, 0.005];

    % Upper bounds
    upper_bounds = [0.30, 0.15, 400, 0.50, ...
                    0.20, 0.25, 1.00, 5.5, 2.0, ...
                    1.00, 10.0, ...
                    0.50, 10.0, ...
                    55, 300, 15, 10, 30, 0.05];

    % Optimized flags (1 = estimated, 0 = fixed)
    optimized = [0, 0, 0, 0, ...       % PK params fixed
                 1, 1, 0, 0, 0, ...    % k_syn, k_deg estimated
                 1, 0, ...             % Imax_deg estimated
                 0, 0, ...             % Emax_syn, EC50_syn fixed
                 1, 0, 0, 0, 1, 0];   % Pain_0, Pmax estimated

    sources = {
        'Persiani 2005', 'Persiani 2005', 'Setnikar 1993', 'Assumed', ...
        'Estimated', 'Estimated', 'Assumed (mild OA)', 'Reginster 2001', 'Assumed', ...
        'Estimated', 'In vitro data', ...
        'Fixed (in vitro)', 'Fixed (in vitro)', ...
        'GAIT 2006', 'Assumed', 'Assumed', 'Assumed', 'Estimated', 'Calibrated'
    };

    %% =============================================================================
    %% 2. CREATE PARAMETER INFO TABLE
    %% =============================================================================

    n_params = length(param_names);

    param_info = table(param_names', descriptions', units', defaults', ...
                       lower_bounds', upper_bounds', logical(optimized'), sources', ...
                       'VariableNames', {'Name', 'Description', 'Units', 'Default', ...
                                        'Lower', 'Upper', 'Optimized', 'Source'});

    %% =============================================================================
    %% 3. CREATE DEFAULT PARAMETER STRUCTURE
    %% =============================================================================

    default_params = struct();
    for i = 1:n_params
        default_params.(param_names{i}) = defaults(i);
    end

    %% =============================================================================
    %% 4. CREATE PARAMETER BOUNDS STRUCTURE
    %% =============================================================================

    param_bounds = struct();
    param_bounds.lower = struct();
    param_bounds.upper = struct();

    for i = 1:n_params
        param_bounds.lower.(param_names{i}) = lower_bounds(i);
        param_bounds.upper.(param_names{i}) = upper_bounds(i);
    end

    %% =============================================================================
    %% 5. DISPLAY PARAMETER TABLE
    %% =============================================================================

    fprintf('=============================================================\n');
    fprintf('GLUCOSAMINE QSP MODEL - PARAMETER DEFINITIONS\n');
    fprintf('=============================================================\n\n');

    fprintf('Total parameters: %d\n', n_params);
    fprintf('Fixed parameters: %d\n', sum(~optimized));
    fprintf('Optimized parameters: %d\n\n', sum(optimized));

    fprintf('Parameters to be optimized:\n');
    opt_idx = find(optimized);
    for i = 1:length(opt_idx)
        idx = opt_idx(i);
        fprintf('  %s: [%.3f, %.3f] (default: %.3f)\n', ...
                param_names{idx}, lower_bounds(idx), upper_bounds(idx), defaults(idx));
    end

    fprintf('\n=============================================================\n');

end

%% =============================================================================
%% HELPER FUNCTIONS
%% =============================================================================

function opt_param_names = get_optimized_params(param_info)
    % Get names of parameters to be optimized
    opt_idx = param_info.Optimized;
    opt_param_names = param_info.Name(opt_idx);
end

function bounds = get_parameter_bounds_array(param_info)
    % Get bounds as arrays for optimization
    opt_idx = param_info.Optimized;
    bounds.lower = param_info.Lower(opt_idx)';
    bounds.upper = param_info.Upper(opt_idx)';
    bounds.names = param_info.Name(opt_idx);
end
