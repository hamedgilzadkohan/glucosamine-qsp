%% =============================================================================
%% M02_model_equations.m
%% Model Equations Test and Demonstration
%% =============================================================================
%%
%% Glucosamine QSP Model for Knee Osteoarthritis
%%
%% This script runs a verification simulation using the core
%% simulate_glucosamine() function to confirm the model equations produce
%% expected outputs.
%%
%% Dependencies:
%%   - M01_model_structure.m (parameter definitions)
%%   - simulate_glucosamine.m (core simulation function)
%%
%% =============================================================================

fprintf('=============================================================\n');
fprintf('GLUCOSAMINE QSP MODEL - EQUATIONS TEST\n');
fprintf('=============================================================\n\n');

% Load default parameters
[param_info, default_params, param_bounds] = M01_model_structure();

fprintf('Running test simulation (3 years, 1500 mg/day sulfate)...\n\n');

sim_p = simulate_glucosamine(default_params, 1500, 156, 'sulfate', false);
sim_t = simulate_glucosamine(default_params, 1500, 156, 'sulfate', true);

fprintf('Results at Year 3 (Week 156):\n');
fprintf('  Placebo JSW change: %.3f mm\n', sim_p.JSW_change(end));
fprintf('  Treatment JSW change: %.3f mm\n', sim_t.JSW_change(end));
fprintf('  Treatment effect: %.3f mm\n', sim_t.JSW_change(end) - sim_p.JSW_change(end));
fprintf('  Placebo pain: %.1f\n', sim_p.Pain(end));
fprintf('  Treatment pain: %.1f\n', sim_t.Pain(end));

fprintf('\n=============================================================\n');
fprintf('Model equations test complete.\n');
fprintf('=============================================================\n');
