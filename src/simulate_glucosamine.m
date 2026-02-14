function results = simulate_glucosamine(params, dose_mg, duration_weeks, formulation, treatment)
% SIMULATE_GLUCOSAMINE Simulate glucosamine treatment effects on knee OA
%
% Core simulation function for the Glucosamine QSP Model.
%
% Drug effect mechanisms:
%   1. Degradation inhibition (Imax_deg, IC50_deg) - primary anti-catabolic
%   2. Synthesis stimulation (Emax_syn, EC50_syn) - minor anabolic
%
% Inputs:
%   params          - Structure with model parameters (from M01_model_structure)
%   dose_mg         - Daily dose in mg (default: 1500)
%   duration_weeks  - Treatment duration in weeks (default: 156 = 3 years)
%   formulation     - 'sulfate' or 'HCl' (default: 'sulfate')
%   treatment       - Logical, true for treatment, false for placebo
%
% Outputs:
%   results - Table with columns: week, time_years, GAG, JSW, JSW_change,
%             Pain, Pain_change
%
% Note: R_syn, Vd_apparent, and gamma_JSW are defined in the parameter
% structure (M01) for documentation but are not used in this implementation.
% The PK uses a reference-scaling approach (see Section 2.7 of manuscript).
% gamma_JSW = 1.0 is equivalent to the linear JSW = JSW_0 * (GAG/GAG_0).
%
% Example:
%   [~, params, ~] = M01_model_structure();
%   sim = simulate_glucosamine(params, 1500, 156, 'sulfate', true);

if nargin < 2, dose_mg = 1500; end
if nargin < 3, duration_weeks = 156; end
if nargin < 4, formulation = 'sulfate'; end
if nargin < 5, treatment = true; end

% -------------------------------------------------------------------------
% 1. Pharmacokinetics
% -------------------------------------------------------------------------

if strcmpi(formulation, 'sulfate')
    F_bio = params.F_sulfate;
elseif strcmpi(formulation, 'HCl')
    F_bio = params.F_HCl;
else
    F_bio = params.F_sulfate;
end

if ~treatment
    F_bio = 0;
end

% Effective concentration at target site
C_eff_reference = 0.30;   % ug/mL at target
F_reference = 0.22;
dose_reference = 1500;

C_eff = C_eff_reference * (F_bio / F_reference) * (dose_mg / dose_reference);

% -------------------------------------------------------------------------
% 2. Drug Effects
% -------------------------------------------------------------------------

% Synthesis stimulation (minor anabolic effect)
if isfield(params, 'Emax_syn')
    Emax_syn = params.Emax_syn;
else
    Emax_syn = 0.15;
end

if isfield(params, 'EC50_syn')
    EC50_syn = params.EC50_syn;
else
    EC50_syn = 2.0;
end

if C_eff > 0
    E_syn = 1 + Emax_syn * C_eff / (EC50_syn + C_eff);
else
    E_syn = 1.0;
end

% Degradation inhibition (primary anti-catabolic effect)
if C_eff > 0
    I_drug = params.Imax_deg * C_eff / (params.IC50_deg + C_eff);
    I_drug = max(0, min(0.99, I_drug));
else
    I_drug = 0;
end

% -------------------------------------------------------------------------
% 3. Initialize state vectors
% -------------------------------------------------------------------------

weeks = (0:duration_weeks)';
n = length(weeks);

GAG = zeros(n, 1);
JSW = zeros(n, 1);
Pain = zeros(n, 1);

GAG(1) = params.GAG_0;
JSW(1) = params.JSW_0;

% Convert annual rates to weekly
k_syn_week = params.k_syn / 52;
k_deg_week = params.k_deg / 52;

% -------------------------------------------------------------------------
% 4. GAG Dynamics (weekly Euler integration)
% -------------------------------------------------------------------------

for i = 2:n
    % dGAG/dt = k_syn * E_syn - k_deg * (1 - I_drug) * GAG
    effective_k_syn = k_syn_week * E_syn;
    effective_k_deg = k_deg_week * (1 - I_drug);

    dGAG = effective_k_syn - effective_k_deg * GAG(i-1);
    GAG(i) = GAG(i-1) + dGAG;

    % Physiological bounds
    GAG(i) = max(0.50, min(1.05, GAG(i)));
end

% -------------------------------------------------------------------------
% 5. Clinical Endpoints
% -------------------------------------------------------------------------

for i = 1:n
    % JSW from GAG: linear relationship (gamma_JSW = 1.0)
    JSW(i) = params.JSW_0 * (GAG(i) / params.GAG_0);

    % Pain score
    time_days = weeks(i) * 7;

    % Structural pain from cartilage loss
    GAG_loss = max(0, params.GAG_0 - GAG(i));
    pain_structural = params.alpha_struct * GAG_loss;

    % Placebo effect (exponential onset, no waning)
    placebo_effect = params.Pmax * (1 - exp(-params.kpl * time_days));

    % Direct drug effect on pain
    if C_eff > 0
        drug_pain_effect = params.beta_drug * C_eff / (params.EC50_pain + C_eff);
    else
        drug_pain_effect = 0;
    end

    Pain(i) = params.Pain_0 + pain_structural - placebo_effect - drug_pain_effect;
    Pain(i) = max(0, min(100, Pain(i)));
end

% -------------------------------------------------------------------------
% 6. Return results
% -------------------------------------------------------------------------

time_years = weeks / 52;
JSW_change = JSW - JSW(1);
Pain_change = Pain - Pain(1);

results = table(weeks, time_years, GAG, JSW, JSW_change, Pain, Pain_change, ...
    'VariableNames', {'week', 'time_years', 'GAG', 'JSW', 'JSW_change', 'Pain', 'Pain_change'});

end
