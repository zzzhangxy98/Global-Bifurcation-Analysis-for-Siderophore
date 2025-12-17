function F = ODE_siderophore_iron_partition(X, opt)
% ODE_SIDEROPHORE_IRON_PARTITION Saddle-node dynamics compatible version
% 
% Inputs:
%   X:   State matrix (N_spe + N_sid + 1) × n_points
%        - X(1:N_spe, :):             Biomass concentrations M_i
%        - X(N_spe+1:N_spe+N_sid, :): Siderophore concentrations R_j
%        - X(N_spe+N_sid+1, :):       Free iron concentration R_iron
%   opt: Structure containing all parameters:
%        - alpha_sid_pro: (N_spe × N_sid) allocation matrix for siderophore production
%        - v_sid_rec:     (N_spe × N_sid) receptor allocation matrix
%        - e:             Siderophore production rate (scalar or N_sid × 1 vector)
%        - u:             Iron uptake rate (scalar or N_sid × 1 vector)
%        - migr:          Migration rate (scalar)
%        - gamma:         Growth conversion efficiency (scalar)
%        - d:             Dilution/death rate (scalar)
%        - R_sup:         Iron supply concentration (scalar)
%
% Output:
%   F:   Derivative matrix (same size as X)
%
% Example usage:
%   opt.alpha_sid_pro = [0.3,0,0; 0,0.4,0; 0,0,0.5];
%   opt.v_sid_rec = [0.3,0.7,0; 0,0.3,0.7; 0.7,0,0.3];
%   opt.e = 10;
%   opt.u = 1;
%   opt.migr = 0;
%   opt.gamma = 1;
%   opt.d = 0.1;
%   opt.R_sup = 1;
%   F = ODE_siderophore_iron_partition(X, opt);

%% Extract parameters from opt structure
alpha_sid_pro = opt.alpha_sid_pro;  % N_spe × N_sid matrix
v_sid_rec = opt.v_sid_rec;          % N_spe × N_sid matrix

% Set default values if fields don't exist
if isfield(opt, 'e')
    e = opt.e;
else
    e = 10;
end

if isfield(opt, 'u')
    u = opt.u;
else
    u = 10;
end

if isfield(opt, 'migr')
    migr = opt.migr;
else
    migr = 1e-8;
end

if isfield(opt, 'gamma')
    gamma = opt.gamma;
else
    gamma = 1;
end

if isfield(opt, 'd')
    d = opt.d0;
else
    d = 0.1;
end

if isfield(opt, 'R_sup')
    R_sup = opt.R_sup;
else
    R_sup = 1;
end

%% Determine dimensions
N_spe = size(v_sid_rec, 1);  % Number of species
N_sid = size(v_sid_rec, 2);  % Number of siderophore types
n_points = size(X, 2);        % Number of state points (for parallel computation)

%% Compute growth allocation
% alpha_i0 = 1 - sum_j(alpha_ij)
alpha_growth = 1 - sum(alpha_sid_pro, 2);  % N_spe × 1 vector

%% Convert scalar parameters to vectors if needed
if isscalar(e)
    e = e * ones(N_sid, 1);
end
if isscalar(u)
    u = u * ones(N_sid, 1);
end

%% Extract state variables from X
M = X(1:N_spe, :);                      % Biomass: N_spe × n_points
R = X(N_spe+1:N_spe+N_sid, :);         % Siderophores: N_sid × n_points
Fe = X(N_spe+N_sid+1, :);              % Free iron: 1 × n_points

%% Compute iron uptake rates
% J_j = u_j * R_j * Fe (Equation 6, mass action kinetics)
J = (u .* R) .* Fe;  % Element-wise: N_sid × n_points

%% Compute biomass derivatives (Equation 1)
% dM_i/dt = migr + (gamma * alpha_i0 * sum_j(v_ij * J_j) - d) * M_i
iron_uptake_per_species = v_sid_rec * J;  % N_spe × n_points
growth_rate = gamma * alpha_growth .* iron_uptake_per_species - d;  % N_spe × n_points
dMdt = migr + growth_rate .* M;  % N_spe × n_points

%% Compute siderophore derivatives (Equation 2)
% dR_j/dt = -d * R_j + epsilon_j * sum_i(alpha_ij * M_i)
siderophore_production = alpha_sid_pro' * M;  % N_sid × n_points
dRdt = -d * R + e .* siderophore_production;  % N_sid × n_points

%% Compute free iron derivative (Equation 3)
% dFe/dt = d * (R_sup - Fe) - sum_i(M_i * sum_j(v_ij * J_j))
iron_consumption = sum(M .* iron_uptake_per_species, 1);  % 1 × n_points
dFedt = d * (R_sup - Fe) - iron_consumption;  % 1 × n_points

%% Assemble output matrix
F = [dMdt; dRdt; dFedt];  % (N_spe + N_sid + 1) × n_points

end