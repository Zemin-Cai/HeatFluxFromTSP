function [q_s_matrix, theata_ps_matrix, t] = HFComputation_Numerical_FiniteBase_TemperatureDependent(ThermalDiffus_p, ThermalDiffus_b,...
    ThermalConduc_p, ThermalConduc_b, ThermalDiffus_p_Func, ThermalDiffus_b_Func, ThermalConduc_p_Func,...
        ThermalConduc_b_Func, Theata_ps, L, ymax_b, tmax, ny_p, ny_b, nt, tolerance, IterTimes, h_c, T_ini)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_Numerical_FiniteBase_TemperatureDependent: recover the heat flux from the inverse problem 
%             based on the numerical method for finite base
%             (the thermal conductivity and diffusivity are sensitive with the temperature change)
% 
% Syntax:  HFComputation_Numerical_VariableParam(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThremalConduc_b: the thermal conductivity of the base layer
%   ThermalDiffus_p_Func: the thermal diffusivity function of the polymer
%                         layer(independent variable is temperature 'T')
%   ThermalDiffus_b_Func: the thermal diffusivity function of the base layer
%   ThermalConduc_p_Func: the thermal conductivity function of the polymer layer
%   ThermalConduc_b_Func: the thermal conductivity function of the base layer
%   Theata_ps: the temparature change in the polymer surface (a col vector)
%   L: length of the space domain in the polymer layer
%   ymax_b: length of the space domain in the base layer (positive value)
%   tmax: maximum time for the simulation
%   ny_p: number of mesh points in y direction of polymer layer
%   ny_b: number of mesh points in y direction of base layer
%   nt: number of the step in time direction
%   tolerance: tolerable error
%   IterTimes: the iterational times
%   T_ini: initial temperature in both the polymer layer and vase (in K)
%
% Output arguments:
%   q_s_matrix: the computed heat flux (a matrix, the ith column is the recovered result in the ith iteration)
%   theata_ps_matrix: the recomputed heat flux from the computed heat flux
%                     in each iteration
%   t: the time vector
%
% Zemin Cai 2008.09.11
%  Modified on 2017.08.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 17
    error('Not enough input arguments!!!');
end
%%
% compute the unit heat flux response matrix
UnitHFMatrix = UnitHFMatrixComputation_FiniteBase_TemperatureDependent(ThermalDiffus_p, ThermalDiffus_b,...
    ThermalConduc_p, ThermalConduc_b, ThermalDiffus_p_Func, ThermalDiffus_b_Func, ThermalConduc_p_Func,...
    ThermalConduc_b_Func, L, ymax_b, tmax, ny_p, ny_b, nt, h_c, T_ini);
UnitHFMatrix_inv = inv(UnitHFMatrix);

%% Recover the heat flux
q_s_matrix = [];
k = 1;
disp(sprintf('Iter %d', k));
HF = UnitHFMatrix_inv*Theata_ps(2:end);          % the initial heat flux
HF = [0;HF];
q_s_matrix = [q_s_matrix HF];

% Recompute the surface temperature using the recovered heat flux
theata_ps_matrix = [];
[theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS_FiniteBase_TemperatureDependent( ThermalDiffus_p,...
        ThermalDiffus_b, ThermalConduc_p, ThermalConduc_b, ThermalDiffus_p_Func, ThermalDiffus_b_Func,...
        HF, ThermalConduc_p_Func, ThermalConduc_b_Func, L, ymax_b, tmax, ny_p, ny_b, nt, h_c, T_ini);
theata_ps_recomputed = theata_p(end, :)';
theata_ps_matrix = [theata_ps_matrix theata_ps_recomputed];
Terror = Theata_ps - theata_ps_recomputed;

% iteration
while(k<IterTimes & norm(Terror(floor(end/5):end))>tolerance)
%while(norm(Terror(floor(end/5):end))>tolerance)
    HFError = UnitHFMatrix_inv*Terror(2:end);
    HFError = [0;HFError];
    HF = HF + HFError;
    q_s_matrix = [q_s_matrix HF];
    
    [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS_FiniteBase_TemperatureDependent( ThermalDiffus_p,...
        ThermalDiffus_b, ThermalConduc_p, ThermalConduc_b, ThermalDiffus_p_Func, ThermalDiffus_b_Func,...
        HF, ThermalConduc_p_Func, ThermalConduc_b_Func, L, ymax_b, tmax, ny_p, ny_b, nt, h_c, T_ini);
    theata_ps_recomputed = theata_p(end, :)';
    theata_ps_matrix = [theata_ps_matrix theata_ps_recomputed];
    
    k = k + 1;
    Terror = Theata_ps - theata_ps_recomputed;
    disp(sprintf('Iter %d', k));
end