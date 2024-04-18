function [q_s, t] = HFComputation_FourierLaw( ThermalDiffus_p, ThermalConduc_p, L, Theata_ps, t)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_FourierLaw: calculate the heat flux directly using the 
%                           Fourier Law
% 
% Syntax:  HFComputation_FourierLaw(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   L: length of the space domain in the polymer layer
%   Theata_ps: the temparature change in the polymer surface (a col vector)
%   t: the time vector
%
% Output arguments:
%   q_s: the computed heat flux in the selected t point (a vector)
%   t: the time vector
%
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
% Zemin Cai 2008.07.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
FourierLawResult(1) = 0;
for k = 2:max(size(t))
    FourierLawResult(k) = ThermalConduc_p/L*Theata_ps(k);
end
q_s = FourierLawResult';
