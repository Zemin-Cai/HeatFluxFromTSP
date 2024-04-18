function [UnitHFMatrix] = UnitHFMatrixComputation(ThermalDiffus_p, ThermalDiffus_b,...
    ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UnitHFMatrixComputation: Compute the unit heat flux matrix based on the 
%                          fixed polymer layer and base layer to construct
%                          the coefficient matrix
%
% Syntax:  UnitHFMatrixComputation(...)
%
% Input arguments:
%  ThermalDiffus_p: the thermal diffusivity of the polymer layer
%  ThermalDiffus_b: the thermal diffusivity of the base layer
%  ThermalConduc_p: the thermal conductivity of the polymer layer
%  ThermalConduc_b: the thermal conductivity of the base layer
%  L: the thickness of the polymer layer
%  ymax_b: length of the space domain in the base layer (positive value)
%  tmax: maximum time for the simulation
%  ny_p: number of mesh points in y direction of polymer layer
%  ny_b: number of mesh points in y direction of base layer
%  nt: number of the step in time direction
%
% Output arguments:
%  UnitHFMatrix: the computed unit heat flux matrix. the size of this
%                matrix is (n-1)*(n-1), where n is the length of t vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% parameter setting

UnitHF_Func = '(t>=0).*(1) + (t<0).*(0)';
[theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
    ThermalDiffus_b, UnitHF_Func, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt);

Theata_ps_Unit = theata_p(end,:)';
n = size(Theata_ps_Unit,1);
phi = zeros(n, n);

for i = 1:n
    phi(i, i) = Theata_ps_Unit(i);
end
UnitHFMatrix = phi(2:n, 2:n);
