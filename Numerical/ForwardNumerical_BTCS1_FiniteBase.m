function [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS1( ThermalDiffus_p,...
    ThermalDiffus_b, HeatfluxVec, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt, h_c)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ForwardNumerical_BTCS: Using the Backward Time, Centered Space Scheme to solve 
%                        the forward heat transfer problem (the heat flux value are given)
% 
% Syntax:  ForwardNumerical_BTCS(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   HeatfluxVec: the heat flux vector (a vector which saved the value)
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalConduc_b: the thermal conductivity of the base layer
%   L: length of the space domain in the polymer layer
%   ymax_b: length of the space domain in the base layer (positive value)
%   tmax: maximum time for the simulation
%   ny_p: number of mesh points in y direction of polymer layer
%   ny_b: number of mesh points in y direction of base layer
%   nt: number of the step in time direction
%
% Output arguments:
%   theata_p: matrix of solutions of the polymer layer. theata_p(:,k) is
%             theata_p(y) at t = t(k)
%   y_p: location of finite difference nodes of the polymer layer
%   err_p: L2 norm of error evaluated at the spatial nodes on last time
%          step of the polymer layer
%   theata_b: matrix of solutions of the base layer. theata_b(:,k) is
%             theata_b(y) at t = t(k)
%   y_b: location of finite difference nodes of the base layer
%   err_b: L2 norm of error evaluated at the spatial nodes on last time
%          step of the base layer
%   t: values of time at which solution is obtained
% 
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
% Reference:
%   T. Liu et.al., "Remote Surface Temperature and Heat Transfer Mapping for
%   a Waverider Model at Mach 10 Using Fluorescent Paint", 18th AIAA 
%   Aerospace Ground Testing Conference, June 20-23, 1994.
%
% Zemin Cai 2008.05.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
%close all;
%% parameters initialize
if nargin < 1
    ThermalDiffus_p = 9.7e-8;    % assume that it is a insulate layer whose 
                                 % Thermal Conductivity is 0.15, specific
                                 % heat is 1090 while the density is 1420
end
if nargin < 2
    ThermalDiffus_b = 8.36e-5;   % assume that it is aluminum with thermal 
                                 % conductivity 204, specific heat 904 and 
                                 % density 2700     
end
if nargin < 3
    HeatFluxVec = linspace(3000, 3000, 501);
end
if nargin < 4
    ThermalConduc_p = 0.15;
end
if nargin < 5
    ThermalConduc_b = 204;
end
if nargin < 6
    L = 50e-6;
end
if nargin < 7
    ymax_b = sqrt(204*0.01*500^2/19.6);
end
if nargin < 8
    tmax = 5;
end
if nargin < 9
    ny_p = 20;
end
if nargin < 10
    ny_b = 500;
end
if nargin < 11
    nt = 501;
end

% mesh dissection 
dy_p = L/(ny_p - 1);
dy_b = ymax_b/(ny_b-1);
dt = tmax/(nt - 1);


% initialization
y_p = linspace(0, L, ny_p)';
theata_p = zeros(ny_p, nt);
err_p = 0;
r_p = ThermalDiffus_p*dt/dy_p.^2;
%r_p = vpa(ThermalDiffus_p*dt/dy_p.^2);
r2_p = 1 + 2*r_p;

y_b = linspace(0, ymax_b, ny_b)';
theata_b = zeros(ny_b, nt);
err_b = 0;
r_b = ThermalDiffus_b*dt/dy_b.^2;
%r_b = vpa(ThermalDiffus_b*dt/dy_b.^2);
r2_b = 1 + 2*r_b;

t = linspace(0, tmax, nt)';
HeatFlux = HeatfluxVec;

theata_p(:, 1) = 0;
theata_b(:, 1) = 0;
theata_b(ny_b, :) = 0;

%Omiga = vpa((ThermalConduc_b*dy_p)/(ThermalConduc_p*dy_b));
Omiga = (ThermalConduc_b*dy_p)/(ThermalConduc_p*dy_b);

%% computation
% Coefficient matrix construction
Coef_matrix = zeros((ny_p+ny_b), (ny_p+ny_b));
Coef_matrix(1, 1) = 1 + Omiga;         % first row
Coef_matrix(1, 2) = -1;
Coef_matrix(1, ny_p+2) = -Omiga;
for i = 2:(ny_p-1)
    Coef_matrix(i, i-1) = -r_p;
    Coef_matrix(i, i) = r2_p;
    Coef_matrix(i, i+1) = -r_p;
end
Coef_matrix(ny_p, ny_p-1) = -1;
Coef_matrix(ny_p, ny_p) = 1;
Coef_matrix(ny_p+1, 2) = -1;
Coef_matrix(ny_p+1, ny_p+1) = 1 + Omiga;
Coef_matrix(ny_p+1, ny_p+2) = -Omiga;
for i = (ny_p+2):(ny_p+ny_b-1)
    Coef_matrix(i, i-1) = -r_b;
    Coef_matrix(i, i) = r2_b;
    Coef_matrix(i, i+1) = -r_b;
end
Coef_matrix(ny_p+ny_b, ny_p+ny_b-1) = -ThermalConduc_b/dy_b;
Coef_matrix(ny_p+ny_b, ny_p+ny_b) = h_c + ThermalConduc_b/dy_b;

for k = 2:nt
    % update the right hand side
    b = zeros(ny_p+ny_b, 1);
    for i = 2:(ny_p-1)
        b(i) = theata_p(i,k-1);
    end
    b(ny_p) = dy_p*HeatFlux(k)/ThermalConduc_p;
    %b(ny_p) = vpa(dy_p*HeatFlux(k)/ThermalConduc_p);
    for i = (ny_p+2):(ny_p+ny_b-1)
        b(i) = theata_b(i-ny_p,k-1);
    end
    
    theata = linsolve(Coef_matrix, b);
    theata_p(:, k) = theata(1:ny_p);
    theata_b(:, k) = theata((ny_p+1):end);
end
