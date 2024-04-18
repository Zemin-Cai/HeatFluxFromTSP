function [theata_ps_Finite_analytical, t_temperature_analytical, Delta] = SurfaceTemperatureComputation_Analytical_FiniteBase(ThermalDiffus_p,...
    ThermalConduc_p, ThermalDiffus_b, ThermalConduc_b, Heatflux_Original, L_polymer, L_base, t, h_c)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SurfaceTemperatureComputation_Analytical_FiniteBase: 
%    Surface temperature calculation from heat flux for the two layers conduction
%    problem of finite-thick base using our analytical solution
% 
% Syntax:  SurfaceTemperatureComputation_Analytical_FiniteBase(ThermalDiffus_p,...
%    ThermalConduc_p, ThermalDiffus_b, ThermalConduc_b, Heatflux_Original, L_polymer, L_base, tmax, ny_p, ny_b, nt, h_c);
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThermalConduc_b: the thermal conductivity of the base layer
%   Heatflux_Original: the heat flux used, a vector
%   L_polymer: length of the space domain in the polymer layer
%   L_base: thickness of the base
%   tmax: maximum time for the simulation
%   ny_p: number of mesh points in y direction of polymer layer
%   ny_b: number of mesh points in y direction of base layer
%   nt: number of the step in time direction
%   h_c: an empirical value
%
% Output arguments:
%   theata_ps_Finite_analytical: the calculation result of the surface temperature
%   t_temperature_analytical: values of time at which solution is obtained
% 
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
% Zemin Cai 2017.06.13
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
    %HeatfluxFunc = '0.7 + 0.3*sin(pi*t)';      % a function of t
    HeatfluxFunc = '(t>0).*(25000) + (t<=0).*(0)';
    %HeatfluxFunc = '34000*(1-exp(-5*t)+18000*t+8000*t.^2)'; 
                                                % use the analytical result
                                                % by Prof.Liu in 1994's
                                                % paper
end
if nargin < 4
    ThermalConduc_p = 0.15;
end
if nargin < 5
    ThermalConduc_b = 204;
end
if nargin < 6
    L_p = 5.0e-5;
end
if nargin < 7
    L_b = 0.002;
end
if nargin < 8
    tmax = 6;
end
if nargin < 9
    ny_p = 20;
end
if nargin < 10
    ny_b = 500;
end
if nargin < 11
    nt = 80;
end
if nargin < 12
    errPlots = 1;
end
if nargin < 13
    h_c = 0;
end

%%
if max(size(t)) ~= max(size(Heatflux_Original))
    error('the size of heat flux vector and t vector should match!');
end

dt = (t(end)-t(1))/(max(size(t))-1);
%%
% computation
epsilon = (ThermalConduc_p*sqrt(ThermalDiffus_b))/(ThermalConduc_b*sqrt(ThermalDiffus_p));
epsilon_ba = (1 - epsilon)/(1 + epsilon)

if(epsilon_ba > 0.9)
    %torrence = 1.0e-332;
    torrence = 1.0e-64;
else
    torrence = 1.0e-64;
end

h_c_head = h_c*sqrt(ThermalDiffus_b)/ThermalConduc_b
%Delta = pi/4;
%Delta = pi/512;
Delta = pi/1024;
%Delta = pi/2048;
%Delta = pi/4096;

ComputationalTime = cputime;
for k = 2:max(size(t))
    t_value = (k-1)*dt;
    
    % base on the derivation on 2017.06.12
    E_Func = ['exp(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_polymer) './sqrt(' num2str(ThermalDiffus_p) ').*cos((pi-' num2str(Delta) ')/2))'];
    m_Func = ['((cthe+' num2str(h_c_head) '.*i.*exp(i.*' num2str(Delta) '/2).*sqrt(' num2str(t_value) ')).*exp(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_base) './sqrt(' num2str(ThermalDiffus_b) ').*i.*exp(-i.*' num2str(Delta) '))./(cthe-' num2str(h_c_head) '.*i.*exp(i.*' num2str(Delta) './2).*sqrt(' num2str(t_value) ')))'];
    m_head_Func = ['((' num2str(epsilon_ba) '-' m_Func ')./(1-' num2str(epsilon_ba) '.*' m_Func '))'];
    A_Func = ['(abs(' m_head_Func ').*' E_Func ')'];
    Psai_Func = ['(angle(' m_head_Func '))'];
    Alpha_Func = ['(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_polymer) './sqrt(' num2str(ThermalDiffus_p) ').*sin((pi-' num2str(Delta) ')./2))'];
    Beta_Func = ['(' Psai_Func '+' Alpha_Func ')'];
    IntegralFuncStr = ['(1-' A_Func '.*' A_Func '-2.*' A_Func '.*sin(' Beta_Func ').*i).*exp(-cthe.^2.*exp(-i.*' num2str(Delta) ')-i.*' num2str(Delta) './2)./(1+' A_Func '.*' A_Func '+2.*' A_Func '.*cos(' Beta_Func '))'];
    U_Fun(k-1) = 2/sqrt(pi)*real(IntegralComputation_Infinite_CloseZero(IntegralFuncStr, 'cthe', torrence));
   
end

theata_ps_Finite_analytical(1) = 0;
for k =2:max(size(t))
    sum = 0;
    for j = 2:k
        t_value1 = t(k) - t(j);
        t_value2 = t(k) - t(j-1);
        if(t_value1 ~= 0)
            U_Value1 = U_Fun(k-j);
        else
            U_Value1 = 0;
        end
        if(t_value2 ~= 0)
            U_Value2 = U_Fun(k-j+1);
        else
            U_Value2 = 0;
        end
        sum = sum + (Heatflux_Original(j) + Heatflux_Original(j-1))*(U_Value1 + U_Value2)./(sqrt(t_value1)+sqrt(t_value2));
        %sum = sum + (Heatflux_Original(j))*(U_Value1 + U_Value2)./(sqrt(t_value1)+sqrt(t_value2));
        %sum = sum + (Heatflux_Original(j) - Heatflux_Original(j-1))*(U_Value1 + U_Value2)./(sqrt(t_value1)+sqrt(t_value2));
    end
    theata_ps_Finite_analytical(k) = sqrt(ThermalDiffus_p)/(ThermalConduc_p*sqrt(pi))*sum;
end
theata_ps_Finite_analytical = theata_ps_Finite_analytical';
t_temperature_analytical = t;
ComputationalTime = cputime - ComputationalTime;
