
%% This program calculates the relevant paramters of the impinging jet.  


clear all
close all

%T_off=71.3+273.15; % case 2, initial temperature (K)
T_off=43.4+273.15; % case 3
%T_off=33.8+273.15; % case 4


ThermalDiffus_p = 9.8*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.25;                 % W/m-K

ThermalDiffus_b = 8.3579*10^(-5);            % Al base
ThermalConduc_b = 204;

% ThermalDiffus_b = 4.051*10^(-6);            % Stainless Steel base
% ThermalConduc_b = 16;

% ThermalDiffus_b = 1.17*10^(-5);            % Steel base
% ThermalConduc_b = 81;


L_p = 40e-6; % m
  
T_total=23+273.15; % Total temperature (K)
p_total=25*6894.757; % Pa, total pressure
R=287; % J/kg-K, universal gas constant for air
rho_total=p_total/(R*T_total); % total density (kg/m^2)

% ambient conditions
p_amb=101325; % Pa
T_amb=23+273.15; % K
p_ratio=p_total/p_amb;

p_amb_psi=101325/6894.757;

% Conditions at jet exit
M_sonic=1; % Mach No. at the sonic condition
T_sonic=T_total*0.833; % temperature (K)
p_sonic=p_total*0.528; % pressure (Pa)
rho_sonic=rho_total*0.634; % density (kg/m^2)

gamma=1.4;
a_sonic=(gamma*R*T_sonic)^0.5; % speed of sound (m/s) at the sonic condition
u_sonic=M_sonic*a_sonic; % velocity at the sonic condition

% Sutherland formula for air viscosity
mu0=18.268*10^(-6); % kg/m-s
T0=291.15; % refrence temparature (K)
C=120; %K

mu_sonic=mu0*((T0+C)/(T_sonic+C))*(T_sonic/T0)^(3/2); % dynamic viscosity of air (kg/m-s)

D=6.35*10^(-3); % exit diameter (m)
Re_D=u_sonic*D*rho_sonic/mu_sonic;






