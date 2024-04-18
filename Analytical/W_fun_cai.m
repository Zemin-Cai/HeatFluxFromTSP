function [W, t] = W_fun_cai( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_polymer, L_base, t, h_c)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W_fun_cai: evaulate the W function
% 
% Syntax:  W_fun_cai(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThremalConduc_b: the thermal conductivity of the base layer
%   L_polymer and L_base: length of the space domain in the polymer layer
%    and base
%   Theata_ps: the temparature change in the polymer surface (a col vector)
%   t: the time vector
%
% Output arguments:
%   W: the value of the function W (a vector)
%   t: the time vector
%
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
%close all;
%addpath(GetAbsolutePath('../Numerical'));
%addpath(GetAbsolutePath('../Core funcs'));


% computation
epsilon = (ThermalConduc_p*sqrt(ThermalDiffus_b))/(ThermalConduc_b*sqrt(ThermalDiffus_p));
epsilon_ba = (1 - epsilon)/(1 + epsilon);

if(epsilon_ba > 0.9)
    %torrence = 1.0e-332;
    torrence = 1.0e-64;
else
    torrence = 1.0e-64;
end



dt = (t(end)-t(1))/(max(size(t))-1);
%%
% computation
epsilon = (ThermalConduc_p*sqrt(ThermalDiffus_b))/(ThermalConduc_b*sqrt(ThermalDiffus_p));
epsilon_ba = (1 - epsilon)/(1 + epsilon);

h_c_head = h_c*sqrt(ThermalDiffus_b)/ThermalConduc_b;
Delta = pi/480;

for k = 2:max(size(t))
    t_value = (k-1)*dt;
    
    % base on the derivation on 2017.06.10
    E_Func = ['exp(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_polymer) './sqrt(' num2str(ThermalDiffus_p) ').*cos((pi-' num2str(Delta) ')/2))'];
    m_Func = ['((cthe+' num2str(h_c_head) '.*i.*exp(i.*' num2str(Delta) '/2).*sqrt(' num2str(t_value) ')).*exp(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_base) './sqrt(' num2str(ThermalDiffus_b) ').*i.*exp(-i.*' num2str(Delta) '))./(cthe-' num2str(h_c_head) '.*i.*exp(i.*' num2str(Delta) './2).*sqrt(' num2str(t_value) ')))'];
    m_head_Func = ['((' num2str(epsilon_ba) '-' m_Func ')./(1-' num2str(epsilon_ba) '.*' m_Func '))'];
    A_Func = ['(abs(' m_head_Func ').*' E_Func ')'];
    Psai_Func = ['(angle(' m_head_Func '))'];
    Alpha_Func = ['(-cthe./sqrt(' num2str(t_value) ').*2.*' num2str(L_polymer) './sqrt(' num2str(ThermalDiffus_p) ').*sin((pi-' num2str(Delta) ')./2))'];
    Beta_Func = ['(' Psai_Func '+' Alpha_Func ')'];
    IntegralFuncStr = ['(1-' A_Func '.*' A_Func '+2.*' A_Func '.*sin(' Beta_Func ').*i).*exp(-cthe.^2.*exp(-i.*' num2str(Delta) ')-i.*' num2str(Delta) './2)./(1+' A_Func '.*' A_Func '-2.*' A_Func '.*cos(' Beta_Func '))'];
    %W_Fun(k-1) = 2/sqrt(pi)*real(IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence));
    W(k-1) = 2/sqrt(pi)*real(IntegralComputation_Infinite_CloseZero(IntegralFuncStr, 'cthe', torrence));
        
end



