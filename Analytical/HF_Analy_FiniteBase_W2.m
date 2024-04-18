function [q_s, t, epsilon_ba, Delta] = HFComputation_Analytical_FiniteBase_W2( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_polymer, L_base, Theata_ps, t, h_c, W_Fun)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the heat flux with
% finite-thickness base (Fast convergence version)
% 
% Syntax:  HFComputation_Analytical_FiniteBase_FastConvergence(...)
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
%   q_s: the computed heat flux in the selected t point (a vector)
%   t: the time vector
%   epsilon_ba: to verify its cagegory
%
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
if max(size(t)) ~= max(size(Theata_ps))
    error('the size of Theata_ps vector and t vector should match!');
end


%%
% computation
epsilon = (ThermalConduc_p*sqrt(ThermalDiffus_b))/(ThermalConduc_b*sqrt(ThermalDiffus_p));
epsilon_ba = (1 - epsilon)/(1 + epsilon);

if(epsilon_ba > 0.9)
    %torrence = 1.0e-332;
    torrence = 1.0e-64;
else
    torrence = 1.0e-64;
end

%h_c = 250;                         % convective heat transfer coefficient
%h_c = 0;
h_c_head = h_c*sqrt(ThermalDiffus_b)/ThermalConduc_b;

Delta = pi/480;
% [W_Fun, t] = W_fun_cai( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%     ThermalConduc_b, L_polymer, L_base, t, h_c);

q_s(1) = 0;
for k =2:max(size(t))
    sum = 0;
    for j = 2:k
        t_value1 = t(k) - t(j);
        t_value2 = t(k) - t(j-1);
        if(t_value1 ~= 0)
            W_Value1 = W_Fun(k-j);
        else
            W_Value1 = 0;
        end
        if(t_value2 ~= 0)
            W_Value2 = W_Fun(k-j+1);
        else
            W_Value2 = 0;
        end
        sum = sum + (Theata_ps(j) - Theata_ps(j-1))*(W_Value1 + W_Value2)./(sqrt(t_value1)+sqrt(t_value2));
    end
    q_s(k) = ThermalConduc_p/sqrt(ThermalDiffus_p*pi)*sum;
end
q_s = q_s';










