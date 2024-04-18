function [q_s, t, Delta] = HF_Analy_FiniteBase_W1( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_polymer, L_base, Theata_ps, t, h_c, N_interp, W_Fun)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fuction calculates the heat flux with
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
%   N-interp: number for interpolation
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

Delta = pi/480;

tmax = max(max(t));
t_old=t;


N=N_interp;
dtq=tmax/N;
tq=[0:dtq:tmax];

Theata_ps_q=interp1(t_old,Theata_ps,tq);
W_Fun_q=interp1(t_old,W_Fun,tq);

Theata_ps=Theata_ps_q;
W_Fun=W_Fun_q;
t=tq;


% [W_Fun, t] = W_fun_nathan( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%     ThermalConduc_b, L_polymer, L_base, t, h_c);


q_s_interp(1) = 0;
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
    q_s_interp(k) = ThermalConduc_p/sqrt(ThermalDiffus_p*pi)*sum;
end
q_s_interp = q_s_interp';

q_s=interp1(t,q_s_interp,t_old); 
t=t_old;









