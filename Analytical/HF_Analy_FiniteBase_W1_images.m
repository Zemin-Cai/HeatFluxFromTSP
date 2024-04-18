function [q_s, t, Delta] = HF_Analy_FiniteBase_W1_images( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_polymer, L_base, Theata_ps, t, h_c, W_Fun)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This fuction calculates the heat flux by using the analytical method for
% the finite-thickness base.  Summation is done directly for a sequence of images.  
% 
% Syntax:  [q_s, t, Delta] = HF_Analy_FiniteBase_W1_images(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThremalConduc_b: the thermal conductivity of the base layer
%   L_polymer and L_base: length of the space domain in the polymer layer
%    and base
%   Theata_ps: the temparature change in the polymer surface (3D array)
%   t: the time vector
%   W_Fun: W function (3D array)
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

[n1,n2,n3]=size(Theata_ps);


for k =2:max(size(t))
    sum = zeros(n1,n2);
    for j = 2:k
        t_value1 = t(k) - t(j);
        t_value2 = t(k) - t(j-1);
        if(t_value1 ~= 0)
            W_Value1 = W_Fun(:,:,k-j);
        else
            W_Value1 = zeros(n1,n2);
        end
        if(t_value2 ~= 0)
            W_Value2 = W_Fun(:,:,k-j+1);
        else
            W_Value2 = zeros(n1,n2);
        end
        sum = sum + (Theata_ps(:,:,j) - Theata_ps(:,:,j-1)).*(W_Value1 + W_Value2)./(sqrt(t_value1)+sqrt(t_value2));
    end
    q_s_interp(:,:,k) = ThermalConduc_p/sqrt(ThermalDiffus_p*pi)*sum;
end

q_s=q_s_interp;









