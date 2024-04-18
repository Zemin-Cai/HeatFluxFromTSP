function [q_s, t, epsilon_ba] = HF_Analy_SemiInfBase( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L, Theata_ps, t)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the heat flux directly for the general 
%                            case using analytical method
% 
% Syntax:  HFComputation_Analytical(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThremalConduc_b: the thermal conductivity of the base layer
%   L: length of the space domain in the polymer layer
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

dt = (t(end)-t(1))/(max(size(t))-1);
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

% if dt is fixed (constant), we can have fast computation
for k = 2:max(size(t))
    t_value = (k-1)*dt;
    IntegralFuncStr = strcat(['(exp(-cthe.^2)./(1+', num2str(epsilon_ba), '^2-2*', num2str(epsilon_ba), '*cos(2*',...
           num2str(L), './sqrt(', num2str(ThermalDiffus_p*t_value), ').*cthe)))']);
%     IntegralFuncStr = strcat(['(exp(-' num2str(ThermalDiffus_p*t_value/4*L^2) '*cthe.^2)./(1+', num2str(epsilon_ba), '^2-2*', num2str(epsilon_ba), '*cos(cthe)))'])
    Whead_Value(k-1) = 2/sqrt(pi)*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
    %Whead_Value(k-1) = (1-epsilon_ba^2)/pi*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
    %Whead_Value(k-1) = 2*L*(1-epsilon_ba^2)/(pi*sqrt(ThermalDiffus_p*t_value))*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
end


%Whead_Value'
q_s(1) = 0;
for k = 2:max(size(t))
    sum = 0;
    for i = 2:k
        t_value1 = t(k) - t(i);
        if(t_value1 ~= 0)
            Whead_Value1 = Whead_Value(k-i);
        else
            Whead_Value1 = 0;
        end
        t_value2 = t(k) - t(i-1);
        if(t_value2 ~= 0)
            Whead_Value2 = Whead_Value(k-i+1);
        else
            Whead_Value2 = 0;
        end
        
        % sum
        %a=Theata_ps(i) - Theata_ps(i-1)
        %b=Whead_Value1 + Whead_Value2
        %c=sqrt(t_value1) + sqrt(t_value2)
        sum = sum + (Theata_ps(i) - Theata_ps(i-1))*(Whead_Value1 + Whead_Value2)/(sqrt(t_value1) + sqrt(t_value2));
        %(Theata_ps(i) - Theata_ps(i-1))*(Whead_Value1 + Whead_Value2)/2;
        %sum = sum + (Theata_ps(i) - Theata_ps(i-1))*(Whead_Value1 + Whead_Value2)/2;
    end
    %sum
    %ThermalConduc_p*(1 - epsilon_ba^2)/sqrt(pi*ThermalDiffus_p)
    q_s(k) = ThermalConduc_p*(1 - epsilon_ba^2)/sqrt(pi*ThermalDiffus_p)*sum;
    %q_s(k) = ThermalConduc_p/L*sum;
end
q_s = q_s';


