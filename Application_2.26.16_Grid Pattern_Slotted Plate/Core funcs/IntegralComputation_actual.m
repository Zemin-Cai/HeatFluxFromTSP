%%
function IntegralResult = IntegralComputation_actual(t_value, h_c_head, ThermalDiffus_p, ThermalDiffus_b, L_base, L_polymer, epsilon_ba)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegralComputation with int function
%2017.06.07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
syms cthe;
m_func = (cthe*i - h_c_head.*sqrt(t_value)).*exp(-i*sqrt(ThermalDiffus_p./ThermalDiffus_b).*L_base./L_polymer)./(cthe*i + h_c_head.*sqrt(t_value));
m_head_func = (epsilon_ba - m_func)./(1 - epsilon_ba.*m_func);
A_func = abs(m_head_func);
Beta_func = 2.*L_polymer.*cthe./sqrt(ThermalDiffus_p.*t_value) - angle(m_head_func);
f = (1 - A_func.^2).*exp(-cthe.^2)./(1 + A_func.^2 - 2.*A_func.*cos(Beta_func));
IntegralResult = int(f, cthe, 0, inf);

