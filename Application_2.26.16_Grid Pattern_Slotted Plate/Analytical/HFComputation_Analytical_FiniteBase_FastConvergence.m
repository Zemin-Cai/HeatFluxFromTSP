function [q_s, t, epsilon_ba, Delta] = HFComputation_Analytical_FiniteBase_FastConvergence( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_polymer, L_base, Theata_ps, t, h_c)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_Analytical_FiniteBase_FastConvergence: calculate the heat flux with
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
% Zemin Cai 2017.06.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
%close all;
%addpath(GetAbsolutePath('../Numerical'));
%addpath(GetAbsolutePath('../Core funcs'));
%%
if nargin < 1
    ThermalDiffus_p = 9.7e-8;             % polymer layer Mylar
end
if nargin < 2
    ThermalConduc_p = 0.15;
end
if nargin < 3
    L_polymer = 1.0e-4;
end
if nargin < 4                             
    % just for test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Case 1: Mylar+Al
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ThermalDiffus_b = 8.36e-5;      % Al base
    ThermalConduc_b = 204;
    HeatfluxFunc = '(t>0).*(25000) + (t<=0).*(0)';       % in order to test this
    
    ny_p = 20;
    ny_b = 500;
    %k = 1;                   % you can change the factor here
    %ymax_b = L*k*100;
    ymax_b = sqrt(ThermalDiffus_b*0.01*ny_b^2/19.6);
    L_base = sqrt(ThermalDiffus_b*0.01*ny_b^2/19.6);
    tmax = 8;
    nt = 500;
    h_c = 0;
    
    [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
        ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L_polymer, ymax_b, tmax, ny_p, ny_b, nt);
    Theata_ps = theata_p(end,:)';
    figure;
    plot(t, Theata_ps, '-.r');
    title('The Temperature used to recover the heat flux');
    xlabel('Time');
    ylabel('Temperature');
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Case 2: Mylar+Steel
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ThermalDiffus_b = 4.05e-6;      % Steel base
%     ThermalConduc_b = 16.3;
%     HeatfluxFunc = '(t>0).*(25000) + (t<=0).*(0)';       % in order to test this
%     
%     ny_p = 20;
%     k = 1;                   % you can change the factor here
%     ymax_b = L*k*100;
%     ny_b = 500;
%     tmax = 3;
%     nt = 300;
%     
%     [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
%     ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt);
%     Theata_ps = theata_p(end,:)';
%     figure;
%     plot(t, Theata_ps, '-.r');
%     Title('The Temperature used to recover the heat flux');
%     xlabel('Time');
%     ylabel('Temperature');


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Case 3: Mylar+Nylon
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ThermalDiffus_b = 1.26e-7;      % Nylon base 
%     ThermalConduc_b = 0.24;
%     
%     %HeatfluxFunc = '(t>0 & t<3).*(25000) +(t>3).*(50000) + (t<=0).*(0)';       % in order to test this
%     HeatfluxFunc = '(t>0).*(25000) + (t<=0).*(0)';
% 
%     ny_p = 20;
%     ymax_b = L*100;
%     ny_b = 500;
%     tmax = 3;
%     nt = 300;
%     %nt = 1500;                         % if L is small, choose larger 'nt'
%     %nt = 3500;
%     
%     % Generate the surface temperature data for testing    
%     [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
%     ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt);
%     Theata_ps = theata_p(end,:)';
%     
%     figure;
%     plot(t, Theata_ps, '-.r');
%     Title('The generated surface temparature data');
%     xlabel('time');
%     ylabel('theata');

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Case 4: Mylar+Macor
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ThermalDiffus_b = 7.53e-7;         % Macor, a material used in the reference 
%     ThermalConduc_b = 1.5;
%     
%     %HeatfluxFunc = '(t>0 & t<3).*(25000) +(t>3).*(50000) + (t<=0).*(0)';       % in order to test this
%     HeatfluxFunc = '(t>0).*(25000) + (t<=0).*(0)';
%     
%     ny_p = 20;
%     ymax_b = L*100;
%     ny_b = 500;
%     tmax = 3;
%     nt = 300;
%     %nt = 3000;                        % if L is small, choose larger 'nt'
%     
%     % Generate the surface temperature data for testing    
%     [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
%     ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt);
%     Theata_ps = theata_p(end,:)';
%     
%     figure;
%     plot(t, Theata_ps, '-.r');
%     Title('The generated surface temparature data');
%     xlabel('time');
%     ylabel('theata');
end

if max(size(t)) ~= max(size(Theata_ps))
    error('the size of Theata_ps vector and t vector should match!');
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
% if(epsilon_ba > 0.9)
%     %torrence = 1.0e-332;
%     torrence = 1.0e-32;
% else
%     torrence = 1.0e-32;
% end

%L_base = 5.0e-2;
%L_base = 5.0e+10;
%h_c = 250;                         % convective heat transfer coefficient
%h_c = 0;
h_c_head = h_c*sqrt(ThermalDiffus_b)/ThermalConduc_b;
%Delta = pi/4;
%Delta = pi/8;
%Delta = pi/16;
%Delta = pi/64;
Delta = pi/480;

ComputationalTime = cputime;
for k = 2:max(size(t))
    t_value = (k-1)*dt;
    
%     m_Func = ['(cthe.*i-' num2str(h_c_head) '.*sqrt(' num2str(t_value) ')).*exp(-i.*sqrt(' num2str(ThermalDiffus_p) '/' num2str(ThermalDiffus_b) ')*' num2str(L_base) '/' num2str(L_polymer) ')./(cthe.*i+' num2str(h_c_head) '.*sqrt(' num2str(t_value) '))'];
%     %m_Func = exp(-i*sqrt(ThermalDiffus_p/ThermalDiffus_b)*(L_base/L_polymer))
%     %m_Func = num2str(m_Func)
%     m_head_Func = ['((' num2str(epsilon_ba) '-' m_Func ')./(1-' num2str(epsilon_ba) '.*' m_Func '))'];
%     A_Func = ['(abs(' m_head_Func '))'];
%     Beta_Func = ['2.*' num2str(L_polymer) '.*cthe./sqrt(' num2str(ThermalDiffus_p) '.*' num2str(t_value) ')-angle(' m_head_Func ')'];
%     IntegralFuncStr = ['(1-' A_Func '.*' A_Func ').*exp(-cthe.^2)./(1+' A_Func '.*' A_Func '-2.*' A_Func '.*cos(' Beta_Func '))'];
%     
%     W_Fun(k-1) = 2/sqrt(pi)*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
    
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
    W_Fun(k-1) = 2/sqrt(pi)*real(IntegralComputation_Infinite_CloseZero(IntegralFuncStr, 'cthe', torrence));
    
%     m_Func = ['(cthe.*i-' num2str(h_c_head) '.*sqrt(' num2str(t_value) ')).*exp(-i.*sqrt(' num2str(ThermalDiffus_p) '/' num2str(ThermalDiffus_b) ')*' num2str(L_base) '/' num2str(L_polymer) ')./(cthe.*i+' num2str(h_c_head) '.*sqrt(' num2str(t_value) '))'];
%     %m_Func = exp(-i*sqrt(ThermalDiffus_p/ThermalDiffus_b)*(L_base/L_polymer))
%     %m_Func = num2str(m_Func)
%     m_head_Func = ['((' num2str(epsilon_ba) '-' m_Func ')./(1-' num2str(epsilon_ba) '.*' m_Func '))'];
%     A_Func = ['(abs(' m_head_Func '))'];
%     Beta_Func = ['2.*' num2str(L_polymer) '.*cthe./sqrt(' num2str(ThermalDiffus_p) '.*' num2str(t_value) ')-angle(' m_head_Func ')'];
%     IntegralFuncStr = ['(1-' A_Func '.*' A_Func ').*exp(-cthe.^2)./(1+' A_Func '.*' A_Func '-2.*' A_Func '.*cos(' Beta_Func '))'];
%     
%     W_Fun(k-1) = 2/sqrt(pi)*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
    

% ¡¡¡¡W_Fun(k-1) = 2/sqrt(pi)*IntegralComputation_new(t_value, h_c_head, ThermalDiffus_p, ThermalDiffus_b, L_base, L_polymer, epsilon_ba);
%     W_Fun(k-1) = 2/sqrt(pi)*IntegralComputation_actual(t_value, h_c_head, ThermalDiffus_p, ThermalDiffus_b, L_base, L_polymer, epsilon_ba);
    
end
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
ComputationalTime = cputime - ComputationalTime;
