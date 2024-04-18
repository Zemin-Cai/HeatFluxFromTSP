function [q_s, t, epsilon_ba] = HFComputation_GeneralCase( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L, Theata_ps, t)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_GeneralCase: calculate the heat flux directly for the general 
%                            case
% 
% Syntax:  HFComputation_GeneralCase(...)
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
% Zemin Cai 2008.06.25
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
    L = 1.0e-4;
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
    tmax = 3;
    nt = 300;
    
    [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS( ThermalDiffus_p,...
        ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L, ymax_b, tmax, ny_p, ny_b, nt);
    Theata_ps = theata_p(end,:)';
    figure;
    plot(t, Theata_ps, '-.r');
    Title('The Temperature used to recover the heat flux');
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
epsilon_ba = (1 - epsilon)/(1 + epsilon);

if(epsilon_ba > 0.9)
    %torrence = 1.0e-332;
    torrence = 1.0e-64;
else
    torrence = 1.0e-64;
end

ComputationalTime = cputime;
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

% figure;
% % plot(t(1:end-1), Whead_Value, '-');
% % hold on;
% % plot(t(1:end-1), 200*sqrt(t(1:end-1)), '-r');
% %t(2:end)
% %Whead_Value
% plot(t(2:end), Whead_Value'./t(2:end));

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

ComputationalTime = cputime - ComputationalTime;
%fprintf('Heat flux computational time is: %d', ComputationalTime);
%% error computation

%% plot result


%% save data
% SaveDir = 'C:\Zemin Cai\code\Forward\Finite Difference\Results\';
% FileName = ['HFGeneralCase_Mylar_Al_' num2str(L) '_' num2str(ny_p) '_' num2str(ymax_b) '_' num2str(ny_b) '_' num2str(tmax) '_' num2str(nt) '.mat'];
% %FileName = ['HFGeneralCase_Mylar_Nylon_' num2str(L) '_' num2str(ny_p) '_' num2str(ymax_b) '_' num2str(ny_b) '_' num2str(tmax) '_' num2str(nt) '.mat'];
% %FileName = ['HeatFlux_Mylar_Macor_' num2str(L) '_' num2str(ny_p) '_' num2str(ymax_b) '_' num2str(ny_b) '_' num2str(tmax) '_' num2str(nt) '.mat'];
% if not(exist(SaveDir))
%     mkdir(SaveDir);
% end
% Result = q_s;
% % Result_New = q_s_New;
% % Result_Cook_Mylar = q_s_Cook_Mylar;
% % Result_Cook_Nylon = q_s_Cook_Nylon;
% % SurfaceTemperature = Theata_ps;
% save([SaveDir FileName],'t','Result','ComputationalTime');
% disp('The result have been saved!');
