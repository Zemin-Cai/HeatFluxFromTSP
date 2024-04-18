function [q_s, x, z, t] = HFComputation_main_3D( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L, Theata_ps, x, z, t)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_main_3D: main routine to compute heat flux concerning
%                        effect of lateral heat conduction (3D case)
% 
% Syntax:  HFComputation_main_3D(...)
%
% Input arguments:
%   ThermalDiffus_p: the thermal diffusivity of the polymer layer
%   ThermalConduc_p: the thermal conductivity of the polymer layer
%   ThermalDiffus_b: the thermal diffusivity of the base layer
%   ThremalConduc_b: the thermal conductivity of the base layer
%   L: length of the space domain in the polymer layer
%   Theata_ps: the temparature change in the polymer surface (a 3D array: 
%              Theata_ps(x,z,t))
%   x: the x vector
%   z: the z vector
%   t: the time vector
%
% Output arguments:
%   q_s: the computed heat flux in the selected t point (a 3D array: q_s(x,z,t))
%   x: the x vector
%   z: the z vector
%   t: the time vector
%
% Zemin Cai 2008.08.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% initialization
if(nargin<1)
    ThermalDiffus_p = 9.7e-8;
end
if(nargin<2)
    ThermalConduc_p = 0.15;
end
if(nargin<3)
    ThermalDiffus_b = 8.36e-5;
end
if(nargin<4)
    ThermalConduc_b = 204;
end
if(nargin<5)
    L = 100e-6;
end
if(nargin<6)
    Theata_ps = rand(10,20,120);
end
if(nargin<7)
    x = linspace(0,1,10);
end
if(nargin<8)
    z = linspace(0,3,20);
end
if(nargin<9)
    t = linspace(0,5,120);
end

x_length = max(size(x));
z_length = max(size(z));
t_length = max(size(t));
dx = (x(end)-x(1))/(x_length-1);
dz = (z(end)-z(1))/(z_length-1);
dt = (t(end)-t(1))/(t_length-1);

q_s(:, :, 1) = 0;

%%
% computation
epsilon = (ThermalConduc_p*sqrt(ThermalDiffus_b))/(ThermalConduc_b*sqrt(ThermalDiffus_p));
epsilon_ba = (1 - epsilon)/(1 + epsilon);
torrence = 1.0e-64;
for k = 2:t_length
    t_value = (k-1)*dt;
    IntegralFuncStr = strcat(['(exp(-cthe.^2)./(1+', num2str(epsilon_ba), '^2-2*', num2str(epsilon_ba), '*cos(2*',...
           num2str(L), './sqrt(', num2str(ThermalDiffus_p*t_value), ').*cthe)))']);
    Whead_Value(k-1) = 2/sqrt(pi)*IntegralComputation_Infinite(IntegralFuncStr, 'cthe', torrence);
end
%for i = 1:x_length
for i = 20:20
    for j = 1:z_length
        disp(sprintf('Now processing the (%d,%d)th pixel point', j, i));
        for k = 2:t_length
            
            Tao_length = k;
            for Tao_index = 1:(Tao_length-1)
                % smooth the surface temperature
                % inner integration
                for z_pie_index = 1:z_length
                    g1 = 1/(4*pi*ThermalDiffus_p*t(k-Tao_index+1)).*exp(-((x(i)-x(:)).^2+(z(j)-z(z_pie_index)).^2)./(4*ThermalDiffus_p*t(k-Tao_index+1)));
                    Theata_ps_func = Theata_ps(:,z_pie_index,Tao_index);
                    FuncVec = g1.*Theata_ps_func;
                    InnerInteg(z_pie_index) = IntegralComputation_MultiRules(FuncVec, x(1), x(end), dx); 
                end
                % outer integration
                Integral_Result(Tao_index) = IntegralComputation_MultiRules(InnerInteg, z(1), z(end), dz);
            end
            Smoothed_Theata_ps1 = Integral_Result;
            
            for Tao_index = 1:(Tao_length-1)
                % smooth the surface temperature
                % inner integration
                for z_pie_index = 1:z_length
                    g2 = 1/(4*pi*ThermalDiffus_p*t(k-Tao_index+1).^2).*(1+((x(i)-x(:)).^2+(z(j)-z(z_pie_index)).^2)./(4*ThermalDiffus_p*t(k-Tao_index+1))).*exp(-((x(i)-x(:)).^2+(z(j)-z(z_pie_index)).^2)./(4*ThermalDiffus_p*t(k-Tao_index+1)));
                    Theata_ps_func = Theata_ps(:,z_pie_index,Tao_index);
                    FuncVec = g2.*Theata_ps_func;
                    InnerInteg(z_pie_index) = IntegralComputation_MultiRules(FuncVec, x(1), x(end), dx); 
                end
                % outer integration
                Integral_Result(Tao_index) = IntegralComputation_MultiRules(InnerInteg, z(1), z(end), dz);
            end
            Smoothed_Theata_ps2 = Integral_Result;
            
            if k == 2
                sum1 = Whead_Value(1)/sqrt(t(2))*Smoothed_Theata_ps1(1);
                sum2 = Whead_Value(1)/sqrt(t(2))*Smoothed_Theata_ps2(1)*dt;
            else
                sum1 = 0;
                sum2 = 0;
                for Tao_index = 2:(k-1)
                    t_value1 = t(k) - t(Tao_index);
                    if(t_value1 ~= 0)
                        Whead_Value1 = Whead_Value(k-Tao_index);
                    else
                        Whead_Value1 = 0;
                    end
                    t_value2 = t(k) - t(Tao_index-1);
                    if(t_value2 ~= 0)
                        Whead_Value2 = Whead_Value(k-Tao_index+1);
                    else
                        Whead_Value2 = 0;
                    end

                    % sum
                    sum1 = sum1 + (Smoothed_Theata_ps1(Tao_index) - Smoothed_Theata_ps1(Tao_index-1))*(Whead_Value1 + Whead_Value2)/(sqrt(t_value1) + sqrt(t_value2));
                    sum2 = sum2 + (Smoothed_Theata_ps2(Tao_index) + Smoothed_Theata_ps2(Tao_index-1))*(Whead_Value1 + Whead_Value2)*dt/(2*(sqrt(t_value1)+sqrt(t_value2)));
                end
            end
            q_s(i, j, k) = ThermalConduc_p*(1 - epsilon_ba^2)/sqrt(pi*ThermalDiffus_p)*(sum1 + sum2);
        end
    end
end
