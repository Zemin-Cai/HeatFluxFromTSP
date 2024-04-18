function [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS_FiniteBase_TemperatureDependent( ThermalDiffus_p_normal,...
    ThermalDiffus_b_normal, ThermalConduc_p_normal, ThermalConduc_b_normal, ThermalDiffus_p_Func, ThermalDiffus_b_Func,...
    q_s, ThermalConduc_p_Func, ThermalConduc_b_Func, L_p, L_b, tmax, ny_p, ny_b, nt, h_c, T_ini)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ForwardNumerical_BTCS_FiniteBase_TemperatureDependent: BTCS method for variable parameters, which
%                 means that the ThermalConductivity and Thermal Diffusivity are
%                 sensitive with the surface temperature change.
% 
% Syntax:  ForwardNumerical_BTCS_FiniteBase_TemperatureDependent(...)
%
% Input arguments:
%   ThermalDiffus_p_normal: for original setting. the thermal diffusivity value of
%                    the polymer layer under normal temperature
%   ThermalDiffus_b_normal: thermal diffusivity value of the base layer under
%                    normal temperature;
%   ThermalConduc_p_normal: thermal conductivity value of the polymer layer under
%                    normal temperature;
%   ThermalConduc_b_normal: thermal conductivity value of the base layer under
%                    normal temperature;
%   ThermalDiffus_p_Func: the thermal diffusivity function of the polymer
%                         layer(independent variable is temperature 'T')
%   ThermalDiffus_b_Func: the thermal diffusivity function of the base layer
%   HeatfluxFunc: the heat flux function qs(t)
%   ThermalConduc_p_Func: the thermal conductivity function of the polymer layer
%   ThermalConduc_b_Func: the thermal conductivity function of the base layer
%   L: length of the space domain in the polymer layer
%   L_b: length of the space domain in the base layer (positive value)
%   tmax: maximum time for the simulation
%   ny_p: number of mesh points in y direction of polymer layer
%   ny_b: number of mesh points in y direction of base layer
%   nt: number of the step in time direction
%   T_ini: initial temperature for both the polymer layer and base (in K)
%
% Output arguments:
%   theata_p: matrix of solutions of the polymer layer. theata_p(:,k) is
%             theata_p(y) at t = t(k)
%   y_p: location of finite difference nodes of the polymer layer
%   err_p: L2 norm of error evaluated at the spatial nodes on last time
%          step of the polymer layer
%   theata_b: matrix of solutions of the base layer. theata_b(:,k) is
%             theata_b(y) at t = t(k)
%   y_b: location of finite difference nodes of the base layer
%   err_b: L2 norm of error evaluated at the spatial nodes on last time
%          step of the base layer
%   t: values of time at which solution is obtained
% 
% Units of the physical quantities:
%   Thermal Conductivity k: W/m.*C   (??????)
%   Specific Heat C: W.s/kg.*C  (????????)
%   The Density: kg/m3   (??????)
%   Thermal Diffusivity a: m2/s  (?????)
%
% Reference:
%   T. Liu et.al., "Remote Surface Temperature and Heat Transfer Mapping for
%   a Waverider Model at Mach 10 Using Fluorescent Paint", 18th AIAA 
%   Aerospace Ground Testing Conference, June 20-23, 1994.
%
% Modified by Zemin Cai 2017.08.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
%close all;
%addpath(GetAbsolutePath('../Core funcs'));
%% parameters initialize
if nargin < 1
    ThermalDiffus_p_normal = 9.7e-8;
end
if nargin < 2
    ThermalDiffus_b_normal = 8.36e-5;
end
if nargin < 3
    ThermalConduc_p_normal = 0.15;
end
if nargin < 4
    ThermalConduc_b_normal = 204;
end
if nargin < 5
    ThermalDiffus_p_Func = '(0.2e-8*T+0.1e-8)';    % assume that it is a insulate layer whose 
                                         % Thermal Conductivity is 0.15, specific
                                         % heat is 1090 while the density is 1420
end
if nargin < 6
    ThermalDiffus_b_Func = '(0.2e-5*T+0.1e-5)';   % assume that it is aluminum with thermal 
                                        % conductivity 204, specific heat 904 and 
                                        % density 2700     
end
if nargin < 7
    HeatfluxFunc = '((t>=0).*(3000) + (t<0).*(0))';
end
if nargin < 8
    ThermalConduc_p_Func = '(0.1*T+0.05)';
end
if nargin < 9
    ThermalConduc_b_Func = '(10.5*T+5)';
end
if nargin < 10
    L = 50e-6;
end
if nargin < 11
    L_b = 0.002;
end
if nargin < 12
    tmax = 5;
    %tmax = 0.5;
end
if nargin < 13
    ny_p = 20;
end
if nargin < 14
    ny_b = 500;
end
if nargin < 15
    nt = 501;
    %nt = 51;
end

% mesh dissection 
dy_p = L_p/(ny_p - 1);
dy_b = L_b/(ny_b-1);
dt = tmax/(nt - 1);

% initialization
y_p = linspace(0, L_p, ny_p)';
theata_p = zeros(ny_p, nt);
err_p = 0;
%r_p = ThermalDiffus_p*dt/dy_p.^2;
%r_p = vpa(ThermalDiffus_p*dt/dy_p.^2);
%if(r_p > 1/2)
    %error('the algorithm will be unstable, return!');
%    disp('Warning! the algorithm may be unstable!');
%end
%r2_p = 1 + 2*r_p;

y_b = linspace(0, L_b, ny_b)';
theata_b = zeros(ny_b, nt);
err_b = 0;
%r_b = ThermalDiffus_b*dt/dy_b.^2;
%r_b = vpa(ThermalDiffus_b*dt/dy_b.^2);
%if(r_b > 1/2)
    %error('the algorithm will be unstable, return!');
%    disp('Warning! the algorithm may be unstable!');
%end
%r2_b = 1 + 2*r_b;

t = linspace(0, tmax, nt)';
%HeatFlux = FuncComputation(HeatfluxFunc, 't', t);
HeatFlux = q_s;

theata_p(:, 1) = 0;
theata_b(:, 1) = 0;
%theata_b(ny_b, :) = 0;
%Omiga = (ThermalConduc_b*dy_p)/(ThermalConduc_p*dy_b);

%% computation
% T_p = zeros(ny_p, 1);
% T_b = zeros(ny_b, 1);

T_p = T_ini*ones(ny_p, 1);
T_b = T_ini*ones(ny_b, 1);



for k = 2:nt
    %disp(sprintf('The %dth time',k));
    ThermalDiffus_p = FuncComputation(ThermalDiffus_p_Func, 'T', T_p);
    ThermalConduc_p = FuncComputation(ThermalConduc_p_Func, 'T', T_p);
    ThermalDiffus_b = FuncComputation(ThermalDiffus_b_Func, 'T', T_b);
    ThermalConduc_b = FuncComputation(ThermalConduc_b_Func, 'T', T_b);
    
%     if(T_p==0)
%         ThermalDiffus_p = ones(ny_p, 1).*ThermalDiffus_p_normal;
%         ThermalConduc_p = ones(ny_p, 1).*ThermalConduc_p_normal;
%     else
%         ThermalDiffus_p = FuncComputation(ThermalDiffus_p_Func, 'T', T_p);
%         ThermalConduc_p = FuncComputation(ThermalConduc_p_Func, 'T', T_p);
%     end
%     if(T_b==0)
%         ThermalDiffus_b = ones(ny_b, 1).*ThermalDiffus_b_normal;
%         ThermalConduc_b = ones(ny_b, 1).*ThermalConduc_b_normal;
%     else
%         ThermalDiffus_b = FuncComputation(ThermalDiffus_b_Func, 'T', T_b);
%         ThermalConduc_b = FuncComputation(ThermalConduc_b_Func, 'T', T_b);
%     end
    r_p = ThermalDiffus_p.*dt./(dy_p.^2);  
    r2_p = 1 + 2.*r_p;
    r_b = ThermalDiffus_b.*dt./(dy_b.^2);
    r2_b = 1 + 2.*r_b;
    Omiga = (ThermalConduc_b(1)*dy_p)/(ThermalConduc_p(1)*dy_b);
    
    % Coefficient matrix construction
    Coef_matrix = zeros((ny_p+ny_b), (ny_p+ny_b));
    Coef_matrix(1, 1) = 1 + Omiga;         % first row
    Coef_matrix(1, 2) = -1;
    Coef_matrix(1, ny_p+2) = -Omiga;
    for i = 2:(ny_p-1)
        Coef_matrix(i, i-1) = -r_p(i);
        Coef_matrix(i, i) = r2_p(i);
        Coef_matrix(i, i+1) = -r_p(i);
    end
    Coef_matrix(ny_p, ny_p-1) = -1;
    Coef_matrix(ny_p, ny_p) = 1;
    Coef_matrix(ny_p+1, 2) = -1;
    Coef_matrix(ny_p+1, ny_p+1) = 1 + Omiga;
    Coef_matrix(ny_p+1, ny_p+2) = -Omiga;
    for i = (ny_p+2):(ny_p+ny_b-1)
        Coef_matrix(i, i-1) = -r_b(i-ny_p);
        Coef_matrix(i, i) = r2_b(i-ny_p);
        Coef_matrix(i, i+1) = -r_b(i-ny_p);
    end
    Coef_matrix(ny_p+ny_b, ny_p+ny_b-1) = -ThermalConduc_b(1)/dy_b;
    Coef_matrix(ny_p+ny_b, ny_p+ny_b) = h_c + ThermalConduc_b(1)/dy_b;
    
    % update the right hand side
    b = zeros(ny_p+ny_b, 1);
    for i = 2:(ny_p-1)
        b(i) = theata_p(i,k-1);
    end
    b(ny_p) = dy_p*HeatFlux(k)/ThermalConduc_p(ny_p);
    for i = (ny_p+2):(ny_p+ny_b-1)
        b(i) = theata_b(i-ny_p,k-1);
    end
    theata = linsolve(Coef_matrix, b);
    theata_p(:, k) = theata(1:ny_p);
    theata_b(:, k) = theata((ny_p+1):end);
    T_p = theata_p(:, k);
    T_b = theata_b(:, k);
    
    T_p = theata_p(:, k)+T_ini;
    T_b = theata_b(:, k)+T_ini;
    
    
end

% T_p = 0;                                % The normal temperature;
% T_b = 0;
% for k = 2:nt
%     disp(sprintf('The %dth time',k));
%     if(T_p==0)
%         ThermalDiffus_p = ThermalDiffus_p_normal;
%         ThermalConduc_p = ThermalConduc_p_normal;
%     else
%         ThermalDiffus_p = FuncComputation(ThermalDiffus_p_Func, 'T', T_p);
%         ThermalConduc_p = FuncComputation(ThermalConduc_p_Func, 'T', T_p);
%     end
%     if(T_b==0)
%         ThermalDiffus_b = ThermalDiffus_b_normal;
%         ThermalConduc_b = ThermalConduc_b_normal;
%     else
%         ThermalDiffus_b = FuncComputation(ThermalDiffus_b_Func, 'T', T_b);
%         ThermalConduc_b = FuncComputation(ThermalConduc_b_Func, 'T', T_b);
%     end
%     r_p = vpa(ThermalDiffus_p*dt/dy_p.^2);
%     r2_p = 1 + 2*r_p;
%     r_b = vpa(ThermalDiffus_b*dt/dy_b.^2);
%     r2_b = 1 + 2*r_b;
%     Omiga = vpa((ThermalConduc_b*dy_p)/(ThermalConduc_p*dy_b));
%     
%     % Coefficient matrix construction
%     Coef_matrix = zeros((ny_p+ny_b), (ny_p+ny_b));
%     Coef_matrix(1, 1) = 1 + Omiga;         % first row
%     Coef_matrix(1, 2) = -1;
%     Coef_matrix(1, ny_p+2) = -Omiga;
%     for i = 2:(ny_p-1)
%         Coef_matrix(i, i-1) = -r_p;
%         Coef_matrix(i, i) = r2_p;
%         Coef_matrix(i, i+1) = -r_p;
%     end
%     Coef_matrix(ny_p, ny_p-1) = -1;
%     Coef_matrix(ny_p, ny_p) = 1;
%     Coef_matrix(ny_p+1, 2) = -1;
%     Coef_matrix(ny_p+1, ny_p+1) = 1 + Omiga;
%     Coef_matrix(ny_p+1, ny_p+2) = -Omiga;
%     for i = (ny_p+2):(ny_p+ny_b-1)
%         Coef_matrix(i, i-1) = -r_b;
%         Coef_matrix(i, i) = r2_b;
%         Coef_matrix(i, i+1) = -r_b;
%     end
%     Coef_matrix(ny_p+ny_b, ny_p+ny_b) = 1;
%     
%     % update the right hand side
%     b = zeros(ny_p+ny_b, 1);
%     for i = 2:(ny_p-1)
%         b(i) = theata_p(i,k-1);
%     end
%     b(ny_p) = vpa(dy_p*HeatFlux(k)/ThermalConduc_p);
%     for i = (ny_p+2):(ny_p+ny_b-1)
%         b(i) = theata_b(i-ny_p,k-1);
%     end
%     theata = linsolve(Coef_matrix, b);
%     theata_p(:, k) = theata(1:ny_p);
%     theata_b(:, k) = theata((ny_p+1):end);
%     T_p = theata_p(ny_p, k);
%     T_b = theata_b(1, k);
% end

% % just for test
% theata_ps = theata_p(ny_p, :);
% theata_numeri_p_lt = theata_p(:,nt);             % polymer layer's temperature distribution at last time
% theata_numeri_b_lt = theata_b(:,nt); 
% fprintf('\nthe numerical result of the polymer layer at last time: t = %8.6f is:\n',t(nt));
% fprintf('%8.6f %8.30f\n',[VecInverseShow(y_p) VecInverseShow(theata_numeri_p_lt)]');
% fprintf('\nthe numerical result of the base layer at last time: t = %8.6f\n',t(nt));
% fprintf('%8.6f %8.30f\n',[-y_b theata_numeri_b_lt]');
% fprintf('\nthe numerical result of the whole object at last time: t = %8.6f\n',t(nt));
% fprintf('%8.6f %8.30f\n',[[VecInverseShow(y_p); -y_b] [VecInverseShow(theata_numeri_p_lt); theata_numeri_b_lt]]');
% fprintf('\nthe numerical result on the surface at different time: L = %8.6f\n', L);
% fprintf('%8.6f %8.30f\n',[t theata_ps']');
% 
% figure;
% plot(t, theata_ps, '-');
% xlabel('Time (Sec)');
% ylabel('Temperature Change (K)');
% hold on;
% plot(t, theata_b(1, :), '--r');
% legend('Surface Temperature','Interface Temperature');