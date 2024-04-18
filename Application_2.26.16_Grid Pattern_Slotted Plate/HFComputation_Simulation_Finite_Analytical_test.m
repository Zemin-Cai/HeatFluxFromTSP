%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HFComputation_Simulation_Finite_Analytical_test:
% recover the heat flux at the polymer surface on a base with a finite 
% thickness using our analytical method for simulation
%
% Zemin Cai 2017/6/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc;
close all;
addpath(GetAbsolutePath('Analytical'));
addpath(GetAbsolutePath('Numerical'));
addpath(GetAbsolutePath('Core funcs'));

% parameters
ThermalDiffus_p = 9.7e-8;             % Mylar insulate layer
ThermalConduc_p = 0.15;

L_polymer = 20e-6;
L_base = 10e-3;                      % with a finite-thickness base
ny_p = 20;
ny_b = 500;
tmax = 3;
nt = 301;
h_c = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% heat flux used for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HeatfluxFunc_index = 2;                   % 1 - 4
switch HeatfluxFunc_index
    case 1
        HeatfluxFunc = '(t>=0).*(25000) + (t<0).*(0)';
    case 2
        HeatfluxFunc = '(t>=0 & t<=1.5).*(25000) +(t>1.5).*(50000) + (t<0).*(0)';
    case 3
        HeatfluxFunc = '(t>0 & t<=1.5).*(50000/1.5.*t) + (t>1.5).*((50000/1.5).*(3-t)) + (t<=0).*(0)';
    case 4
        HeatfluxFunc = '(t>=0 & t<=1).*(25000*(2*t-t.^2)) + (t>1).*(25000) + (t<0).*(0)';
        tmax = 5;
        nt = 501;    
end
t = linspace(0, tmax, nt)';
Heatflux_Original = FuncComputation(HeatfluxFunc, 't', t);

% switches
Test1On_Or_Off = 1;            % turn on; Mylar + Al
Test2On_Or_Off = 0;            %          Mylar + Steel
Test3On_Or_Off = 0;            %          Mylar + Nylon
Test4On_Or_Off = 0;            %          Mylar + Macor
Test5On_Or_Off = 0;            %          Mylar + ArtificialBase
Test6On_Or_Off = 0;            %          K. Asai's data


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 1: Mylar + Al
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Test1On_Or_Off)
    disp('The Mylar+Al case is testing......');

    ThermalDiffus_b = 8.36e-5;
    ThermalConduc_b = 204;
    
    disp(sprintf('\tSimulate the surface temperature with finite base ...'));
    %ny_b = 3000;            % 1.5meters
    %ny_b = 7000;            % 3 meters
    ny_b = 500;              % 10 micro meters
    [theata_p_Finite, y_p, err_p, theata_b_Finite, y_b, err_b, t_temperature] = ForwardNumerical_BTCS_FiniteBase( ThermalDiffus_p,...
        ThermalDiffus_b, HeatfluxFunc, ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt, h_c);
    Theata_ps_Finite = theata_p_Finite(end,:)';
    Theata_Interface_Finite = theata_b_Finite(1,:)';
    
   % plot
   % Original Heat Flux
    figure;
    plot(t_temperature, Heatflux_Original, '-.');
    xlabel('Time (sec)');
    ylabel('Original Heat Flux');
    %axis([0 3 0 26000]);
    title('The original heat flux used for simulation');
     
    figure;
    plot(t_temperature, Theata_ps_Finite, '-b');
    xlabel('Time (Sec)');
    ylabel('Temperature Change (K)');
    hold on;
    plot(t_temperature, Theata_Interface_Finite, '--r');
    grid on;
    legend('Tsurface', 'Tinterface');
    title('The simulated surface temperature with finite-thick base');    
    
    disp(sprintf('\tInverse heat transfer solution for polymer layer ...'));   
    disp(sprintf('\tRecovering with finite base...'));
    generaltime = cputime; 
    [q_s_generalcase_fb, t_generalcase_fb, epsilon_ba, Delta] = HFComputation_Analytical_FiniteBase_FastConvergence(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
        ThermalConduc_b, L_polymer, L_base, Theata_ps_Finite, t_temperature, h_c); 
    generaltime = cputime - generaltime;
    disp(sprintf('\t\tTime used: %d s', generaltime));
    
    % plot
    % reconvered heat flux with finite-thick base
    figure;
    plot(t_generalcase_fb, q_s_generalcase_fb, '--b');
    xlabel('Time (Sec)');
    ylabel('Heat Flux (W/m^2)');
    hold on;
    plot(t_generalcase_fb, Heatflux_Original, '--r');
    %axis([0 5 0 8500]);
    grid on;
    legend(['Recoverd Heat Flux (finite-thick base) with \delta = PI/' num2str(pi/Delta)], 'Original Heat Flux');
 
    
    q_data=[t_generalcase_fb Heatflux_Original t_generalcase_fb q_s_generalcase_fb];
    T_data=[t_temperature Theata_ps_Finite t_temperature Theata_Interface_Finite];
    
%     dlmwrite('T_step_20mic_10mm_delta_pid480_h100.dat',T_data);
%     dlmwrite('q_step_20mic_10mm_delta_pid480_h100.dat',q_data);
% % %     
    
    
    
%     % save data
%     ResultsSaveDir = GetAbsolutePath('../Results');
%     SaveDir = [ResultsSaveDir '\Simulation\HeatFluxRecover\'];
%     if(not(exist(SaveDir)))
%         mkdir(SaveDir);
%     end
%     FileName = ['HeatFlux_Mylar_Al_' num2str(L_polymer) '_' num2str(ny_p) '_' num2str(ymax_b) '_' num2str(ny_b) '_' num2str(tmax) '_' num2str(nt) '.mat'];
%     t = t_temperature;
%     q_s_AsGeneralCase = q_s_generalcase;
% %     q_s_AsSpecialCase = q_s_highconduc;
%     SurfaceTemperature = Theata_ps;
%     save([SaveDir FileName],'epsilon_ba','t','SurfaceTemperature','q_s_AsGeneralCase');
%     %save([SaveDir FileName],'epsilon_ba','t','SurfaceTemperature','q_s_AsSpecialCase','q_s_AsGeneralCase');
%     disp('The result have been saved!');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test 2: Mylar + Steel (General case method and high-conductivity method)
% coming soon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
