
% This code calculates heat flux from TSP images of a 7-deg Al cone in Run 15 (Chris Ward).
% The image registration is applied to images since the model moves
% considerably in the strat-up of the tunnel

clear all
close all

%T_off=71.3+273.15; % case 2, initial temperature (K)
T_off=43.4+273.15; % case 3
%T_off=33.8+273.15; % case 4

No_images=25;

for i=1:No_images
    file_name=strcat('dT_heated3_',num2str(i),'.dat');
    X =load(file_name);
    
    dT=single(X);

    dT_matrix(:,:,i)= dT;
     
    figure(2);
    a=[-20 5];
    imagesc(dT,a);
    colormap('jet')
    colorbar
    axis image;

    i
end


% give time
maxframe=25;
minframe=1;
framerate=5; % f/s
maxtime=maxframe/framerate;
mintime=minframe/framerate;
t1=linspace(mintime-mintime,maxtime-mintime,maxframe-minframe+1);


% show typical temperature image
figure(1);
a=[-20 5];
imagesc(dT_matrix(:,:,25),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('Temperature Change (K), t = 1.6 s in Run 11');



ThermalDiffus_p = 9.8*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.25;                 % W/m-K

ThermalDiffus_b = 8.3579*10^(-5);            % Al base
ThermalConduc_b = 204;

% ThermalDiffus_b = 4.051*10^(-6);            % Stainless Steel base
% ThermalConduc_b = 16;

% ThermalDiffus_b = 1.17*10^(-5);            % Steel base
% ThermalConduc_b = 81;


% give the coating thickness
% L = 250e-6; % meter

% calculate heat flux at that point using the 1D analytical inverse method
% [q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%    ThermalConduc_b, L, Theata_ps, t1);

L_p = 40e-6; % m
  
% select a typical point
%xc=247;
%yc=99; 
xc=[200;247;300;320;247];
yc=[50;99;120;130;99];


qs_add=[];
qs_accum=[];
Ts_add=[];
Ts_accum=[];

for i=1:length(xc)
    Theata_ps=squeeze(dT_matrix(yc(i),xc(i),:));

    a = 1;
    b = [1/4 1/4 1/4 1/4];
    Theata_ps_2 = filter(b,a,Theata_ps); 
    
    Ts=Theata_ps_2+T_off;

    [q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
        ThermalConduc_b, L_p, Theata_ps_2, t1);
 
    qs=q_s_run;
    
    qs_add=[qs_add qs];
    qs_accum=[qs_accum; qs];
    
    Ts_add=[Ts_add Ts];
    Ts_accum=[Ts_accum; Ts];
    
       
end

Ts_cases=Ts_add;
qs_cases=qs_add;


[coef_T,r]=polyfit(Ts_cases(10:end,1),qs_cases(10:end,1),1);
T_cal=[300:2:320];
qs_cal=polyval(coef_T,T_cal);


[coef_q,r]=polyfit(qs_cases(10:end,1),Ts_cases(10:end,1),1);
T_aw=polyval(coef_q,0);


figure(3);
plot(Ts_accum,qs_accum,'o',T_cal,qs_cal,'-');
xlabel('T (K)');
ylabel('q_s (W/m^2)');
grid;


figure(4); 
plot(t1,Ts,'o-');
xlabel('Time (s)');
ylabel('Temperature (K)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 5 310 318]);


% show the time history of heat flux at that point
figure(5);
plot(t1,qs_cases(:,1),'-o');
xlabel('Time (s)')
ylabel('Heat Flux (W/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
axis([0 5 -25000 1000]);


D=6.35*10^(-3); % m
k_air=0.0257; % W/m-K

h=qs(:,1)./(Ts-T_aw);
Nu=h*D/k_air;

figure(6);
plot(t1,Nu,'-o');
xlabel('Time (s)')
ylabel('Nu');
%legend('x/L_c = 0.59, Run 1');
grid;
axis([0 5 0 1200]);





data=[Ts_accum qs_accum];

% dlmwrite('Ts_qs_case4.dat', data,'delimiter', '\t');




Ts_qs_case2=load('Ts_qs_case2.dat');
Ts_qs_case3=load('Ts_qs_case3.dat');
Ts_qs_case4=load('Ts_qs_case4.dat');

[coef_2,r]=polyfit(Ts_qs_case2(10:25,1),Ts_qs_case2(10:25,2),1);
[coef_3,r]=polyfit(Ts_qs_case3(10:25,1),Ts_qs_case3(10:25,2),1);
[coef_4,r]=polyfit(Ts_qs_case4(10:25,1),Ts_qs_case4(10:25,2),1);

T_cal=[290:2:350];
qs_cal_2=polyval(coef_2,T_cal);
qs_cal_3=polyval(coef_3,T_cal);
qs_cal_4=polyval(coef_4,T_cal);

[coef_q_2,r]=polyfit(Ts_qs_case2(10:25,2),Ts_qs_case2(10:25,1),1);
[coef_q_3,r]=polyfit(Ts_qs_case3(10:25,2),Ts_qs_case3(10:25,1),1);
[coef_q_4,r]=polyfit(Ts_qs_case4(10:25,2),Ts_qs_case4(10:25,1),1);

T_aw_2=polyval(coef_q_2,0);
T_aw_3=polyval(coef_q_3,0);
T_aw_4=polyval(coef_q_4,0);


figure(30);
plot(Ts_qs_case2(:,1),Ts_qs_case2(:,2),'o',Ts_qs_case3(:,1),Ts_qs_case3(:,2),'d',Ts_qs_case4(:,1),Ts_qs_case4(:,2),'s');
xlabel('T (K)');
ylabel('q_s (W/m^2)');
axis([290 360 -4*10^4 8000]);
grid;
legend('Case 2','Case 3','Case 4');
hold on;

figure(30);
plot(T_cal,qs_cal_2,'-',T_cal,qs_cal_3,'-',T_cal,qs_cal_4,'-');
hold off;






