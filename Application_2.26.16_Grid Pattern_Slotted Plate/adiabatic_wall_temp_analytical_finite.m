
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
     
    figure(1);
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
figure(2);
a=[-6 0];
imagesc(dT_matrix(:,:,10),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('Temperature Change (K), t = 2 s');



ThermalDiffus_p = 9.7*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.15;                 % W/m-K

% ThermalDiffus_b = 8.3579*10^(-5);            % Al base
% ThermalConduc_b = 204;

ThermalDiffus_b = 6.903*10^(-5);            % Al 6061 alloy base
ThermalConduc_b = 167;


% ThermalDiffus_b = 4.051*10^(-6);            % Stainless Steel base
% ThermalConduc_b = 16;

% ThermalDiffus_b = 1.17*10^(-5);            % Steel base
% ThermalConduc_b = 81;


% give the coating thickness
% L = 250e-6; % meter

% calculate heat flux at that point using the 1D analytical inverse method
% [q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%    ThermalConduc_b, L, Theata_ps, t1);


 L_p = 120e-6; % m
 L_b=9.525e-3; %m 
 h_c=0;
  
% select a typical point
%xc=247;
%yc=99; 
xc=[247;300;320;247;200;50];
yc=[99;120;130;99;50;100];


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

        [q_s_run, t_generalcase_fb, epsilon_ba, Delta] = HFComputation_Analytical_FiniteBase_FastConvergence(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
        ThermalConduc_b, L_p, L_b, Theata_ps_2, t1, h_c); 
    
%     [q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%         ThermalConduc_b, L_p, Theata_ps_2, t1);
 
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
plot(Ts_accum,qs_accum/1000,'o',T_cal,qs_cal/1000,'-');
xlabel('T (K)');
ylabel('q_s (kW/m^2)');
grid;
axis([308 320 -8 2]);


figure(4); 
subplot(3,1,1);
plot(t1,Ts_cases(:,1),'o-k');
%xlabel('Time (s)');
ylabel('T (K)');
%legend('x/L_c = 0.59, Run 1');
grid;
axis([0 5 306 320]);


% show the time history of heat flux at that point
subplot(3,1,2);
plot(t1,qs_cases(:,1)/1000,'-ok');
%xlabel('Time (s)')
ylabel('q_s (kW/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
axis([0 5 -8 1]);


D=6.35*10^(-3); % m
k_air=0.0257; % W/m-K

h=qs_cases(:,1)./(Ts_cases(:,1)-T_aw);
Nu=h*D/k_air;
Nu_avg=mean(Nu(10:25));

subplot(3,1,3);
plot(t1,Nu,'-ok');
xlabel('Time (s)')
ylabel('Nu');
%legend('x/L_c = 0.59, Run 1');
grid;
axis([0 5 -10 350]);


data=[Ts_accum qs_accum];

% dlmwrite('Ts_qs_case4.dat', data,'delimiter', '\t');


T_off=43.4+273.15; % case 3
%T_off=33.8+273.15; % case 4
factor_fourier=0.72;
Tb=337; % K
Tref=337;


qs_fourier=ThermalConduc_p*(Ts_cases(:,1)-Ts_cases(1,1))/L_p;

figure(40); 
plot(t1,qs_cases(:,1)/1000,'o-',t1,qs_fourier/1000,'--');
%plot(t1,qs_fourier/1000,'--');
xlabel('Time (s)');
ylabel('q_s (kW/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 5 -8 1]);



x1=[1:length(squeeze(dT_matrix(yc,:,1)))];
y1=[1:length(squeeze(dT_matrix(:,xc,1)))];

sx=0.0765; % mm/pixel
D=(1/4)*25.4; % nozzle diameter of 1/4''
x0=378;
y0=93;

x2=-(x1-x0)*sx/D;
y2=(y1-y0)*sx/D;



figure(20);
a=[-6 0];
imagesc(x2,y2,dT_matrix(:,:,10),a);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Temperature Change (K), t = 2 s');
axis([0 4.5 -0.88 0.88]);


figure(21);
a=[-8 1];
qs_fourier_field=ThermalConduc_p*(dT_matrix(:,:,10)-dT_matrix(20,20,1))/L_p;
imagesc(x2,y2,qs_fourier_field/1000,a);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Heat Flux (kW/m^2), t = 2 s');
axis([0 4.5 -0.88 0.88]);


figure(22);
a=[0 350];

D=6.35*10^(-3); % m
k_air=0.0257; % W/m-K
Ts_field=dT_matrix(:,:,10)+T_off;
h=qs_fourier_field./(Ts_field-T_aw);
Nu_field=h*D/k_air;

imagesc(x2,y2,Nu_field,a);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Nu, t = 2 s');
axis([0 4.5 -0.88 0.88]);












