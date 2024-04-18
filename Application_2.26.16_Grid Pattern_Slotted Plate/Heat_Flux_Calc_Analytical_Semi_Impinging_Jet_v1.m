
% This code calculates heat flux from TSP images of a 7-deg Al cone in Run 15 (Chris Ward).
% The image registration is applied to images since the model moves
% considerably in the strat-up of the tunnel

clear all
close all

%T_aw=346.3; % K for Case 2
T_aw=317.04; % K for Case 3
%T_aw=307.4; % K for Case 4


T_off=71.3+273.15; % case 2, initial temperature (K)
T_off=43.4+273.15; % case 3
%T_off=33.8+273.15; % case 4

factor_fourier=0.72;


No_images=25;

for i=1:No_images
    file_name=strcat('dT_heated3_',num2str(i),'.dat');
    X =load(file_name);
    
    dT=single(X);

    dT_matrix(:,:,i)= dT;
     
    figure(20);
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


% select a typical point
xc=247;
yc=99;

%xc=200;
%yc=150;

Theata_ps=squeeze(dT_matrix(yc,xc,:));

a = 1;
b = [1/4 1/4 1/4 1/4];
Theata_ps_1 = filter(b,a,Theata_ps);

% show the time history of temperature at that point
figure(2); 
plot(t1,Theata_ps_1,'o-');
xlabel('Time (s)');
ylabel('Temperature Change (K)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 3 -12 10]);


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

 L_p = 120e-6; % m
 
[q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_p, Theata_ps_1, t1);
 
qs=q_s_run;

Tb=337; % K
Tref=337;
qs_fourier=ThermalConduc_p*(Theata_ps_1-(Tb-Tref))/L_p;

Ts=Theata_ps+T_off;

D=6.35*10^(-3); % m
h=qs./(Ts-T_aw);

k_air=0.0257; % W/m-K
Nu=h*D/k_air;


figure(3); 
plot(t1,Ts,'o-');
xlabel('Time (s)');
ylabel('Temperature (K)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 3 -12 10]);


% show the time history of heat flux at that point
figure(4);
plot(t1,qs,'o-',t1,factor_fourier*qs_fourier,'--r');
xlabel('Time (s)');
ylabel('Heat Flux (W/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 3 -20000 30000]);


figure(5);
plot(t1,Nu,'o-');
xlabel('Time (s)');
ylabel('Nu');
%axis([0 5 -500 500]);
%legend('x/L_c = 0.59, Run 1');
grid;



data_out1=[t1' Ts];
data_out2=[t1' qs];
data_out3=[t1' Nu];

% dlmwrite('Ts_t_analy_infinite.dat',data_out1);
% dlmwrite('qs_t_analy_infinite.dat',data_out2);
% dlmwrite('Nu_t_analy_infinite.dat',data_out3);




% calculate heat flux distribution along the x-axis using the 1D analytical
% inverse method
j=1;
for k=1:5:370
    Theata_ps=squeeze(dT_matrix(yc,k,:));
%     [q_s, t, epsilon_ba] = HFComputation_GeneralCase( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
%     ThermalConduc_b, L, Theata_ps, t1);
    a = 1;
    b = [1/4 1/4 1/4 1/4];
    Theata_ps_2 = filter(b,a,Theata_ps);    

    [q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L_p, Theata_ps_2, t1);
 
    qs=q_s_run;
    
    qs_fourier=ThermalConduc_p*(Theata_ps_2-(Tb-Tref))/L_p;
       
%    T_off=71.3+273.15;
    Ts=Theata_ps+T_off;
    h=qs./(Ts-T_aw);
    Nu=h*D/k_air;
        
    q_s_x0(j)=mean(qs(15:16)); % avereged value over a timespan after the transient stage
    q_s_x0_fourier(j)=mean(qs_fourier(15:end));
    
    Nu_x(j)=mean(Nu(15:16));
    

    j=j+1
end




x=[1:length(q_s_x0)];
figure(6);
plot(x,q_s_x0,'-r',x, factor_fourier*q_s_x0_fourier,'--b');
%plot(x,q_s_x0_fourier,'--b');
xlabel('x (m)');

ylabel('q (W/m^2)');
legend('General Solution', 'Discrete Fourier Law');
grid;
%axis([0 0.4 0 8000]);0


x=[1:length(Nu_x)];
figure(7);
plot(x,Nu_x,'-r');
xlabel('x (m)');
ylabel('Nu)');
legend('General Solution', 'Discrete Fourier Law');
grid;
%axis([0 0.4 0 8000]);0







No_images=25;

for i=1:No_images
%     file_name=strcat('dT_heated3_',num2str(i),'.dat');
%     %file_name=strcat('dT_room_temp',num2str(i),'.dat');
%     X =load(file_name);
    
    dT=dT_matrix(:,:,i);
    qs=factor_fourier*ThermalConduc_p*(dT-(Tb-Tref))/L_p;
    
    
%     T_off=71.3+273.15;
    Ts=dT+T_off;
    
    h=qs./(Ts-T_aw);
    Nu_field=h*D/k_air; 
    
    
%     file_qs=strcat('qs_heated4_',num2str(i),'.dat');
%     dlmwrite(file_qs,qs);
% 
%     file_Nu=strcat('Nu_heated4_',num2str(i),'.dat');
%     dlmwrite(file_Nu,Nu_field);
 
    i
end


figure(200);
a=[-1*10^4 -2*10^3];
imagesc(qs,a);
colormap('pink');
colorbar
axis image;


figure(201);
a=[230 250];
imagesc(Nu_field,a);
colormap('pink');
colorbar
axis image;

figure(203);
plot(Nu_field(99,:));

figure(204);
plot(qs(99,:));


data_out4=[x' q_s_x0'];
data_out5=[x' Nu_x'];

% dlmwrite('qs_x_analy_infinite.dat',data_out4);
% dlmwrite('Nu_x_analy_infinite.dat',data_out5);














