
% This code calculates heat flux from TSP images of a 7-deg Al cone in Run 15 (Chris Ward).
% The image registration is applied to images since the model moves
% considerably in the strat-up of the tunnel

clear all
close all

No_images=25;

for i=1:No_images
    file_name=strcat('dT_heated4_',num2str(i),'.dat');
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
figure(3);
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
Theata_ps=squeeze(dT_matrix(yc,xc,:));

a = 1;
b = [1/4 1/4 1/4 1/4];
Theata_ps_1 = filter(b,a,Theata_ps);

% show the time history of temperature at that point
figure(4); 
plot(t1,Theata_ps,'+-');
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
L = 250e-6; % meter

% calculate heat flux at that point using the 1D analytical inverse method
[q_s_run, t, epsilon_ba] = HFComputation_GeneralCase(ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L, Theata_ps, t1);
q_s=q_s_run;

% show the time history of heat flux at that point
figure(5);
plot(t1,q_s,'-');
xlabel('Time (s)')
ylabel('Heat Flux (W/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 3 -20000 30000]);


% calculate heat flux distribution along the x-axis using the 1D analytical
% inverse method
j=1;
for k=1:5:370
    Theata_ps=squeeze(dT_matrix(yc,k,:));
    [q_s, t, epsilon_ba] = HFComputation_GeneralCase( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
    ThermalConduc_b, L, Theata_ps, t1);
    q_s_x0(j)=mean(q_s(15:end)); % avereged value over a timespan after the transient stage

    j=j+1
end


L = 250e-6; % meter
Tb=337; % K
Tref=337;

j=1;
for k=1:5:370
    Theata_ps=squeeze(dT_matrix(yc,k,:));
    a = 1;
    b = [1/4 1/4 1/4 1/4];
    Theata_ps_1 = filter(b,a,Theata_ps); 
    q_s_x_fourier(j)=ThermalConduc_p*(mean(Theata_ps_1(15:end))-(Tb-Tref))/L;
    
    j=j+1
end


x=[1:length(q_s_x0)];
figure(6);
plot(x,q_s_x0,'-r',x, 1.05*q_s_x_fourier(1:74),'--');
xlabel('x (m)');
ylabel('q (W/m^2)');
legend('General Solution', 'Discrete Fourier Law');
grid;
%axis([0 0.4 0 8000]);0



No_images=25;

for i=1:No_images
    file_name=strcat('dT_heated2_',num2str(i),'.dat');
    %file_name=strcat('dT_room_temp',num2str(i),'.dat');
    X =load(file_name);
    
    qs=1.05*ThermalConduc_p*(X-(Tb-Tref))/L;
    
    %file_output=strcat('qs_heated2_',num2str(i),'.dat');
    %dlmwrite(file_output,qs);
    
     
    figure(200);
    a=[-10^4 10];
    imagesc(qs,a);
    colormap('pink');
    colorbar
    axis image;

    i
end














%data0=[x' q_s_x0'];
%save q_TSP_Ward_run15.dat data0 -ascii

%data2=[x' q_s_x_fourier'];
%save('q_fourier_run1.dat', 'data2','-ascii');
%dlmwrite('q_fourier_run1.dat', data2,'delimiter', '\t');











