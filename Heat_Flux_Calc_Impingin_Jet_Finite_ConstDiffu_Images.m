%% This program calculates the heat flux by using the analytical inverse heat transfer solution
%% for a finite base in the imping jet. As an example, the analytical method uses the summation of a sequence of temperature images
%% directly without interpolation for fast comutation.  


clear all
close all

addpath(GetAbsolutePath('Core Funcs'));
addpath(GetAbsolutePath('Analytical'));
addpath(GetAbsolutePath('Numerical'));



%% read surface temperature files and form a 3D matrix
No_images=25;
for i=1:No_images
    file_name=strcat('dT_heated4_',num2str(i),'.dat');
    X =load(file_name);    
    dT=single(X);

    dT_matrix(:,:,i)= dT;
end

%% show typical temperature image
figure(1);
a=[-20 5];
imagesc(dT_matrix(:,:,25),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('Typical Temperature Change (K), t = 1.6 s in Run 11');



%% define time series
minframe=1;
maxframe=25;
framerate=5; % f/s
maxtime=maxframe/framerate;
mintime=minframe/framerate;
t=linspace(mintime-mintime,maxtime-mintime,maxframe-minframe+1);


%% give the thermal properties of the polymer and base material
ThermalDiffus_p = 9.7*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.15;                 % W/m-K

ThermalDiffus_b = 6.903*10^(-5);            % Al 6061 alloy base
ThermalConduc_b = 167;


%% give the thicknesses of the polymer coating and base, and heat transfer coefficient at the base bottom 
 L_p = 120e-6;      % polymer coating thickness (m)
 L_b=9.525e-3;      % base thickness (m) 
 h_c=0;             % heat transfer coefficient at the base bottom


 %% give relevant temperatues 
Tb=43.4+273.15;         %  initial base temperature (K) for TSP measurements
T_off=43.4+273.15;      % wind-off temperature (K) in case 3, which is equal to Tb
T_aw=317.1;             % adiabatic wall temperature (K) determined for Case 3

%% jet parameters
 D=6.35*10^(-3);        % jet diameter (m)
 k_air=0.0257;          % air thermal conductivity (W/m-K)

 
 %% calulate the heat flux by using the analytical method in the image plane:
 
[W_Fun, t] = W_fun_nathan( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
              ThermalConduc_b, L_p, L_b, t, h_c);
 
[row,col,nt]=size(dT_matrix);
 
[W_Fun_seq, t] = W_fun_nathan( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
                  ThermalConduc_b, L_p, L_b, t, h_c);
for i=1:row
        for j=1:col
             for k=1:nt
                   W_Fun(i,j,k)=W_Fun_seq(k);
             end
         end
 end
                          
PreTime=cputime;          
  
[qs_matrix, t, Delta] = HF_Analy_FiniteBase_W1_images( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
             ThermalConduc_b, L_p, L_b, dT_matrix, t, h_c, W_Fun);
                    
Time_Anal=cputime-PreTime;


%% save the calculated heat flux data in sequence
for i=1:nt
    q2=qs_matrix(:,:,i);
    figure(20);
    a=[-6000 0];
    imagesc(q2);
    colormap('jet');
    colorbar;
    axis image;
    
    file_output=strcat('q_cal_cone_',num2str(i),'.dat');
    %dlmwrite(file_output,q2);  
end




















