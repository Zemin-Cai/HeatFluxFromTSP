
% This code calculates heat flux from TSP images of a 7-deg Al cone in Run 15 (Chris Ward).
% The image registration is applied to images since the model moves
% considerably in the strat-up of the tunnel

clear all
close all


ThermalDiffus_p = 9.7*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.15;                 % W/m-K

ThermalDiffus_b = 6.903*10^(-5);            % Al 6061 alloy base
ThermalConduc_b = 167;

 L_p = 120e-6; % m
 L_b=9.525e-3; %m 
 h_c=0;




No_images=25;

for i=1:No_images
    
    file_name=strcat('dT_heated3_',num2str(i),'.dat');
    X =load(file_name);
    
    dT=single(X);

    dT_matrix(:,:,i)= dT;
    
    
    qs_fourier_field=ThermalConduc_p*(dT_matrix(:,:,i)-dT_matrix(10,10,1))/L_p;

    T_off=43.4+273.15; % case 4
    T_aw=317.046;
    D=6.35*10^(-3); % m
    k_air=0.0257; % W/m-K
    Ts_field=dT_matrix(:,:,i)+T_off;
    h=qs_fourier_field./(Ts_field-T_aw);
    Nu_field=h*D/k_air;
    
    X=qs_fourier_field;
    X_Nu=Nu_field;
    
    
%     file_name=strcat('qs_heated4_',num2str(i),'.dat');
%     %file_name=strcat('dT_room_temp',num2str(i),'.dat');
%     X =load(file_name);
%     
%     file_name_Nu=strcat('Nu_heated4_',num2str(i),'.dat');
%     %file_name=strcat('dT_room_temp',num2str(i),'.dat');
%     X_Nu =load(file_name_Nu);
    
    Nu_matrix(:,:,i)=single(X_Nu);
    
    x_ref=10;
    y_ref=93;
    
    dq_ref(i)=single(X(y_ref,x_ref));
    %dq=single(X)-dq_ref(i);
    dq=single(X);
    
%     mask_size=6;
%     std=0.61*mask_size;
%     h=fspecial('gaussian',mask_size,std);
%     dT=imfilter(dT,h);
    
    
    dq_max(i)=0.98*max(max(abs(dq(50:150,50:300))));

    dq_matrix(:,:,i)= dq;
    dq_matrix_norm(:,:,i)= dq/dq_max(i);
     
    figure(2);
    a=[-25000 1000];
    imagesc(dq,a);
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




figure(3);
subplot(3,1,1);
a=[-10000 1000];
imagesc(dq_matrix(:,:,5),a);
colormap('jet');
colorbar;
axis image;
%xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 0.6 s');

subplot(3,1,2);
imagesc(dq_matrix(:,:,10),a);
colormap('jet');
colorbar;
axis image;
%xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 1.8 s');


subplot(3,1,3);
imagesc(dq_matrix(:,:,20),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 3.8 s');






















