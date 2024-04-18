
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
    dq=single(X)-dq_ref(i);
    
%     mask_size=6;
%     std=0.61*mask_size;
%     h=fspecial('gaussian',mask_size,std);
%     dT=imfilter(dT,h);
    
    
    dq_max(i)=0.98*max(max(abs(dq(50:150,50:300))));

    dq_matrix(:,:,i)= dT;
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


% show typical temperature image
figure(3);
a=[-8000 1000];
imagesc(dq_matrix(:,:,10),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('Heat Flux (W/m^2), t = 2 s');


% select a typical point
xc=247;
yc=93;
q_ps=squeeze(dq_matrix(yc,xc,:));

a = 1;
b = [1/4 1/4 1/4 1/4];
q_ps_1 = filter(b,a,q_ps);

% show the time history of temperature at that point
figure(4); 
plot(t1,q_ps,'+-');
xlabel('Time (s)');
ylabel('Heat Flux (W/m^2)');
%legend('x/L_c = 0.59, Run 1');
grid;
%axis([0 3 -12 10]);



% show the time history of temperature at that point
figure(5);
x1=[1:length(squeeze(dq_matrix_norm(yc,:,1)))];
y1=[1:length(squeeze(dq_matrix_norm(:,xc,1)))];

sx=0.0765; % mm/pixel
D=(1/4)*25.4; % nozzle diameter of 1/4''
x0=378;
y0=93;

x2=-(x1-x0)*sx/D;
y2=(y1-y0)*sx/D;

q_cut = squeeze(dq_matrix_norm(yc,:,5));
plot(x2,q_cut,'-');
%set(gca,'XDir','reverse');
axis([0 4.5 -1.2 0.2]);
hold on;

for i=4:2:maxframe
    q_cut = squeeze(dq_matrix_norm(yc,:,i));
    plot(x2,q_cut,'-');
    %set(gca,'XDir','reverse');
end
grid;
xlabel('x/D');
ylabel('Normalized Heat Flux');
hold off;
%axis([0 370 -1.3 0.5]);

q_add=zeros(size(dq));
for i=15:25
    q_add=q_add+dq_matrix_norm(:,:,i);
end


q_norm_avg=q_add/(25-15+1);


% show typical temperature image
figure(6);
a=[-1 0.5];
h=imagesc(x2,y2,q_norm_avg,a);
colormap('pink');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Normalized Heat Flux');
set(gca,'XDir','reverse');
axis([0 2.25 -0.44 0.44]);
hold on;

figure(61);
a=[-1.2 0.3];
h=imagesc(x2,y2,q_norm_avg,a);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Normalized Heat Flux');
%set(gca,'XDir','reverse');
axis([0 4.5 -0.88 0.88]);






x_ref_2=-(x_ref-x0)*sx/D;
y_ref_2=(y_ref-y0)*sx/D;

figure(6);
plot(x_ref_2,y_ref_2,'ok');
hold off;



figure(7)
plot(t1,dq_max/1000,'-ok',t1,dq_ref/1000,'-sr');
grid;
xlabel('Time (s)');
ylabel('Heat Flux (kW/m^2)');
axis([0 5 -10 10]);
legend('Max Value of |q_s|','q_s at Ref Location');


scale=imread('Square Grid_0.13_000001.tif');
figure(8);
imagesc(scale);
colormap(gray);
axis image;


Nu_add=zeros(size(X_Nu));
for i=15:25
    Nu_add=Nu_add+Nu_matrix(:,:,i);
end


Nu_avg=Nu_add/(25-15+1);

Nu_max=max(max(abs(Nu_avg)));

Nu_avg_norm=Nu_avg./Nu_max;

% show typical temperature image
figure(9);
a=[0.85 1.01];
h=imagesc(x2,y2,Nu_avg_norm,a);
%h=imagesc(x2,y2,Nu_avg_norm);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Normalized Nusselt Number');
%set(gca,'XDir','reverse');
axis([0 4.5 -0.88 0.88]);



% dlmwrite('Nu_avg_case3.dat', Nu_avg,'delimiter', '\t');
% % 
% 
% 


% 








