
% This code calculates heat flux from TSP images of a 7-deg Al cone in Run 15 (Chris Ward).
% The image registration is applied to images since the model moves
% considerably in the strat-up of the tunnel

clear all
close all

No_images=25;

for i=1:No_images
    file_name=strcat('dT_heated3_',num2str(i),'.dat');
    %file_name=strcat('dT_room_temp',num2str(i),'.dat');
    X =load(file_name);
    
    x_ref=10;
    y_ref=93;
    
    dT_ref(i)=single(X(y_ref,x_ref));
    dT=single(X)-dT_ref(i);
    
%     mask_size=6;
%     std=round(0.61*mask_size);
%     %h=fspecial('gaussian',mask_size,std);
%     h=ones(std,std)/std^2;
%     dT=imfilter(dT,h);
    
    
    dT_max(i)=0.98*max(max(abs(dT(50:150,50:300))));

    dT_matrix(:,:,i)= dT;
    dT_matrix_norm(:,:,i)= dT/dT_max(i);
     
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
subplot(3,1,1);
a=[-10 2];
imagesc(dT_matrix(:,:,5),a);
colormap('jet');
colorbar;
axis image;
%xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 0.6 s');

subplot(3,1,2);
imagesc(dT_matrix(:,:,10),a);
colormap('jet');
colorbar;
axis image;
%xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 1.8 s');


subplot(3,1,3);
imagesc(dT_matrix(:,:,20),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('t = 3.8 s');



Ts_t=load('Ts_t_analy_finite.dat');
qs_t_1=load('qs_t_analy_finite.dat');
qs_t_2=load('qs_t_num_finite.dat');
qs_t_3=load('qs_t_analy_infinite.dat');


figure(100);
plot(Ts_t(:,1),Ts_t(:,2),'+-k');
xlabel('Time (s)');
ylabel('Surface Temperature Change (K)');
%legend('x/L_c = 0.59, Run 1');
grid;


figure(101);
plot(qs_t_1(:,1),qs_t_1(:,2),'+-b',qs_t_2(:,1),qs_t_2(:,2),'o-r',qs_t_3(:,1),qs_t_3(:,2),'--k');
xlabel('Time (s)');
ylabel('Heat Flux (W/m^2)');
legend('Analytical Method for Finite Base', 'Numerical Method for Finite Base','Analytical Method for Semi-Infinite Base');
grid;



qs_x_1=load('qs_x_analy_finite.dat');
qs_x_2=load('qs_x_num_finite.dat');
qs_x_3=load('qs_x_analy_infinite.dat');

figure(102);
plot(qs_x_1(:,1),qs_x_1(:,2),'+-k',qs_x_2(:,1),qs_x_2(:,2),'o-k',qs_x_3(:,1),qs_x_3(:,2),'--k');
xlabel('x (pixels)');
ylabel('Heat Flux (W/m^2)');
legend('Analytical Method for Finite Base', 'Numerical Method for Finite Base','Analytical Method for Semi-Infinite Base');
grid;
axis([0 80 -9000 0]);























