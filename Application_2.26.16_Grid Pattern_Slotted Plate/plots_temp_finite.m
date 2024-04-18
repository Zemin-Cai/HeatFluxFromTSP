
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
yc=93;
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



% show the time history of temperature at that point
figure(5);
x1=[1:length(squeeze(dT_matrix_norm(yc,:,1)))];
y1=[1:length(squeeze(dT_matrix_norm(:,xc,1)))];

sx=0.0765; % mm/pixel
D=(1/4)*25.4; % nozzle diameter of 1/4''
x0=378;
y0=93;

x2=-(x1-x0)*sx/D;
y2=(y1-y0)*sx/D;

Theata_cut = squeeze(dT_matrix_norm(yc,:,5));
plot(x2,Theata_cut,'-');
%set(gca,'XDir','reverse');
axis([0 4.5 -1.2 0.2]);
hold on;

for i=4:2:maxframe
    Theata_cut = squeeze(dT_matrix_norm(yc,:,i));
    plot(x2,Theata_cut,'-');
    %set(gca,'XDir','reverse');
end
grid;
xlabel('x/D');
ylabel('Normalized Temperature Change');
hold off;
%axis([0 370 -1.3 0.5]);



temp_add=zeros(size(dT));
for i=15:25
    temp_add=temp_add+dT_matrix_norm(:,:,i);
end


temp_norm_avg=temp_add/(25-15+1);


% show typical temperature image
figure(6);
a=[-1 0.5];
h=imagesc(x2,y2,temp_norm_avg,a);
colormap('pink');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Normalized Temperature Change');
set(gca,'XDir','reverse');
axis([0 4 -0.44 0.44]);
hold on;


% show typical temperature image
figure(60);
a=[-1.2 0.5];
h=imagesc(x2,y2,temp_norm_avg,a);
colormap('jet');
colorbar;
axis image;
xlabel('x/D');
ylabel('y/D');
title('Normalized Temperature Change');
%set(gca,'XDir','reverse');
axis([0 4.5 -0.88 0.88]);





x_ref_2=-(x_ref-x0)*sx/D;
y_ref_2=(y_ref-y0)*sx/D;

figure(6);
plot(x_ref_2,y_ref_2,'ok');
hold off;



figure(7)
plot(t1,dT_max,'-ok',t1,dT_ref,'-sr');
grid;
xlabel('Time (s)');
ylabel('Temperature Change (K)');
axis([0 5 -8 8]);
legend('Max Value of |\DeltaT|','T at Ref Location');


scale=imread('Square Grid_0.13_000001.tif');
figure(8);
imagesc(scale);
colormap(gray);
axis image;



%data0=[x' q_s_x0'];
%save q_TSP_Ward_run15.dat data0 -ascii

%data2=[x' q_s_x_fourier'];
%save('q_fourier_run1.dat', 'data2','-ascii');
%dlmwrite('q_fourier_run1.dat', data2,'delimiter', '\t');











