
%% This program reads TSP images and convert them into temperature images.    

clear all
close all

%% read raw background image for correction of ambient light effect
I_back_0=imread('background_000025.tif');

%% read sequence of background images and calculate the averaged one
% I_back_add=zeros(size(I_back_0));
% n1=1;
% n2=1;
% for i=n1:1:n2
%     if i<10
%         imagefile=strcat('background_00000',num2str(i),'.tif');
%     elseif (i>=10) && (i<100)
%         imagefile=strcat('background_0000',num2str(i),'.tif');
%     elseif (i>=100) && (i<1000)
%         imagefile=strcat('background_000',num2str(i),'.tif');
%     else
%         imagefile=strcat('background_00',num2str(i),'.tif');   
%     end
%     
%     I_back =imread(imagefile);
%     
%     I_back_add=I_back_add+double(I_back);
% end
% 
% I_back_avg=double(I_back_add)/(n2-n1+1);

%% In this vase, a simp;ification is 
I_back_avg=I_back_0;


%% read sequence of wind-off TSP images andalculate the averaged one c
% I_off_0=imread('wind-off_ref_000001.tif');
% 
% I_off_add=zeros(size(I_off_0));
% n1=1;
% n2=10;
% for i=n1:1:n2
%     if i<10
%         imagefile=strcat('wind-off_ref_00000',num2str(i),'.tif');
%     elseif (i>=10) && (i<100)
%         imagefile=strcat('wind-off_ref_0000',num2str(i),'.tif');
%     elseif (i>=100) && (i<1000)
%         imagefile=strcat('wind-off_ref_000',num2str(i),'.tif');
%     else
%         imagefile=strcat('wind-off_ref_00',num2str(i),'.tif');   
%     end
%     
%     I_off =imread(imagefile);
%     
%     I_off_add=I_off_add+double(I_off);
% end
% 
% I_off_avg=double(I_off_add)/(n2-n1+1);


%% read raw wind-on TSP images and calculate the averaged one to show the
%% overall pattern and the size of the measurement field
I_on_0=imread('heated run_wind-on_3_000001.tif');
I_on_add=zeros(size(I_on_0));

n1=1;
n2=25;
for i=n1:1:n2
    if i<10
        imagefile=strcat('heated run_wind-on_3_00000',num2str(i),'.tif');
    elseif (i>=10) && (i<100)
        imagefile=strcat('heated run_wind-on_3_0000',num2str(i),'.tif');
    elseif (i>=100) && (i<1000)
        imagefile=strcat('heated run_wind-on_3_000',num2str(i),'.tif');
    else
        imagefile=strcat('heated run_wind-on_3_00',num2str(i),'.tif');   
    end
    
    I_on =imread(imagefile);
    
    I_on_add=I_on_add+double(I_on);
end

I_on_avg=double(I_on_add)/(n2-n1+1);


%% In this case of the impinging jet, the initial TSP image of the heated plate before jet impingement
%% as the wind-off image
I_on_0=imread('heated run_wind-on_3_000001.tif');
I_on_0=double(I_on_0);

I_off_avg=I_on_0;


%% show the background image, averaged wind-off image, and intensity ratio image
figure(1);
imagesc(I_back_avg);
colormap(gray);
title('Raw Background Image');


figure(2);
imagesc(I_off_avg);
colormap(gray);
title('Raw Wind-off Image');


a=[0.98 1.1];
figure(3);
imagesc(I_on_avg./I_off_avg,a);
colormap('jet')
colorbar;
axis image;
title('Intensity Ratio Image (I_{on}/I_{off})');



%% read raw TSP calibration data and fit the data using a polynomial
data_cal=load('data_TSP_cal.dat');
Tk=data_cal(:,1);      % temperature in K
I=data_cal(:,2);

[coef1,r1]=polyfit(Tk,I,3);
I1=polyval(coef1,Tk);

figure(4);
plot(Tk,I,'o',Tk, I1,'-');
axis([280 340 0.6*10^4 1.6*10^4]); 
xlabel('T (K)');
ylabel('I (counts)');
legend('Data','Polynomial Fit');
grid;
title('TSP_Calibration (Raw Data)');


%% normalized TSP calibration
Tref=273.5+43.4;       % wind-off temperature (K) in this case 
Iref=polyval(coef1,Tref);

x=I/Iref;
y=Tk/Tref;

[coef2,r2]=polyfit(x,y,3);
y1=polyval(coef2,x);

figure(5);
plot(x,y,'o',x,y1,'-');
xlabel('I/I_{ref}');
ylabel('T/T_{ref}');
grid;
legend('Data Nomalized at Ref. Temp. (317 K)', 'Polynomial Fit');


%% select region of interst
y1=510;
y2=1160;
x1=43;
x2=1558;


%% read sequence of wind-on TSP images and convert them to temperature images 
skip=1;
minframe=1;    % the begining frame number
maxframe=25;    %the ending frame number

%% Flow-on case
k=1;
for i=minframe:skip:maxframe
    if i<10
        imagefile=strcat('heated run_wind-on_3_00000',num2str(i),'.tif');
    elseif (i>=10) && (i<100)
        imagefile=strcat('heated run_wind-on_3_0000',num2str(i),'.tif');
    elseif (i>=100) && (i<1000)
        imagefile=strcat('heated run_wind-on_3_000',num2str(i),'.tif');
    else
        imagefile=strcat('heated run_wind-on_3_00',num2str(i),'.tif');   
    end
    
    X =imread(imagefile);
    
    
    %% correcting background noise   
    I_corr=double(X)-double(I_back_avg);
    I_off_corr=double(I_off_avg)-double(I_back_avg);
    
    
    %% image ratio, I_on/I_off
    IoverIref=I_corr(y1:y2,x1:x2)./I_off_corr(y1:y2,x1:x2);
    
    %% converting to temperature by using the calibration relation       
    xx=IoverIref;
    a1=coef2(1);
    a2=coef2(2);
    a3=coef2(3);
    a4=coef2(4);
    
    factor=a1*xx.^3+a2*xx.^2+a3*xx+a4;
    
    %% temperature change (in K) from reference temperature
    dT=factor.*Tref-Tref; 
    
    
    %% scale-down and filtering images
    scale_factor=0.5; % percentage of original image
    size_filter=4; % in pixels

    [dT2,dT2] = pre_processing_a(dT,dT,scale_factor,size_filter);
    
    %% save temperature data
%     file_output=strcat('dT_heated3_',num2str(i),'.dat');
%     dlmwrite(file_output,dT2);
    
    %% show temperature image at each moment     
    figure(20);
    a=[-15 2];
    imagesc(dT2,a);
    colormap('jet');
    colorbar;
    axis image;
    title('Temperature Change (K)');
    
    
    k=k+1;
    i
end








