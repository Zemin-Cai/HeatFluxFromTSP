
%Code to read appropriate frames from the multi-tif files from PCO-camera, 
%apply 2-D filter, find Temperature, and trace temperature in time for a given point 

clear all
close all

% Flow-off without heating
I_back_0=imread('background_000025.tif');

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

I_back_avg=I_back_0;


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

 %Tref=273.5+33.8; % K, wind-off temperature in Case 4
 Tref=273.5+43.4; % K, wind-off temperature in Case 3
 %Tref=273.5+71.3; % K, wind-off temperature in Case 2

I_on_0=imread('heated run_wind-on_3_000001.tif');
I_on_0=double(I_on_0);

I_on_add=zeros(size(I_on_0));

n1=1;
n2=25;
for i=n1:1:n2
    if i<10
        imagefile=strcat('heated run_wind-on_3_00000',num2str(i),'.tif');
    elseif (i>=10) && (i<100)
        imagefile=strcat('heated run_wind-on_3_0000',num2str(i),'.tif');
    elseif (i>=100) && (i<1000)
        imagefile=strcat('wind-on_1_000',num2str(i),'.tif');
    else
        imagefile=strcat('wind-on_1_00',num2str(i),'.tif');   
    end
    
    I_on =imread(imagefile);
    
    I_on_add=I_on_add+double(I_on);
end

I_on_avg=double(I_on_add)/(n2-n1+1);



I_off_avg=I_on_0;


figure(1);
imagesc(I_back_avg);
colormap(gray);


figure(2);
imagesc(I_off_avg);
colormap(gray);

a=[0.98 1.2];
figure(3);
imagesc(I_on_avg./I_off_avg,a);
colormap('jet')
colorbar;
axis image;

ruler=imread('Square Grid_0.13_000001.tif');
figure(30);
imagesc(ruler);
colormap('gray')
colorbar;
axis image;



% read TSP calibration data
 data_cal=load('data_TSP_cal.dat');
 Tk=data_cal(:,1);
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
 legend('Data Nomalized at 317 K', 'Polynomial Fit');




skip=1;
minframe=1;    %Tunnel start occurs at (around) this frame number
maxframe=25;    %Tunnel un-start occurs at (around) this frame number

y1=510;
y2=1160;
x1=43;
x2=1558;


% Flow-on case
k=1;
for i=minframe:skip:maxframe
    if i<10
        imagefile=strcat('heated run_wind-on_3_00000',num2str(i),'.tif');
    elseif (i>=10) && (i<100)
        imagefile=strcat('heated run_wind-on_3_0000',num2str(i),'.tif');
    elseif (i>=100) && (i<1000)
        imagefile=strcat('wind-on_1_000',num2str(i),'.tif');
    else
        imagefile=strcat('wind-on_1_00',num2str(i),'.tif');   
    end
    
    X =imread(imagefile);
       
    I_corr=double(X)-double(I_back_avg);
    I_off_corr=double(I_off_avg)-double(I_back_avg);

    IoverIref=I_corr(y1:y2,x1:x2)./I_off_corr(y1:y2,x1:x2);
    
           
    xx=IoverIref;
    a1=coef2(1);
    a2=coef2(2);
    a3=coef2(3);
    a4=coef2(4);
    
    factor=a1*xx.^3+a2*xx.^2+a3*xx+a4;;
    dT=factor.*Tref-Tref;
    
    
    level_wavelet=2; % 0 = no downsampling, 1 = downsampling by 2, 2 = downsampling by 2^2=4
    size_filter=4; % in pixels

    [dT2,dT2] = pre_processing(dT,dT,level_wavelet,size_filter);
    
    
    %dT_matrix(1:20:(y2-y1),1:20:(x2-x1),k)= dT2(1:20:(y2-y1),1:20:(x2-x1));
    %dT_matrix_1(1:(y2-y1),1:(x2-x1),k)= dT2(1:(y2-y1),1:(x2-x1));
    
%     file_output=strcat('dT_heated2_',num2str(i),'.dat');
%     dlmwrite(file_output,dT2);
    
         
    figure(20);
    a=[-20 5];
    imagesc(dT2,a);
    
    %a=[0.98 1.2];
    %imagesc(IoverIref,a);
    colormap('jet')
    colorbar;
    axis image;
    
    
    k=k+1;
    i
end




% close all;
% clear all;
% 
% I_off =imread('wind-on_2_000001.tif');
% I_on =imread('wind-on_2_000005.tif');
% 
% I_off=double(I_off);
% I_on=double(I_on);
% 
% figure(100);
% a=[0.9 1.2];
% 
% imagesc(I_on./I_off,a);
% colorbar;
% size image;



