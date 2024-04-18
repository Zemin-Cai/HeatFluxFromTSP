function [HeatfluxMatrixes] = SimulatingGaussianHFField(tmax, nt, ImageWidth, ImageHeight)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating Gaussian Heat Flux Field
% 
% q_s(x, y, t) = A(t).*exp(-((x-x0).^2./(2sigma_x^2) +
% (y-y0).^2./(2sigma_y^2)));
%
% Input arguments:
%  tmax: the maximun time you want to concern
%  nt: the number of time points
%  ImageWidth: the heat flux image(field)'s width
%  ImageHeight: the heat flux image's height
% 
% Output argument:
%  HeatfluxMatrixes: the simulated heat flux field. Its size is
%                    Image_Height*Image_Width*nt
% Zemin Cai, 2009.1.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters setting
if nargin < 1
    tmax = 5;
end
if nargin < 2
    nt = 501;
end
if nargin < 3
    ImageWidth  = 20;
end
if nargin < 4
    ImageHeight = 16;
end
    
x_0 = floor(ImageWidth/2);
y_0 = floor(ImageHeight/2);

x = linspace(1, ImageWidth, ImageWidth);
y = linspace(1, ImageHeight, ImageHeight);
%t = [0:0.1:5];
t = linspace(0, tmax, nt);
%sigma_x = 1/sqrt(2)*(max(x)-x_0);
%sigma_y = 1/sqrt(2)*(max(y)-y_0);
sigma_x = 1/2*(max(x)-x_0);
sigma_y = 1/2*(max(y)-y_0);

%Amp = 5000;                     % the amplitude. The maximum heat flux 
%Amp_func = '((t<0).*(0) + (t>=0 & t<=0.1).*(2000) + (t>0.1).*(5000))';
Amp_func = '((t<0).*(0) + (t>=0 & t<=0.1).*(2000) + (t>0.1).*(50000))';
Amp = FuncComputation(Amp_func, 't', t);

%%
for i = 1:max(size(x))
    for j = 1:max(size(y))
        HFAtAllPoints(j, i, :) = Amp.*exp(-((x(i)-x_0).^2./(2*sigma_x^2)+(y(j)-y_0).^2./(2*sigma_y^2)));
    end
end
HeatfluxMatrixes = HFAtAllPoints;


%%
% show
% 1. plot one point's heat flux
close all;
figure; 
for i = 1:max(size(t))
    HFAtOnePoint(i) = HFAtAllPoints(y_0, x_0, i);
    %HFAtOnePoint(i) = HFAtAllPoints(1, 2, i);
end
plot(t, HFAtOnePoint', '-b');
xlabel('Time (Sec)');
ylabel('Heat Flux (W/m^2)');
axis([0 1 0 52000]);

% 2. show average heat flux image
sum = 0;
for i = 2:max(size(t))
    sum = sum + HFAtAllPoints(:, :, i);
end
AverageHFImage = sum./(max(size(t))-1);

figure;
%lims = [0 5000];
lims = [0 50000];
imagesc(AverageHFImage, lims);
xlabel('x (pixel)');
ylabel('z (pixel)');
%text(25, 20, 'W/m^2');
text(22.5,17, 'W/m^2');
colorbar;
%truesize;
