% This program reads the x and y coordinates from a plot image
close all;
clear all;


jet_flow_conditions;


ReD_corr=[2000:1000:200000];
%Nu_corr=0.475*(ReD_corr/2).^0.5;
Pr=0.71;
Nu_corr=0.585*(ReD_corr).^0.5*(Pr)^(0.4);



data1=load('Nu_ReD_crafton_10deg.dat');
data2=load('Nu_ReD_crafton_20deg.dat');
data3=load('Nu_ReD_goldstein_30deg.dat');

data4=[157000 247;157000 242;157000 248];

T_aw_cal=288.6; % K
T_aw_exp=317.18; % K

% factor=T_aw_cal/T_aw_exp;
factor=1;

figure(1);
plot(data1(:,1), data1(:,2)*factor,'ok',data2(:,1), data2(:,2)*factor,'>k',data3(:,1), data3(:,2)*factor,'dk',...
    data4(:,1), data4(:,2),'sk',ReD_corr, Nu_corr,'-');
grid;
xlabel('Re_D');
ylabel('max(Nu)');
axis([0 180000 0 400]);
legend('Crafton et al., 10 deg, H/D = 3.8','Crafton et al., 20 deg, H/D = 4.5','Goldstein et al., 30 deg, H/D = 4',...
    'Present Case, 15 deg, H/D = 1.5','Theory');







