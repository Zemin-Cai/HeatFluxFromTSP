
close all;
clear all;


gamma=1.4;
R=287; % J/kg-K

T0=296.2; % K
p0=172368.925; %Pa
p_atm=101325; % Pa
p_ratio=p0/p_atm;

factor=(2/(gamma-1))*(p_ratio^((gamma-1)/gamma)-1);
Me=factor^(1/2);

%Me=1; % at the nozzle 

ratio_T=1+0.5*(gamma-1)*Me^2;
Te=T0/ratio_T; % in K

ratio_p=ratio_T^(gamma/(gamma-1));
pe=p0/ratio_p;

[mu k cp Pr rho]=air_properties(Te,pe);
% mu (microPa-sec), k (mW/m-K)
vis_e=mu*10^(-6); % kg/m-K
k_e=k/1000; % W/m-K
cp_e=cp; % J/kg-K
rho_e=rho;

a=(gamma*R*Te)^0.5;
ue=a*Me;

r=Pr^0.5;
T_aw_cal=Te+r*ue^2/(2*cp_e);

D=6.35/1000; % m
Re_D=rho_e*ue*D/vis_e;










