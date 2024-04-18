%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [q_s] = function HeatFluxGenerating()
% generating a specific heat flux function
%
% Zemin Cai, 2017/8/8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[q_s] = function HeatFuxGenerating()
addpath(GetAbsolutePath('../Analytical'));
addpath(GetAbsolutePath('../Numerical'));
clc;
close all;
SaveDir = GetAbsolutePath('../../Results/SimulatedSignals');

tmax = 5;
nt = 501;
%q_s(1) = 0;
q_s = normrnd(5000, 2000, [200 1]);
q_s = (q_s-min(q_s(:)))./(max(q_s(:))-min(q_s(:)));
q_s = 3000 + q_s.*4000;

q_s(201:501) = 5000;
q_s = q_s';
%size(q_s)
t = linspace(0, tmax, nt)';

FileName = ['GerneratedRandonSignal.mat'];
FilePath = [SaveDir '\' FileName];
save(FilePath, 't', 'q_s');

figure;
plot(t, q_s, '-r');
xlabel('Time (Sec)');
ylabel('Heat Flux used(W/m^2)');
axis([0 5 2500 7500]);
%title('Heat flus used for simulation')
