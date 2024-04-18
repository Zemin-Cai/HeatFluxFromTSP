function [TemperatureMatrixes_p, TemperatureMatrixes_b, SurfaceTemperatureMatrixes, dx, dy, dz_p, dz_b, dt] = ForwardNumerical_3D_FiniteBase( ThermalConduc_p,...
    ThermalDiffus_p, ThermalConduc_b, ThermalDiffus_b, HeatFluxMatrixes, L_polymer, L_base, Model_Width, Model_Height, tmax, ny_p, ny_b, nt, h_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Numerical Temperature Simulation for 3D case (Lateral effect) on
% finite base
%
% Syntax: ForwardNumerical_3D_FiniteBase(...)
%
% Input arguments:
%  ThermalConduc_p: the thermal conductivity of the polymer layer
%  ThermalDiffus_p: the thermal diffusivity of the polymer layer
%  ThermalConduc_b: the thermal conductivity of the base
%  ThermalDiffus_b: the thermal diffusivity of the base
%  HeatFluxMatrixes: the heat flux field effective on the model surface
%                    HeatFluxMatrixes is a Model_Width*Model_Height*nt
%                    structure
%  L_polymer: the length of the polymer layer, unit is m
%  L_base: the length of the base
%  Model_Width: the width of the model, unit is m. The step length between
%               two grid points in this direction(x-direction) is
%               Model_Width/(Image_Width-1) where Image_Width is the
%               pixels number along the image's row direction. Image_Width
%               can be obtained from the size of 'HeatFluxMatrixes'
%  Model_Height: the height of the model, unit is m. The step length
%                between two grid points in this y-direction is
%                Model_Height/(Image_Height-1) where Image_Height is the
%                pixels number along the image's column direction. 
%  tmax: the length of the time concerned
%  ny_p: the number of grid nodes in the polymer layer
%  ny_b: the number of the grid nodes in the base
%  nt: the number of time/step points
%
% Output arguments:
%  TemperatureMatrixes_p: the temperature at each point of the 3D
%                         structure of the polymer part. The dimension of 
%                         TemperatureMatrixes_p is 
%                         Image_Height*Image_Width*ny_p*nt
%  TemperatureMatrixes_b: the temperature at each point of the 3D structure
%                         of the base part. The dimension is
%                         Image_Height*Image_Width*ny_b*nt
%  SurfaceTemperatureMatrixes: the structure whitch contains the
%                              temperature data of each point
%  dx: the step length along the x_direction: Model_Width/(Image_Width-1)
%  dy: the step length along the y_direction: Model_Height/(Image_Height-1)
%  dz_p: the step length along the z_direction in the polymer layer: dz_p =
%        L/(ny_p-1)
%  dz_b: the step length along the z_direction in the base: dz_b =
%        ymax_b/(ny_b-1)
%  dt: the time step length: dt = tmax/(nt-1)
%
% Zemin Cai. 2017.11.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clc;
%close all;
%clear all;
addpath(GetAbsolutePath('../Core funcs'));
%addpath(GetAbsolutePath('../ImageProcessing'));
addpath(GetAbsolutePath('../Analytical'));

%% parameters initialize
if nargin < 1
    ThermalConduc_p = 0.16;      % the default polymer is PVC
end
if nargin < 2
    ThermalDiffus_p = 8e-8;
end
if nargin < 3
    ThermalConduc_b = 204;       % the defalut base is Al
end
if nargin < 4
    ThermalDiffus_b = 8.36e-5;
end
% if nargin < 3
%     ThermalConduc_b = 0.24;       % the defalut base is Nylon6
% end
% if nargin < 4
%     ThermalDiffus_b = 1.26e-7;
% end
if nargin < 5
    HeatFluxMatrixes = SimulatingGaussianHFField(1, 101, 20, 16);       % simulate a heat flux field of Gaussian distribution
    SimulationMethod = 'Gaussian';
%     HeatFluxMatrixes = SimulatingMutationHFField(1, 101, 20, 16);  
%     SimulationMethod = 'Mutation';
end
if nargin < 6
    L_polymer = 1.0e-5;                         % 0.01mm
end
if nargin < 7
    L_base = 0.002;
end
if nargin < 8
    Model_Width = 0.15;                % the bigger the model, the effect of lateral is more small. The size of the model can not be too small                 
end
if nargin < 9
    Model_Height = 0.12;
end
if nargin < 10
    tmax = 1;
end
if nargin < 11
    ny_p = 20;
end
if nargin < 12
    ny_b = 200;
end
if nargin < 13
    nt = 101;
end
%L_base = sqrt(ThermalDiffus_b*0.01*ny_b^2/19.6);
SaveOrNot = 0;               % if 1, save the calculating result
CompareAndDisplayOrNot = 0;         % if 1, do some comparing and showing work



% mesh dissection 
Image_Width = size(HeatFluxMatrixes, 2);
Image_Height = size(HeatFluxMatrixes, 1);
if(size(HeatFluxMatrixes, 3) ~= nt)
    error('The time length of the heat flux matrixes and nt not match, Please check that!');
end
dt = tmax/(nt - 1);

%Model_Width = 1*(dt*Image_Width);              % OK
%Model_Height = 1*(dt*Image_Height);
% Model_Width = 1.5*(dt*Image_Width);               % PVC+Al Mutation better
% Model_Height = 1.5*(dt*Image_Height);
 Model_Width = 3/4*(dt*Image_Width);              % Gaussian OK
 Model_Height = 3/4*(dt*Image_Height);
disp(sprintf('The Width and Height of the model are: %d*%d', Model_Width, Model_Height));
disp(sprintf('The Width and Height of the image are: %d*%d', Image_Width, Image_Height));

dx = Model_Width/(Image_Width-1);
dy = Model_Height/(Image_Height-1);
dz_p = L_polymer/(ny_p - 1);
dz_b = L_base/(ny_b - 1);


% initialization
% first, initialize the surface temperature matrixes
SurfaceTemperatureMatrixes = zeros(Image_Height, Image_Width, nt);
TemperatureMatrixes_p = zeros(Image_Height, Image_Width, ny_p, nt);
TemperatureMatrixes_b = zeros(Image_Height, Image_Width, ny_b, nt);
for x = 1:2
    for y = 1:Image_Height
        for i = 1:nt
            HFAtOnePoint(i) = HeatFluxMatrixes(y, x, i);
        end
%         [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS1( ThermalDiffus_p,...
%                 ThermalDiffus_b, HFAtOnePoint', ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt);
        [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS_FiniteBase_temp( ThermalDiffus_p,...
            ThermalDiffus_b, HFAtOnePoint', ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt, h_c);
        for i = 1:ny_p
            for j = 1:nt
                TemperatureMatrixes_p(y, x, i, j) = theata_p(i, j);
            end
        end
        for i = 1:ny_b
            for j = 1:nt
                TemperatureMatrixes_b(y, x, i, j) = theata_b(i, j);
            end
        end
    end
end
for x = 3:Image_Width
    for y = 1:2
        for i = 1:nt
            HFAtOnePoint(i) = HeatFluxMatrixes(y, x, i);
        end
%         [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS1( ThermalDiffus_p,...
%                 ThermalDiffus_b, HFAtOnePoint', ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt);
        [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS_FiniteBase_temp( ThermalDiffus_p,...
            ThermalDiffus_b, HFAtOnePoint', ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt, h_c);
        for i = 1:ny_p
            for j = 1:nt
                TemperatureMatrixes_p(y, x, i, j) = theata_p(i, j);
            end
        end
        for i = 1:ny_b
            for j = 1:nt
                TemperatureMatrixes_b(y, x, i, j) = theata_b(i, j);
            end
        end
    end
end

r_x_p = ThermalDiffus_p*dt/dx.^2;
r_x_b = ThermalDiffus_b*dt/dx.^2;
r_y_p = ThermalDiffus_p*dt/dy.^2;
r_y_b = ThermalDiffus_b*dt/dy.^2;
r_z_p = ThermalDiffus_p*dt/dz_p.^2;
r_z_b = ThermalDiffus_b*dt/dz_b.^2;

Omiga = (ThermalConduc_b*dz_p)/(ThermalConduc_p*dz_b);

%% computation
for x = 3:Image_Width
    disp(sprintf('The %dth column', x));
    for y = 3:Image_Height
        % Coefficient matrix construction
        Coef_matrix = zeros((ny_p+ny_b), (ny_p+ny_b));
        Coef_matrix(1, 1) = 1 + Omiga;         % first row
        Coef_matrix(1, 2) = -1;
        Coef_matrix(1, ny_p+2) = -Omiga;
        for i = 2:(ny_p-1)
            Coef_matrix(i, i-1) = -r_z_p;
            Coef_matrix(i, i) = 1 + 2*r_z_p - r_x_p - r_y_p;
            Coef_matrix(i, i+1) = -r_z_p;
        end
        Coef_matrix(ny_p, ny_p-1) = -1;
        Coef_matrix(ny_p, ny_p) = 1;
        Coef_matrix(ny_p+1, 2) = -1;
        Coef_matrix(ny_p+1, ny_p+1) = 1 + Omiga;
        Coef_matrix(ny_p+1, ny_p+2) = -Omiga;
        for i = (ny_p+2):(ny_p+ny_b-1)
            Coef_matrix(i, i-1) = -r_z_b;
            Coef_matrix(i, i) = 1 + 2*r_z_b - r_x_b - r_y_b;
            Coef_matrix(i, i+1) = -r_z_b;
        end
        Coef_matrix(ny_p+ny_b, ny_p+ny_b-1) = -ThermalConduc_b/dz_b;
        Coef_matrix(ny_p+ny_b, ny_p+ny_b) = h_c + ThermalConduc_b/dz_b;
        
        for k = 2:nt
            % update the right hand side
            b = zeros(ny_p+ny_b, 1);
            for i = 2:(ny_p-1)
                b(i) = -2*r_x_p*TemperatureMatrixes_p(y, x-1, i, k) + r_x_p*TemperatureMatrixes_p(y, x-2, i, k) - 2*r_y_p*TemperatureMatrixes_p(y-1, x, i, k) + r_y_p*TemperatureMatrixes_p(y-2, x, i, k) + TemperatureMatrixes_p(y, x, i, k-1);
            end
            b(ny_p) = dz_p*HeatFluxMatrixes(y, x, k)/ThermalConduc_p;
            for i = (ny_p+2):(ny_p+ny_b-1)
                b(i) = -2*r_x_b*TemperatureMatrixes_b(y, x-1, i-ny_p, k) + r_x_b*TemperatureMatrixes_b(y, x-2, i-ny_p, k) - 2*r_y_b*TemperatureMatrixes_b(y-1, x, i-ny_p, k) + r_y_b*TemperatureMatrixes_b(y-2, x, i-ny_p, k) + TemperatureMatrixes_b(y, x, i-ny_p, k-1);
            end

            theata = linsolve(Coef_matrix, b);
            for i = 1:ny_p
                TemperatureMatrixes_p(y, x, i, k) = theata(i);
            end
            for i = 1:ny_b
                TemperatureMatrixes_b(y, x, i, k) = theata(i+ny_p);
            end
        end
    end
end

disp(sprintf('Now constructing the surface temperature matrixes'));
for x = 1:Image_Width
    for y = 1:Image_Height
        for k = 1:nt
            SurfaceTemperatureMatrixes(y, x, k) = TemperatureMatrixes_p(y, x, end, k);
        end
    end
end

if(SaveOrNot)
    disp(sprintf('Now saving the results'));
    SaveDir = [GetAbsolutePath('../../Results/Simulation') '\Numerical_3D\PVC+Al\Gaussian\'];
    %SaveDir = [GetAbsolutePath('../../Results/Simulation') '\Numerical_3D\PVC+Al\Mutation\'];
    %SaveDir = [GetAbsolutePath('../../Results/Simulation') '\Numerical_3D\PVC+Nylon6\Gaussian\'];
    %SaveDir = [GetAbsolutePath('../../Results/Simulation') '\Numerical_3D\PVC+Nylon6\Mutation\'];
    MatSaveDir = [SaveDir 'Mat Files\'];
    if(not(exist(MatSaveDir)))
        mkdir(MatSaveDir);
    end
    ImageSaveDir = [SaveDir 'HFImageAtOneTimePoint\'];
    if(not(exist(ImageSaveDir)))
        mkdir(ImageSaveDir);
    end


    SaveMatName = sprintf('TemperatureResults3D_SimulationGaussianHF_PVC_Al_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults3D_SimulationMutationHF_PVC_Al_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults3D_SimulationGaussianHF_PVC_Nylon6_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults3D_SimulationMutationHF_PVC_Nylon6_%d_%d.mat',Image_Width, Image_Height);
    save([MatSaveDir SaveMatName], 'SurfaceTemperatureMatrixes', 'TemperatureMatrixes_p', 'TemperatureMatrixes_b');
end


%%
if(CompareAndDisplayOrNot)
    disp(sprintf('Now, comparing with the 1D method'));
    disp(sprintf('\tfirst, compute the 1D results'));
    SurfaceTemperatureMatrixes1D = zeros(Image_Height, Image_Width, nt);
    TemperatureMatrixes1D_p = zeros(Image_Height, Image_Width, ny_p, nt);
    TemperatureMatrixes1D_b = zeros(Image_Height, Image_Width, ny_b, nt);
    TemperatureMatrixes1D_p(1:2, 1:2, :, :) = TemperatureMatrixes_p(1:2, 1:2, :, :);
    TemperatureMatrixes1D_b(1:2, 1:2, :, :) = TemperatureMatrixes_b(1:2, 1:2, :, :);
    for x = 1:Image_Width
        for y = 1:Image_Height
            for i = 1:nt
                HFAtOnePoint(i) = HeatFluxMatrixes(y, x, i);
            end
            [theata_p, y_p, err_p, theata_b, y_b, err_b, t] = ForwardNumerical_BTCS1( ThermalDiffus_p,...
                ThermalDiffus_b, HFAtOnePoint', ThermalConduc_p, ThermalConduc_b, L_polymer, L_base, tmax, ny_p, ny_b, nt);
            for i = 1:ny_p
                for j = 1:nt
                    TemperatureMatrixes1D_p(y, x, i, j) = theata_p(i, j);
                end
            end
            for i = 1:ny_b
                for j = 1:nt
                    TemperatureMatrixes1D_b(y, x, i, j) = theata_b(i, j);
                end
            end
        end
    end
    for x = 1:Image_Width
        for y = 1:Image_Height
            for k = 1:nt
                SurfaceTemperatureMatrixes1D(y, x, k) = TemperatureMatrixes1D_p(y, x, end, k);
            end
        end
    end
    SaveMatName = sprintf('TemperatureResults1D_SimulationGaussianHF_PVC_Al_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults1D_SimulationMutationHF_PVC_Al_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults1D_SimulationGaussianHF_PVC_Nylon6_%d_%d.mat',Image_Width, Image_Height);
    %SaveMatName = sprintf('TemperatureResults1D_SimulationMutationHF_PVC_Nylon6_%d_%d.mat',Image_Width, Image_Height);
    save([MatSaveDir SaveMatName], 'SurfaceTemperatureMatrixes1D', 'TemperatureMatrixes1D_p', 'TemperatureMatrixes1D_b');
    
    disp(sprintf('\tShow some results'))
    disp(sprintf('\t\t1. Plot temperature results at some special points'));
    if(strcmp(SimulationMethod, 'Gaussian'))
        position_x = floor(Image_Width/2);
        position_y = floor(Image_Height/2);
        for k = 1:nt
            Theata_ps_1D(k) = SurfaceTemperatureMatrixes1D(position_y, position_x, k);
            Theata_ps_3D(k) = SurfaceTemperatureMatrixes(position_y, position_x, k);
        end
        figure; plot(t, Theata_ps_1D', '-b');
        xlabel('Time (Sec)');
        ylabel('Surface Temperature (K)');
        hold on; plot(t, Theata_ps_3D', '--r');
        legend('1D Simulated Temperature at the Peak Maximum', '3D Simulated Temperature at the Peak Maximum');
        
        for k = 1:nt
            Theata_ps_1D(k) = SurfaceTemperatureMatrixes1D(position_y, floor(Image_Width*7/10), k);
            Theata_ps_3D(k) = SurfaceTemperatureMatrixes(position_y, floor(Image_Width*7/10), k);
        end
        figure; plot(t, Theata_ps_1D', '-b');
        xlabel('Time (Sec)');
        ylabel('Surface Temperature (K)');
        hold on; plot(t, Theata_ps_3D', '--r');
        legend(sprintf('1D Simulated Temperature at the Point (%d, %d)',floor(Image_Width*7/10),position_y), sprintf('3D Simulated Temperature at the Point (%d, %d)',floor(Image_Width*7/10),position_y));
    end
    if(strcmp(SimulationMethod, 'Mutation'))
        position_x = floor(Image_Width/2);
        position_y = floor(Image_Height/2);
        for k = 1:nt
            Theata_ps_1D(k) = SurfaceTemperatureMatrixes1D(position_y, position_x, k);
            Theata_ps_3D(k) = SurfaceTemperatureMatrixes(position_y, position_x, k);
        end
        figure; plot(t, Theata_ps_1D', '-b');
        xlabel('Time (Sec)');
        ylabel('Surface Temperature (K)');
        hold on; plot(t, Theata_ps_3D', '--r');
        legend(sprintf('1D Simulated Temperature at the Point (%d,%d)', position_x, position_y), sprintf('3D Simulated Temperature at the Point (%d, %d)', position_x, position_y));
        
        for k = 1:nt
            Theata_ps_1D(k) = SurfaceTemperatureMatrixes1D(position_y, position_x-1, k);
            Theata_ps_3D(k) = SurfaceTemperatureMatrixes(position_y, position_x-1, k);
        end
        figure; plot(t, Theata_ps_1D', '-b');
        xlabel('Time (Sec)');
        ylabel('Surface Temperature (K)');
        hold on; plot(t, Theata_ps_3D', '--r');
        legend(sprintf('1D Simulated Temperature at the Point (%d, %d)',position_x-1,position_y), sprintf('3D Simulated Temperature at the Point (%d, %d)',position_x-1,position_y));
    end
    
    
    
    disp(sprintf('\t\t2. plot a centerline result at the middle time point'));
    Along_x = linspace(1, Image_Width, Image_Width);
    for x = 1:Image_Width
        Theata_ps_VecAlongX_1D(x) = SurfaceTemperatureMatrixes1D(position_y, x, (nt-5));
        Theata_ps_VecAlongX_3D(x) = SurfaceTemperatureMatrixes(position_y, x, (nt-5));
    end
    figure; plot(Along_x, Theata_ps_VecAlongX_1D, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Surface Temperature (K)');
    hold on; plot(Along_x, Theata_ps_VecAlongX_3D, '--r');
    legend(sprintf('1D Simulated Temperature Along the Center\n Line at the %dth Time Point (z=%d pixel)',(nt-5), position_y), sprintf('3D Simulated Temperature Along the Center\n Line at the %dth Time Point (z=%d pixel)',(nt-5), position_y));
    
    
    disp(sprintf('\t\t3. show average surface temperature image'))
    sum1D = 0;
    sum3D = 0;
    for i = 2:max(size(t))
        sum1D = sum1D + SurfaceTemperatureMatrixes1D(:, :, i);
        sum3D = sum3D + SurfaceTemperatureMatrixes(:, :, i);
    end
    AverageTPImage_1D = sum1D./(max(size(t))-1);
    AverageTPImage_3D = sum3D./(max(size(t))-1);
    
    SaveMatName = sprintf('TemperatureResultsAtOnePointandAverage_PVC_Al_%d_%d.mat',position_x, position_y);
    %SaveMatName = sprintf('TemperatureResultsAtOnePointandAverage_PVC_Nylon6_%d_%d.mat',position_x, position_y);
    save([MatSaveDir SaveMatName], 't', 'Theata_ps_1D', 'Theata_ps_3D', 'AverageTPImage_1D', 'AverageTPImage_3D');
    
    figure;
    lims = [0 5];                       % PVC+Al Gaussian
    %lims = [0 100];                      % PVC+Al Mutation
    %lims = [0 80];                        % PVC+Nylon6 Gaussian
    %lims = [0 1500];                      % PVC+Nylon6 Mutation
    imagesc(AverageTPImage_1D, lims);
    xlabel('x (pixel)');
    ylabel('z (pixel)');
    text(27, -1, 'K');
    colorbar;
    truesize;
    SaveImgName = sprintf('AverageTPImage_1D.jpg');
    saveas(gcf, [SaveDir SaveImgName], 'jpg');
    close;
    
    figure;
    lims = [0 5];                      % PVC+Al Gaussian
    %lims = [0 100];                     % PVC+Al Mutation
    %lims = [0 80];                        % PVC+Nylon6 Gaussian
    %lims = [0 1500];                      % PVC+Nylon6 Mutation
    imagesc(AverageTPImage_3D, lims);
    xlabel('x (pixel)');
    ylabel('z (pixel)');
    text(27, -1, 'K');
    colorbar;
    truesize;
    SaveImgName = sprintf('AverageTPImage_3D.jpg');
    saveas(gcf, [SaveDir SaveImgName], 'jpg');
    close;
    
    %%
    % from the 3D simulated surface temperature, we recover heat flux using our
    % analytical method and compare with the given heat flux
    ComputationalTime_Sum_Analytical = 0;
    ComputationalTime_Sum_Numerical = 0;
    IterTimes = 3;
    tolerance = sqrt(0.01^2*nt);
    for pixel_x = 1:Image_Width
        for pixel_z = 1:Image_Height
            disp(sprintf('\tthe (%d, %d) pixel',pixel_x, pixel_z));
            for k = 1:nt
                Theata_ps_3D(k) = SurfaceTemperatureMatrixes(pixel_z, pixel_x, k);
            end
            if(size(Theata_ps_3D, 1)<size(Theata_ps_3D, 2))
                Theata_ps_3D = Theata_ps_3D';
            end
            
            % Recover heat flux using our analytical method from the 3D surface temperature
            disp(sprintf('\t\tAnalytical recovering...'));
            ComputationalTime = cputime;
            [q_s, t] = HFComputation_GeneralCase( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
                ThermalConduc_b, L_polymer, Theata_ps_3D, t);
            ComputationalTime = cputime - ComputationalTime;
            ComputationalTime_Sum_Analytical = ComputationalTime_Sum_Analytical + ComputationalTime;
            
            % Heat flux matries/images
            HeatFluxMatrixes_Analytical(pixel_z, pixel_x, :) = q_s;
            SaveMatName = sprintf('HeatFluxAtOnePixel_%d_%d_Analytical.mat', pixel_x, pixel_z);
            save([MatSaveDir SaveMatName], 't', 'q_s', 'ComputationalTime');
            
            % Recover heat flux using our Numerical method from the 3D surface temperature
            disp(sprintf('\t\tNumerical recovering...'));
            ComputationalTime = cputime;
            [q_s_matrix, theata_ps_matrix, t] = HFComputation_Numerical(ThermalDiffus_p, ThermalDiffus_b,...
                ThermalConduc_p, ThermalConduc_b, Theata_ps_3D, L_polymer, L_base, tmax, ny_p, ny_b, nt, tolerance, IterTimes);
            q_s = q_s_matrix(:, end);
            ComputationalTime = cputime - ComputationalTime;
            ComputationalTime_Sum_Numerical = ComputationalTime_Sum_Numerical + ComputationalTime;
            
            % Heat flux matries/images
            HeatFluxMatrixes_Numerical(pixel_z, pixel_x, :) = q_s;
            SaveMatName = sprintf('HeatFluxAtOnePixel_%d_%d_Numerical.mat', pixel_x, pixel_z);
            save([MatSaveDir SaveMatName], 't', 'q_s', 'ComputationalTime');
        end
    end
    HeatFluxMatrixes_Simulation = HeatFluxMatrixes;
    
    % average HF image
    Time_begin = 2;
    Time_end = max(size(t));
    sum = 0;
    for i = Time_begin:Time_end
        sum = sum + HeatFluxMatrixes_Simulation(:, :, i);
    end
    AverageImage_Simulation = sum./(Time_end-Time_begin+1);
    sum = 0;
    for i = Time_begin:Time_end
        sum = sum + HeatFluxMatrixes_Analytical(:, :, i);
    end
    AverageImage_Analytical = sum./(Time_end-Time_begin+1);
    sum = 0;
    for i = Time_begin:Time_end
        sum = sum + HeatFluxMatrixes_Numerical(:, :, i);
    end
    AverageImage_Numerical = sum./(Time_end-Time_begin+1);
    % AverageCenterLine
    AverageCenterLine_Simulation = AverageImage_Simulation(position_y, :);
    AverageCenterLine_Analytical = AverageImage_Analytical(position_y, :);
    AverageCenterLine_Numerical = AverageImage_Numerical(position_y, :);
    
    SaveMatName = 'HFAnalyticalAndNumericalFrom3DTemperature.mat';
    save([MatSaveDir SaveMatName], 't', 'HeatFluxMatrixes_Analytical', 'HeatFluxMatrixes_Numerical', 'HeatFluxMatrixes_Simulation', 'AverageImage_Simulation', 'AverageImage_Analytical', 'AverageImage_Numerical', 'AverageCenterLine_Simulation', 'AverageCenterLine_Analytical', 'AverageCenterLine_Numerical', 'ComputationalTime_Sum_Analytical', 'ComputationalTime_Sum_Numerical');
    
    % show and save the heat flux images
    for i = 1:max(size(t))
        if i<=10
            lims = [0 2000];
        else
            lims = [0 50000];
        end
        
        figure;
        imagesc(HeatFluxMatrixes_Analytical(:,:,i), lims);
        xlabel('x (pixel)');
        ylabel('z (pixel)');
        text(25, 20, 'W/m^2');
        colorbar;
        truesize;
        SaveImgName = sprintf('Analytical_HFImage_%d.jpg', i);
        saveas(gcf, [ImageSaveDir SaveImgName], 'jpg');
        close;
        
        figure;
        imagesc(HeatFluxMatrixes_Numerical(:,:,i), lims);
        xlabel('x (pixel)');
        ylabel('z (pixel)');
        text(25, 20, 'W/m^2');
        colorbar;
        truesize;
        SaveImgName = sprintf('Numerical_HFImage_%d.jpg', i);
        saveas(gcf, [ImageSaveDir SaveImgName], 'jpg');
        close;
        
        figure;
        imagesc(HeatFluxMatrixes_Simulation(:,:,i), lims);
        xlabel('x (pixel)');
        ylabel('z (pixel)');
        text(25, 20, 'W/m^2');
        colorbar;
        truesize;
        SaveImgName = sprintf('Simulation_HFImage_%d.jpg', i);
        saveas(gcf, [ImageSaveDir SaveImgName], 'jpg');
        close;
    end
    
    % show the average image
    figure;
    imagesc(AverageImage_Simulation, lims);
    xlabel('x (pixel)');
    ylabel('z (pixel)');
    text(25, 20, 'W/m^2');
    colorbar;
    truesize;
    SaveImgName = sprintf('Simulation_AverageHFImage.jpg');
    saveas(gcf, [SaveDir SaveImgName], 'jpg');
    close;
    
    figure;
    imagesc(AverageImage_Analytical, lims);
    xlabel('x (pixel)');
    ylabel('z (pixel)');
    text(25, 20, 'W/m^2');
    colorbar;
    truesize;
    SaveImgName = sprintf('Analytical_AverageHFImage.jpg');
    saveas(gcf, [SaveDir SaveImgName], 'jpg');
    close;
    
    figure;
    imagesc(AverageImage_Numerical, lims);
    xlabel('x (pixel)');
    ylabel('z (pixel)');
    text(25, 20, 'W/m^2');
    colorbar;
    truesize;
    SaveImgName = sprintf('Numerical_AverageHFImage.jpg');
    saveas(gcf, [SaveDir SaveImgName], 'jpg');
    close;
    
    % plot average center line
    figure; plot([1:1:Image_Width], AverageCenterLine_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], AverageCenterLine_Numerical, '--g');
    legend(sprintf('Simulated Heat Flux(Time-Averaged,z=%d pixel)', position_y), sprintf('Numerical Recovered Heat Flux from the 3D Simulated\n Temperature(Time-Averaged,z=%d pixel)', position_y));
    
    figure; plot([1:1:Image_Width], AverageCenterLine_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], AverageCenterLine_Analytical, '-.r');
    legend(sprintf('Simulated Heat Flux(Time-Averaged,z=%d pixel)', position_y), sprintf('Analytical Recovered Heat Flux from the 3D Simulated\n Temperature(Time-Averaged,z=%d pixel)', position_y));
    
    figure; plot([1:1:Image_Width], AverageCenterLine_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], AverageCenterLine_Analytical, '-.r');
    hold on; plot([1:1:Image_Width], AverageCenterLine_Numerical, '--g')
    legend(sprintf('Simulated Heat Flux(Time-Averaged,z=%d pixel)', position_y), sprintf('Analytical Recovered Heat Flux from the 3D Simulated\n Temperature(Time-Averaged,z=%d pixel)', position_y), sprintf('Numerical Recovered Heat Flux from the 3D Simulated\n Temperature(Time-Averaged,z=%d pixel)', position_y));
    
    % plot the center line at the 96th time point
    for i = 1:Image_Width
        HeatFlux_VecAlongX_Numerical(i) = HeatFluxMatrixes_Numerical(position_y, i, (nt-5));
        HeatFlux_VecAlongX_Analytical(i) = HeatFluxMatrixes_Analytical(position_y, i, (nt-5));
        HeatFlux_VecAlongX_Simulation(i) = HeatFluxMatrixes_Simulation(position_y, i, (nt-5));
    end
    figure; plot([1:1:Image_Width], HeatFlux_VecAlongX_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], HeatFlux_VecAlongX_Numerical, '--g');
    legend(sprintf('Simulated Heat Flux at the %dth Time Point (z=%d pixel)',(nt-5), position_y), sprintf('Numerical Recovered Heat Flux from the 3D Simulated\n Temperature at the %dth Time Point (z=%d pixel)',(nt-5), position_y));
    
    figure; plot([1:1:Image_Width], HeatFlux_VecAlongX_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], HeatFlux_VecAlongX_Analytical, '-.r');
    legend(sprintf('Simulated Heat Flux at the %dth Time Point (z=%d pixel)',(nt-5), position_y), sprintf('Analytical Recovered Heat Flux from the 3D Simulated\n Temperature at the %dth Time Point (z=%d pixel)',(nt-5), position_y));
    
    figure; plot([1:1:Image_Width], HeatFlux_VecAlongX_Simulation, '-b');
    xlabel('Centerline (pixel)');
    ylabel('Heat Flux (W/m^2)');
    hold on; plot([1:1:Image_Width], HeatFlux_VecAlongX_Analytical, '-.r');
    hold on; plot([1:1:Image_Width], HeatFlux_VecAlongX_Numerical, '--g');
    legend(sprintf('Simulated Heat Flux at the %dth Time Point (z=%d pixel)',(nt-5), position_y), sprintf('Analytical Recovered Heat Flux from the 3D Simulated\n Temperature at the %dth Time Point (z=%d pixel)',(nt-5), position_y), sprintf('Numerical Recovered Heat Flux from the 3D Simulated\n Temperature at the %dth Time Point (z=%d pixel)',(nt-5), position_y));
end

