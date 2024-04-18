%% This program calculates the heat flux by using the analytical inverse solution
%% for a finite base. in three cases:
%% (1) at a selected point, (2) alone a horizontal line, and (3) in an image
%% 


clear all
close all

addpath(GetAbsolutePath('Core Funcs'));
addpath(GetAbsolutePath('Analytical'));
addpath(GetAbsolutePath('Numerical'));



%% read surface temperature files and form a 3D matrix
No_images=25;
for i=1:No_images
    file_name=strcat('dT_heated4_',num2str(i),'.dat');
    X =load(file_name);    
    dT=single(X);

    dT_matrix(:,:,i)= dT;
end

%% show typical temperature image
figure(1);
a=[-20 5];
imagesc(dT_matrix(:,:,25),a);
colormap('jet');
colorbar;
axis image;
xlabel('x (pixels)');
ylabel('y (pixels)');
title('Typical Temperature Change (K), t = 1.6 s in Run 11');



%% define time series
minframe=1;
maxframe=25;
framerate=5; % f/s
maxtime=maxframe/framerate;
mintime=minframe/framerate;
t=linspace(mintime-mintime,maxtime-mintime,maxframe-minframe+1);


%% give the thermal properties of the polymer and base material
ThermalDiffus_p = 9.7*10^(-8);          % polymer layer Mylar, m^2/s
ThermalConduc_p = 0.15;                 % W/m-K

ThermalDiffus_b = 6.903*10^(-5);            % Al 6061 alloy base
ThermalConduc_b = 167;


%% give the thicknesses of the polymer coating and base, and heat transfer coefficient at the base bottom 
 L_p = 120e-6;      % polymer coating thickness (m)
 L_b=9.525e-3;      % base thickness (m) 
 h_c=0;             % heat transfer coefficient at the base bottom


 %% give relevant temperatues 
Tb=43.4+273.15;         %  initial base temperature (K) for TSP measurements
T_off=43.4+273.15;      % wind-off temperature (K) in case 3, which is equal to Tb
T_aw=317.1;             % adiabatic wall temperature (K) determined for Case 3

%% jet parameters
 D=6.35*10^(-3);        % jet diameter (m)
 k_air=0.0257;          % air thermal conductivity (W/m-K)

 
 %% calulate the heat flux in the following cases:
 %% (1) at a selected point, (2) alone a horizontal line, and (3) in an image
 
N_interp=200;
[W_Fun, t] = W_fun_nathan( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
              ThermalConduc_b, L_p, L_b, t, h_c);
 
 [row,col,nt]=size(dT_matrix);
 
 for case_index = 2:2
        switch case_index        
        
            case 1
                
            %% calulate the heat flux at a selected point 
            xc=247;
            yc=99;

            Theata_ps=squeeze(dT_matrix(yc,xc,:));

            a = 1;
            b = [1/4 1/4 1/4 1/4];
            Theata_ps = filter(b,a,Theata_ps);
                                          
            [qs, t, Delta] = HF_Analy_FiniteBase_W1( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
                      ThermalConduc_b, L_p, L_b, Theata_ps, t, h_c, N_interp, W_Fun);

            %% calculate surface temperature and Nusselt number
            qs=qs';
            Ts=Theata_ps+T_off;     % surface temperature (K)
            h=qs./(Ts-T_aw);       % heat transfer coefficient
            Nu=h*D/k_air;          % Nusselt number
            
                      
            %% show the time history of temperature at that point
            figure(2); 
            plot(t,Theata_ps,'o-');
            xlabel('Time (s)');
            ylabel('Temperature Change (K)');
            grid;
            title('History of Temperature Change at a Point (K)'); 
            
            figure(3);
            plot(t,qs,'o-');
            xlabel('Time (s)');
            ylabel('Heat Flux (W/m^2)');
            grid;
            title('Surface Heat Flux (W/m^2)');

            figure(4);
            plot(t,Nu,'o-');
            xlabel('Time (s)');
            ylabel('Nu');
            grid;
            title('Nusselt Number');            
            axis([0 5 0 300]);
            
            
            
            case 2
        
            %% calculate heat flux distribution along the x-axis at a
            %% selcted location
      
            yc=99;
            j=1;
            for k=1:5:col
                Theata_ps=squeeze(dT_matrix(yc,k,:));

                a = 1;
                b = [1/4 1/4 1/4 1/4];
                Theata_ps = filter(b,a,Theata_ps); 
               
                [qs, t, Delta] = HF_Analy_FiniteBase_W1( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
                         ThermalConduc_b, L_p, L_b, Theata_ps, t, h_c, N_interp, W_Fun);
                
                qs=qs';
                Ts=Theata_ps+T_off;
                h=qs./(Ts-T_aw);
                Nu=h*D/k_air;

                qs_x(j)=mean(qs(15:end)); % avereged value over a timespan after the transient stage
                Nu_x(j)=mean(Nu(15:end));

                j=j+1
            end

            x=[1:length(qs_x)];
            figure(5);
            plot(x,qs_x,'-r');
            xlabel('x');
            ylabel('q (W/m^2)');
            legend('Analytical Solution');
            grid;

            x=[1:length(Nu_x)];
            figure(6);
            plot(x,Nu_x,'-r');
            xlabel('x');
            ylabel('Nu)');
            legend('Analytical Solution');
            grid;  
 
            
            case 3
        
            %% calculate heat flux distribution in an image using the 1D analytical
            %% inverse method
  
            i=1;
            for n=1:1:row
                j=1;
                for k=1:1:col
                    Theata_ps=squeeze(dT_matrix(n,k,:));

                    a = 1;
                    b = [1/4 1/4 1/4 1/4];
                    Theata_ps = filter(b,a,Theata_ps); 

                    [qs, t, Delta] = HF_Analy_FiniteBase_W1( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
                             ThermalConduc_b, L_p, L_b, Theata_ps, t, h_c, N_interp, W_Fun);
                    
                    qs_matrix(i,j,:)=qs;
                    j=j+1
                end
                i=i+1
                [i j]            
                
            end
                                    
        end
 end

%% save the calculated heat flux data in sequence
for i=1:201
    q2=qs_matrix(:,:,i);
    figure(20);
    a=[0 30000];
    imagesc(q2,a);
    colormap('jet');
    colorbar;
    axis image;
    
    file_output=strcat('q_cal_cone_',num2str(i),'.dat');
    %dlmwrite(file_output,q2);  
end




















