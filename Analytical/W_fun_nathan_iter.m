function [W, t_vec] = W_fun_nathan_iter( ThermalDiffus_p, ThermalConduc_p, ThermalDiffus_b,...
   ThermalConduc_b, L_polymer, L_base, t_vec, h_c)

%This function calculate the integral for W cited in the paper "Analytical
%heat transfer method ....

% function [W] = IntegrandInverse(t)


% h_c = 0.001;
% k_p = .15; %Conductivity of the polymer
% rho_p = 1420; %Density of the polymer
% C_p = 1090; %Specific heat of the polymer
% a_p = k_p/(C_p*rho_p); %Difusivity of the polymer
% k_b = 204; %Conductivity of the aluminum base
% rho_b = 2700; %Density of the aluminum base
% C_b = 904; %Specific heat of the aluminum base
% a_b = k_b/(C_b*rho_b); %Difusivity of the aluminum baser
% ep = sqrt(k_p*rho_p*C_p/(k_b*rho_b*C_b));
% ep_1 = (1-ep)/(1+ep);
% h_c=0;
% L_b = 12.7e-3; % Thickness of the aluminum base
% L_p = 20e-6; % Thickness of the polymer
% delta = pi()/480;
% 
% 
% tmax = 3;
% nt = 301;
% 
% HeatfluxFunc_index = 2;                   % 1 - 4
% switch HeatfluxFunc_index
%     case 1
%         HeatfluxFunc = '(t>=0).*(25000) + (t<0).*(0)';
%     case 2
%         HeatfluxFunc = '(t>=0 & t<=1.5).*(25000) +(t>1.5).*(50000) + (t<0).*(0)';
%     case 3
%         HeatfluxFunc = '(t>0 & t<=1.5).*(50000/1.5.*t) + (t>1.5).*((50000/1.5).*(3-t)) + (t<=0).*(0)';
%     case 4
%         HeatfluxFunc = '(t>=0 & t<=1).*(25000*(2*t-t.^2)) + (t>1).*(25000) + (t<0).*(0)';
%         tmax = 5;
%         nt = 501;    
% end
% t_vec = linspace(0, tmax, nt)';



k_p = ThermalConduc_p ; %Conductivity of the polymer
a_p = ThermalDiffus_p;  %Difusivity of the polymer
k_b = ThermalConduc_b; %Conductivity of the aluminum base
a_b = ThermalDiffus_b; %Difusivity of the aluminum baser

ep = (ThermalConduc_p.*sqrt(ThermalDiffus_b))./(ThermalConduc_b.*sqrt(ThermalDiffus_p));
ep_1 = (1 - ep)./(1 + ep);

h_c_head = h_c*sqrt(ThermalDiffus_b)./ThermalConduc_b;
h_c = h_c_head;

L_b = L_base; % Thickness of the aluminum base
L_p = L_polymer; % Thickness of the polymer

delta = pi()/480;

for i=1:length(t_vec)
    t=t_vec(i);

    F1 = 1i*h_c*exp(1i*delta*.5)*t.^.5; %for m

    F2 = -2*L_b*1i*exp(-1i*delta)./(t*a_b).^.5; %for m

    F3 = -2*L_p*cos((pi()-delta)*.5)./(t*a_p(i)).^.5; %for E

    F4 = -2*L_p*sin(.5*(pi()-delta))./(t*a_p(i)).^.5; %for alpha



    eqn = @(Z) (1 - (abs((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z.*F2))).*exp(Z.*F3)).^2+...
        2*abs((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z.*F2))).*exp(Z.*F3).*...
        sin(Z.*F4+angle((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z.*F2))))*1i)./...
        (1 + (abs((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z.*F2))).*exp(Z.*F3)).^2-...
        2*abs((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z*F2))).*exp(Z.*F3).*...
        cos(Z.*F4+angle((ep_1(i)-(Z+F1)./(Z-F1).*exp(Z.*F2))./...
        (1-ep_1(i).*(Z+F1)./(Z-F1).*exp(Z.*F2))))).*...
        exp(-Z.^2*exp(-1i*delta)-1i*delta*.5);

        W(i) = 2.0/pi()^.5*real(integral(eqn,0,inf));
end
W(1)=0;

W=W';


