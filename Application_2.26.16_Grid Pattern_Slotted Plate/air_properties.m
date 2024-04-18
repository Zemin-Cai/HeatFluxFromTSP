function  [eta lambda Cp Pr rhoSI]=air_properties(T,p)
%  T in degrees K
%   p in Pa
%
%Viscosity and Thermal Conductivity Equations for
%Nitrogen, Oxygen, Argon, and Air
%E. W. Lemmon1, 2 and R. T Jacobsen3
%International Journal of Thermophysics, Vol. 25, No. 1, January 2004 (© 2004)

% Table I
Tc=132.6312;%(K)
rho_c=10.4477;%(mol dm -3)
Pc=3.78502;%(MPa)
M=28.9586;%(g·mol-1)
e_k=103.3;%(K)
sigma=0.36;%(nm)
Xio=0.11;%(nm)
Gamma=0.055;
qD=0.31;%(nm)
Tref=265.262;%(K)

Runiversal=8.31451; % J • mol–1 • K–1.
R=Runiversal/M*1000;

% Table 2
bi=[
0.431
-0.4623
0.08406
0.005341
-0.00331];
%
%Table 3
Ni=[10.72 1.122 0.002019 -8.876 -0.02916];
ti=[0.2 0.05 2.4 0.6 3.6];
di=[1 4 9 1 8];
li=[0 0 0 1 1];
%
% Table 4
Ni2=[1.308 1.405 -1.036 8.743 14.76 -16.62 3.793 -6.142 -0.3778];
ti2=[0.0 -1.1 -0.3 0.1 0 0.5 2.7 0.3 1.3];
di2=[0 0 0 1 2 3 7 7 11];
li2=[0 0 0 0 0 2 2 2 2];
%

rho=p/(R*T)/M;
rhoSI=rho*M;
%rho	=0.041623719;%	mol/dm^3(mol/L)  or 1.205364631	Kg/m^3

tau=Tc/T;
delta=rho/rho_c;
Tstar=T/e_k;

sum=0;
for i=1:5
    sum=sum+bi(i)*(log(Tstar))^(i-1);
end
omega=exp(sum);
eta0=0.0266958*sqrt(M*T)/(sigma^2*omega);

etaR=0;
for i=1:5
    etaR=etaR+Ni(i)*tau^ti(i)*delta^di(i)*exp(-li(i)*delta^li(i));
end

eta=eta0+etaR;  % viscosity micro Pa-sec

% Thermal Conductivity

lambda0=Ni2(1)*eta0+Ni2(2)*tau^ti2(2)+Ni2(3)*tau^ti2(3);

lambdaR=0;
for i=4:9
    lambdaR=lambdaR+Ni2(i)*tau^ti2(i)*delta^di2(i)*exp(-li2(i)*delta^li2(i));
end

lambda=lambda0+lambdaR;

%
% Cp

%Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon,
%and Oxygen From 60 to 2000 K at Pressures to 2000 MPa
%Eric W. Lemmona?
%J. Phys. Chem. Ref. Data, Vol. 29, No. 3, 2000
% Table 12


Ni3=[3.490888032
2.395525583E-6
7.172111248E-9
-3.115413101E-13
0.223806688
0.791309509
0.212236768
0.197938904
3364.011
2242.45
11580.4];

sum=0.0;
for i=1:4
    sum=sum + Ni3(i)*T^(i-1);
end
u=Ni3(9)/T;
v=Ni3(10)/T;
w=Ni3(11)/T;

Cp=R*(sum+Ni3(5)/T^(3/2)+Ni3(6)*u^2*exp(u)/(exp(u)-1)^2+Ni3(7)*v^2*exp(v)/(exp(v)-1)^2)+2/3*Ni3(8)*w^2*exp(-w)/(2/3*exp(-w)+1)^2;

Pr=Cp*eta/lambda/1000;





