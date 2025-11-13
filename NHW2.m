%{
MAE 131A Numerical Homework 2
Written by Jake Sudduth and Juni Mireles, 11/04/2025
<=========================================================================>
Analytical and numerical solutions for a cylinder of copper wire drawn
through a die. There are three sections:
(1) Cooling before annealing. The wire comes out at a known temperature and
is cooled by convection to another known temperature. We must find the
temperature of the air to cool the copper if its moving at a known rate.
(2) Adiabatic enamel. The boundary conditions show that the temperature is
actually constant in this section.
(3) Cooling to room temperature. The wire is again cooled by convection to
another known temperature. Because the temperature is constant is known
after 2, the initial boundary condition is also known here.
<=========================================================================>
%}
clc; close all; clear all;

%Defining all known values
T1k = 900+273; %K
T2k = 582+273; %K
T4k = 40+273; %K
T_inf2 = 20+273; %K

L12 = .5; %m
L23 = .3; %m

u_inf1 = 10; %m/s
u_inf2 = 5; %m/s

r = 2e-3; %Radius of wire, m
Ac = pi*r^2; %Cross sectional area, m^2
P = 2*pi*r;  %Perimeter, m
u = .015; %Speed of wire, m/s

%MATERIAL PROPERTIES
rho = 8960; %kg/m^3
k = 385;    %W/mK
c = 445;    %J/kgK
alpha = k/rho/c;    %Thermal diffusivity
alpha = 1.11*10^-4;
%AIR PROPERTIES
rho_a = 1.225;
kf = .025;
mu = 1.81e-5; %kg/ms
nu = mu/rho;
%==========================================================================
%%
%Now to find the convection coefficients in (1) and (3)
Re1 = rho_a*u_inf1*2*r/mu;
Re2 = rho_a*u_inf2*2*r/mu;

Pr1 = .707; %Standard Pr for air
Pr2 = .707;

Nu1 = .3 + (.62*Re1^.5*Pr1^(1/3)) / (1+(.4/Pr1)^(2/3))^1/4 * (1+(Re1/282000)^(5/8))^.8;
Nu2 = .3 + (.62*Re2^.5*Pr2^(1/3)) / (1+(.4/Pr2)^(2/3))^1/4 * (1+(Re2/282000)^(5/8))^.8;

h1 = Nu1/2/r*kf;
h2 = Nu2/2/r*kf;

%Now calculating lambda values for (1) and (3)
m1 = sqrt(h1*P/k/Ac);
m2 = sqrt(h2*P/k/Ac);

lam11 = .5*( (u_inf1/alpha) + sqrt((u_inf1/alpha)^2 + 4*m1^2));
lam12 = .5*( (u_inf1/alpha) - sqrt((u_inf1/alpha)^2 + 4*m1^2));

lam21 = .5*( (u_inf2/alpha) + sqrt((u_inf2/alpha)^2 + 4*m2^2));
lam22 = .5*( (u_inf2/alpha) - sqrt((u_inf2/alpha)^2 + 4*m2^2));

%Iterate through guessed T_inf1 values until the final boundary condition
%is satisfied





%==========================================================================
%%
%Numerical Solution
%Defining initial variables

dx1 = .01; %m per node in (1)
dx2 = .01;
dx3 = .01;

n1 = linspace(0,L12,L12/dx1); %Vector containing each position in (1) in m
n2 = linspace(0,L23,L23/dx2);
%No n3 here because length is unknown

T1 = zeros(1,length(n1)); %Temperature for each node in (1)
T2 = zeros(1,length(n2));


T_inf1g = 20+273;%K
theta1g = T1k - T_inf1g;

options = optimoptions('fsolve','Display','iter');
Cg = [100;-5];
%C = fzero(@(C)[C(1)+C(2)-theta1g;C(1)*lam11*exp(lam11*L12)+C(2)*lam12*exp(lam12*L12)],Cg);
C = fzero(@(C1,C2)[C1+C2-theta1g;C1*lam11*exp(lam11*L12)+C2*lam12*exp(lam12*L12)],Cg);





%{
C1 = 10;
C2 = T1k-T_inf1g-C1;
the1 = fsolve(@(C) theta1(C,lam11,lam12,L12,theta1g),[C1,C2]);







%n3 = npm1*L34_g;


function F = theta1(C,lam11,lam12,L12,theta1g)
    F(1) = C(1) + C(2) - theta1g;
    F(2) = lam11*exp(lam11*L12)*C(1) + lam12*exp(lam12*L12)*C(2);
end

%}

