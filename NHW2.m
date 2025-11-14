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
u_inf3 = 5; %m/s

r = 2e-3; %Radius of wire, m
Ac = pi*r^2; %Cross sectional area, m^2
P = 2*pi*r;  %Perimeter, m
u = .015; %Speed of wire, m/s

%MATERIAL PROPERTIES
rho1 = 8960; %kg/m^3
k1 = 385;    %W/mK
c1 = 445;    %J/kgK
alpha1 = k1/rho1/c1;    %Thermal diffusivity
%AIR PROPERTIES
rho_a = 1.225;
kf = .025;
mu = 1.81e-5; %kg/ms
nu = mu/rho1;
%==========================================================================
%%
%Now to find the convection coefficients in (1) and (3)
Re1 = rho_a*u_inf1*2*r/mu;
Re3 = rho_a*u_inf3*2*r/mu;

Pr1 = .707; %Standard Pr for air
Pr3 = .707;

Nu1 = .3 + (.62*Re1^.5*Pr1^(1/3)) / (1+(.4/Pr1)^(2/3))^1/4 * (1+(Re1/282000)^(5/8))^.8;
Nu3 = .3 + (.62*Re3^.5*Pr3^(1/3)) / (1+(.4/Pr3)^(2/3))^1/4 * (1+(Re3/282000)^(5/8))^.8;

h1 = Nu1/2/r*kf;
h3 = Nu3/2/r*kf;

%Now calculating lambda values for (1) and (3)
m1 = sqrt(h1*P/k1/Ac);
m3 = sqrt(h3*P/k1/Ac);

lam11 = .5*( (u_inf1/alpha1) + sqrt((u_inf1/alpha1)^2 + 4*m1^2));
lam12 = .5*( (u_inf1/alpha1) - sqrt((u_inf1/alpha1)^2 + 4*m1^2));

lam21 = .5*( (u_inf3/alpha1) + sqrt((u_inf3/alpha1)^2 + 4*m3^2));
lam22 = .5*( (u_inf3/alpha1) - sqrt((u_inf3/alpha1)^2 + 4*m3^2));

%Iterate through guessed T_inf1 values until the final boundary condition
%is satisfied


T_inf1g = 850;%K
theta1g = T1k - T_inf1g;

options = optimoptions('fsolve','Display','iter');
Cg = [100;-5];
%C = fzero(@(C)[C(1)+C(2)-theta1g;C(1)*lam11*exp(lam11*L12)+C(2)*lam12*exp(lam12*L12)],Cg);
%C1 = fzero(@(C1,C2)[C1+C2-theta1g;C1*lam11*exp(lam11*L12)+C2*lam12*exp(lam12*L12)],Cg);





%C1 = 10;
%C2 = T1k-T_inf1g-C1;
%the1 = fsolve(@(C) theta1(C,lam11,lam12,L12,theta1g),[C1,C2]);


%==========================================================================
%%
%Numerical Solution
%Defining initial variables

dx1 = .1/100; %cm per node in (1) converted to per meter
dx2 = 1/100;
dx3 = .1/100;

n1 = linspace(0,L12,L12/dx1); %Vector containing each position in (1) in m
n2 = linspace(0,L23,L23/dx2);
%No n3 here because length is unknown

T1 = zeros(1,length(n1)); %Temperature for each node in (1)

%Solving for the temperature in section 1
l1 = u_inf1*dx1/2/alpha1;

A1 = zeros(length(n1)); %Matrix multiplying the temperature
A1(1,1) = 1;
A1(length(n1),length(n1)) = 1;
A1(length(n1),length(n1)-1) = -2;
A1(length(n1),length(n1)-2) = 1;

%Creating the A matrix for interior nodes, the one multiplying the temperatures
for i = 2:length(n1)-1
    A1(i,i-1) = 1-l1;   %NOTE THESE VALUES MAY NOT BE CORRECT
    A1(i,i) = m1^2*dx1^2;
    A1(i,i+1) = 1+l1;
end

T_pos2 = 850;
while abs(T_pos2 - T2k) > .5
%Creating the C matrix, the one on the other side of A and T
theta1g = T1k - T_inf1g;
C1 = zeros(length(n1),1);
C1(1) = theta1g;

theta1 = A1\C1;
T1 = theta1 + T_inf1g;

figure(1)
plot(n1,T1)
xlabel('Distance (m)')
ylabel('Temperature (K)')
title('Section 1 Temperature Distribution')

T_pos2 = T1(length(n1));
if T_pos2 < T2k
    T_inf1g = T_inf1g + 1;
else
    T_inf1g = T_inf1g - 1;
end
end

%==========================================================================
%%
%Section 1 now solved and plotted in figure 1
%Now we plot section 2
T2 = linspace(T_pos2,T_pos2,length(n2));
T2 = T2.';
figure(2)
plot(n2,T2)
xlabel('Distance (m)')
ylabel('Temperature (K)')
title('Section 2 Temperature Distribution')

%==========================================================================
%%
%Section 2 now solved and plotted in figure 2
%Now we solve and plot section 3

L34 = 1; %Initial guess for the length of section 3
n3 = linspace(0,L34,L34/dx3);
T3 = zeros(1,length(n3));

while abs(T3(length(T3)) - T4k) > 2
n3 = linspace(0,L34,L34/dx3);

l3 = u_inf3*dx3/2/alpha1;

A3 = zeros(length(n3)); %Matrix multiplying the temperature
A3(1,1) = 1;
A3(length(n3),length(n3)) = 1;
A3(length(n3),length(n3)-1) = -2;
A3(length(n3),length(n3)-2) = 1;

%Creating the A matrix for interior nodes, the one multiplying the temperatures
for i = 2:length(n3)-1
    A3(i,i-1) = 1-l3;   %NOTE THESE VALUES MAY NOT BE CORRECT
    A3(i,i) = m3^2*dx3^2;
    A3(i,i+1) = 1+l3;
end%of A matrix for loop

%Creating the C matrix for the third section
theta3 = T1(length(n1)) - T_inf2;
C3 = zeros(length(n3),1);
C3(1) = theta3;

%Calculating temperature in section 3
theta3 = A3\C3;
T3 = theta3 + T_inf2;

figure(3)
plot(n3,T3)
xlabel('Distance (m)')
ylabel('Temperature (K)')
title('Section 3 Temperature Distribution')

if T3(length(T3)) > T4k
    L34 = L34 + .1;
else
    L34 = L34 - .1;
end
end%of while statement

%==========================================================================
%%
%Now that all three sections have been found, we plot them together
n = horzcat(n1,n2+n1(length(n1)),n3+n1(length(n1))+n2(length(n2)));
T = horzcat(T1.',T2.',T3.');
figure(4)
plot(n,T,'r')
xlabel('Distance (m)')
ylabel('Temperature (K)')
title('All Sections Temperature Distribution')





