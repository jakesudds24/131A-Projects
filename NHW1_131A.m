%{
MAE 131A Numerical Homework 1
Written by Jake Sudduth and Juni Mireles, 10/16/2025
<=========================================================================>
Creates a 1D mesh and solves the heat equation for a hollow cylinder with
an outer shell. Heat flux into the center of the tube, heat generation in
the mid portion, and a thin outer metal layer with convection boundary.
Compares with analytical solution and calculates error depending on the
number of nodes
<=========================================================================>
%}
clc; close all; clear all;

%Defining constants
r1 = 5e-3;
r2 = 10e-3;
r3 = 11e-3;
%Inner, outer, and metal radii
k = 1.9;
kg = 25;
%Thermal Conductivity of mid portion and metal
Tinf = 180+273;%K
h = 2000;     %W/m^2/K
qdot = 1.4e6;  %W/m^3
q_in = -2000;  %W/m^2
%Convection temp and coeff, volumetric heat gen in tube, and heat flux to
%inner portion of the tube all in SI units.

%First model is run for the mid portion of the tube, for which the heat
%flux into the center is known and the temperature at the interface with
%the metal is taken to be T2

%Defining number of nodes
n_mid = 10;
%Known temperature at the interface of middle and metal
T2 = 182.1 + 273;
%Known temperature at the metal and convection interface
T3 = 181.93 + 273;
%Distance between nodes
dr = (r2-r1)/n_mid;
%Creating nodes for metal
n_met = (r3-r2)/dr;
%Creating a vector of each node location
r = linspace(r1,r2,n_mid);

%%
%For interior nodes,
n_unk = n_mid - 2;

%Simplifying variable
S = -qdot/k*dr^2;

%Creating A and C matrices for solving equations, [A][T]=[C]
%Check the hand work to understand where these come from
A = zeros(n_mid);
C = zeros(n_mid,1);

%Creating the A matrix for interior nodes, the one multiplying the temperatures
for i = 2:n_mid - 1
    A(i,i) = -2;
    A(i,i-1) = (1 - dr / (2*r(i)) );
    A(i,i+1) = (1 + dr / (2*r(i)) );
end

%Creating the C matrix, the right hand side of the equation
for i = 1:n_mid
    if i == 1, C(i) = -q_in*dr/k;
    elseif i == n_mid, C(i) = T2;
    else, C(i) = S;
    end
end

%Adding in boundary conditions
A(1,1) = -1;
A(1,2) = 1;
A(n_mid,n_mid) = 1;

%Calculating Temperature
T = A\C;  %This looks backwards but is how matlab works

%Now find the temperature in the metal given the temperature and convective
%boundary condition on the outside







%Plotting the temperature over the radius
figure(1)
plot(r*10^3,T-273)
xlabel("Distance (mm)")
ylabel("Temperature (^oC)")

%%
%Plotting the analytical solution

%Temperature distribution in the mid portion
%Creating radius vector
ra = linspace(r1,r2,n_mid);
%Constants found by hand
C1 = -r1/k * (q_in - qdot/2*r1);
C2 = T2 + qdot/4/k * r2^2 - C1*log(r2);
%Temperature distribution from heat equation
Ta = -qdot/4/k*ra.^2 + C1.*log(ra) + C2;

%Temperature distribution in the metal sheet
%Creating radius vector
rm = linspace(r2,r3,n_met);
%Constants
C3 = -r3*h*(T3-Tinf);
C4 = T2 - C3/kg*log(r2);
%Temperature distribution
Tam = C3/kg.*log(rm)+C4;

figure(1)
hold on
plot(ra*10^3,Ta-273)
plot(rm*10^3,Tam-273)
legend('Numerical Solution','Analytical Solution')
title('Numerical and Analytical Temperature Distributions in the Mid Cylinder')
hold off
