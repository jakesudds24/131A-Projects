clc;
clear all;
close all;

%% Properties of Flue Gas (Air)
% at temp of 425K

rhoG = .823; % kg/m^3
mewG = 240.4*10^-7; % Ns/m^2
%veeyG = 26.41*10^-6; % m^2/s
kG = 35.6*10^-3; % W/mK
%alphaG = 38.3*10^-6; % m^2/s
prG = .689;
cpG = 1018; % J/kgK

%% Properties of Water
% at mean temperature of 32.5 C ~305K

rhoW = 995; % kg/m^3
mewW = 769*10^-6; % Ns/m^2
%veeyW = 7.532*10^-7; % m^2/s
kW = .620; % W/mK
%alphaW = 1.46*10^-7; % m^2/s
prW = 5.20;
cpW = 4178; % J/kgK
hfgw = 2426*1000; % J/kg

%% Set Parameters

Tci = 15+273; % k
Tco = 50+273; % k
Thi = 200+273; % k
Tho = 100+273; % k

mwtotal = 5000; % kg/s

%% Inputs
N = 2000; % number of tubes

Di = .025; % m
Do = .05; % m

%% Calcs

mwper = mwtotal / N; % kg/s mdot of water per tube

% find mass flow rate of flue gas
q = mwtotal * cpW * (Tco-Tci); % total heat transfer
mgtotal = q / (cpG * (Thi-Tho)); % total mass flow rate of flue gas
mgper = mgtotal / N; % mass flow rate of flue gas per tube
qper = q / N; % heat transfer rate per tube

% find log mean temperature difference
dT1 = Thi - Tci;
dT2 = Tho - Tco;
LMTD = (dT2-dT1)/log(dT2/dT1);

% for inner tube of water
ReW = (4*mwper)/(mewW*pi*Di);
NuW = .023*(ReW^.8)*(prW^.4); % n = .4 for cooling
hW = (NuW*kW)/Di;

% for outer tube of flue gas (air)
Dh = Do-Di;
ReG = (4*mgper)/(mewG*pi*Dh);
NuG = .023*(ReG^.8)*(prG^.3); % n = .3 for cooling
hG = (NuG*kG)/Dh;

% combine convection coefficients for overall convection coefficient
U = ((hW^-1) + (hG^-1))^-1;

% find length of tubes
As = qper / (U*LMTD); % surface area of tube

L = As / (pi*Di); % length per tube

% pressure drop of water
f = 0;
NuIter = 0;
Counter = 0;
while abs(NuIter-NuW) > 1
    f = f+.001;
    NuIter = ((f/8)*(ReW-1000)*prG) / (1+12.7*((f/8)^.5)*((prW^(2/3))-1));
    Counter = Counter + 1;
end

Ac = .25*pi*(Di^2);
Vw = mwper/(rhoW*Ac);

deltaP = f*(L/Di)*.5*rhoW*(Vw^2);