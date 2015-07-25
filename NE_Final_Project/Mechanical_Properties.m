% NE 472 Mechanical Properties of INCOLOY alloy 800H
% Trey Huffine
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            Must use MKS units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all, clear all, clc;

%% Variables

creep_percent = 1; % Total creep

P1 = 101325;
P2 = 2*P1;
P3 = 3*P1;

Diameter = 5; % diamter of core

time1 = 30*365*24 % time - convering years to hours
time2 = 40*365*24
time3 = 50*365*24

creep_rate1 = creep_percent/time1  % creep per hour
creep_rate2 = creep_percent/
creep_rate3 = creep_percent/time3

stress1 = 6.3*10^6;  % Stress in Pa from Special Metals
stress2 = 10.5*10^6;
stress3 = 4.95*10^6;

% at 1 ATM
thickness1 = P1*Diameter/(2*stress1) % in meters
thickness2 = P1*Diameter/(2*stress2)
thickness3 = P1*Diameter/(2*stress3)

thickin1 = thickness1*37.3701 % compare to thickness in inches
thickin2 = thickness2*37.3701
thickin3 = thickness3*37.3701

% at 2 ATM
thickness1 = P2*Diameter/(2*stress1) % in meters
thickness2 = P2*Diameter/(2*stress2)
thickness3 = P2*Diameter/(2*stress3)

thickin1 = thickness1*37.3701 % compare to thickness in inches
thickin2 = thickness2*37.3701
thickin3 = thickness3*37.3701

% at 3 ATM
thickness1 = P3*Diameter/(2*stress1) % in meters
thickness2 = P3*Diameter/(2*stress2)
thickness3 = P3*Diameter/(2*stress3)

thickin1 = thickness1*37.3701 % compare to thickness in inches
thickin2 = thickness2*37.3701
thickin3 = thickness3*37.3701

