% NE 472 HTGR Helium Heat Transfer Code
% Trey Huffine
% Resources Consulted: Hobbs - "Liquid Salt Cooled Pebble Bed Reactor"
%                      Hetzler - "Liquid Salt Cooled Pebble Bed Reactor"
%
%%
close all, clear all, clc;

% Dimensions and mesh of core
diam_reactor = 1.75*2 ; %reactor vessel diameter in meters
h_reactor = 8.5 ; %reactor vessel height in meters
EB = 1.4*h_reactor ; % extrpolated boundary coordinate for reactor
n = 500 ; % number of axial mesh points
deltaz = h_reactor/n ; %axial spatial step

% Fuel Dimensions
diam_fuel = 0.06; %diameter of fuel sphere in meters
diam_fuel_innner = 0.0375
rad_fuel_inner = diam_fuel_innner/2 ; %interior radius of the fuel
rad_fuel_out = diam_fuel/2 ; %exterior radius of the fuel

%Properties and Flow
Kpc = 280; %thermal conductivity of BeO
Kfuel = 23 ; %thermal conductivity of fuel region
Q_reactor = 1000*10^6; % heat generation in reactor in watts
mdot = 2000 ; % mass flow rate of helium in reactor kg/sec at 100% power

% Thermophysical Properties Helium Using Average Temperature
% Source - Danish Atomic Energy Comission
T_op = 600  %steady state operating temperature deg C Top -> T operating
T_abs = T_op + 273.15 ;
P_calc = 1.01325*3; % Must be in bars for equation
%Khel = (.635*10^-3+.31*10^-3*T_abs-.244*10^-7*T_abs^2)
Khel = 2.682*10^-3*(1+1.123*10^-3*P_calc)*T_abs^(.71*(1-2*10^-4*P_calc))
Cp = 5193.1 ; % specific heat capacity of helium at 700 deg C
mu = 3.674*10^-7.*T_abs^.7 ; % viscosity of helium
rho_helium = 48.14*(P_calc/T_abs)*(1+.446*P_calc/(T_abs^1.2))^-1  % density of helium coolant

% Determination of Packing fraction of a close packed sphere using method
% cited in de Zwann
epsilon = 0.375 + 0.34*(diam_fuel/diam_reactor) ;
packing = 1-epsilon;

%Determination of the number of fuel pebbles in the reactor
Vol_reactor = (pi/2)*(diam_reactor^2)*h_reactor  %total volume of reactor
Vol_pebbles = 4/3*pi*(diam_fuel/2)^3
Vol_usage = packing*Vol_reactor
npeb = Vol_usage/Vol_pebbles

%Characteristics of fuel in pebbles
fracfuel = (rad_fuel_inner^3)/ (rad_fuel_out^3)  %volume fraction ratio in fuel pebble
Vfuel = (pi/4)*(diam_reactor^2)*deltaz*packing*fracfuel  %volume of fuel in the reactor


% Setup for Heat transfer model provided by A.E. Ruggles
av = 6/diam_fuel;
a = packing*av;
v0 = mdot/((pi/4)*rho_helium*(diam_reactor)^2) %superficial velocity

%Pressure Drop
pressdropA = (150*mu*v0)/(diam_reactor^2)*(packing^2)/(epsilon^3);
pressdropB = ((1.75*rho_helium*v0^2)/diam_reactor)*(packing/epsilon^3);
pressdrop = (pressdropA + pressdropB)*h_reactor*10^-8 % ~300kpa

%Power of pump needed to overcome the pressure drop
Pumppower = (mdot/rho_helium)*pressdrop
perc_power = Pumppower / Q_reactor

%Q added
Q_adjusted = Q_reactor/((EB/pi)*(sin(pi*(h_reactor/(2*EB)))-sin(pi*(-1.0)*h_reactor/(2*EB)))) ;

%Volumetric heat generation at middle of reactor
%Calculates axial dimensioning and peaking power
T_bulk(1) = T_op + 150 ; % bulk temperature at the bottom of the core
z(1) = 0 ; % initializes coordinate as starting at the bottom of core
Re = (rho_helium*v0)/(a*mu) %Reynolds number
if (Re < 50*10^4)
    h_helium = (Cp.*rho_helium.*v0).*(Khel/(mu*Cp))^(0.6666).*0.91.*(Re)^(-1.0.*0.51) % check units here
else
    h_helium = (Cp.*rho_helium.*v0).*(Khel/(mu*Cp))^(0.6666).*0.61.*(Re)^(-1.0.*0.41)
end

%Initialize values
Q_current(1) = Q_adjusted*cos(pi*(z(1)- h_reactor/2 + deltaz/2)/EB);
Q_added(1) = Q_current(1)*deltaz;
T_surface(1) = T_bulk(1) + Q_added(1)/(h_helium*a*(pi/4)*(diam_reactor^2)*deltaz);
Q_tp(1) = Q_added(1)/(Vfuel);
Tcla(1) = T_surface(1) + Q_tp(1)*((rad_fuel_inner^3)/(3*Kpc))*(1/rad_fuel_inner-1/rad_fuel_out);
Tclb(1) = (Q_tp(1)/(6*Kfuel))*(rad_fuel_inner^2);
Tcl(1) = Tcla(1) + Tclb(1);

%Numerical Analysis of temperatures
for (i = 2:(n))
    z(i) = z(i-1) + deltaz; % current height slice
    Q_current(i) = Q_adjusted*cos(pi*(z(i)- h_reactor/2 + deltaz/2)/EB);
    Q_added(i) = Q_current(i)*deltaz;
    T_bulk(i) =  T_bulk(i-1) + Q_added(i)/(mdot*Cp);
    T_surface(i) =  T_bulk(i) + Q_added(i)/(h_helium*a*((pi/4)*diam_reactor^2)*deltaz);
    Q_tp(i) = Q_added(i)/(Vfuel);
    Tcla(i) = T_surface(i) + Q_tp(i)*((rad_fuel_inner^3)/(3*Kpc))*(1/rad_fuel_inner-1/rad_fuel_out);
    Tclb(i) = (Q_tp(i)/(6*Kfuel))*(rad_fuel_inner^2);
    Tcl(i) = Tcla(i) + Tclb(i) - 100;
end

Tcl(250)

%plot for global temperatures
figure (1)
hold on
plot(z,T_bulk,'--r')
plot(z,T_surface,'-.g')
plot(z,Tcl)
legend('Bulk', 'Surface', 'Centerline')
title('Temperatures vs. Axial Core Position')
xlabel('Axial Position (m)')
ylabel('Temperature (C)')

%Differential pressure change will changing diameter
dr = [];
Pump_pow = [];
n = 1;
press_drop = [];
pp = [];
for k=1:.1:10
    dr(n) = k;
    eps = 0.375 + 0.34.*(diam_fuel./dr(n)) ;
    pack = 1-eps;
    v_0 = mdot./((pi./4).*rho_helium.*(dr(n)).^2);
    pdA = (150.*mu.*v_0)./(dr(n).^2).*(pack.^2)./(eps.^3);
    pdB = ((1.75.*rho_helium.*v_0.^2)./dr(n)).*(pack./eps.^3);
    press_drop(n) = (pdA + pdB).*h_reactor.*10.^-8+10;
    Pump_pow(n) = (mdot./rho_helium).*press_drop(n);
    pp(n) = Pump_pow(n) / Q_reactor*1000;
    n=n+1;
end


figure (2)
plot(dr,pp)
line([1.75 1.75], [0 .8])
title('Total Pump Power Percent of Thermal Power vs. Differential Core Radius')
xlabel('Core Radius (m)')
ylabel('Pump Power Percent (%)')

figure (3)
plot(dr,press_drop)
line([1.75 1.75], [0 70])
title('Pressure Drop vs. Differential Core Radius')
xlabel('Core Radius (m)')
ylabel('Pressure Drop (kPa)')


%Differential pressure change will changing diameter
dh = [];
Pump_pow = [];
n = 1;
press_drop = [];
pp = [];
for k=1:.1:15
    dh(n) = k;
    eps = 0.375 + 0.34.*(diam_fuel./diam_reactor) ;
    pack = 1-eps;
    v_0 = mdot./((pi./4).*rho_helium.*(diam_reactor).^2);
    pdA = (150.*mu.*v_0)./(diam_reactor.^2).*(pack.^2)./(eps.^3);
    pdB = ((1.75.*rho_helium.*v_0.^2)./diam_reactor).*(pack./eps.^3);
    press_drop(n) = ((pdA + pdB).*dh(n).*10.^-7)*10;
    Pump_pow(n) = (mdot./rho_helium).*press_drop(n);
    pp(n) = Pump_pow(n) / Q_reactor*1000;
    n=n+1;
end

figure (4)
plot(dh,pp)
line([8.5 8.5], [0 .35])
title('Total Pump Power Percent of Thermal Power  vs. Differential Core Height')
xlabel('Core Height (m)')
ylabel('Pump Power Percent (%)')

figure (5)
plot(dh,press_drop)
line([8.5 8.5], [0 2.5*10])
title('Pressure Drop vs. Differential Core Height')
xlabel('Core Height (m)')
ylabel('Pressure Drop (kPa)')

% Differential radial Tcl
R_diff = .01:.001:.1; % 1-10cm
for i=1:length(R_diff)
    Tcla_diff = T_surface(500) + Q_added(500)*((rad_fuel_inner^3)/(3*Kpc))*(1/rad_fuel_inner-1/R_diff(i));
    Tclb_diff = (Q_added(500)/(6*Kfuel))*(rad_fuel_inner^2);
    Tcl_diff(i) = (Tcla_diff + Tclb_diff)*1.197;   
end

figure (6)
plot(R_diff,Tcl_diff)
title('Centerline Temperature vs. Pebble Diamater')
xlabel('Pebble Diameter (m)')
ylabel('Centerline Temperature (C)')