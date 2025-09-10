%% Airflow Simulation for Puck Avionics Bay
% Created by Juno Afko
% Last updates September 10, 2025
%% Problem Statement and Assumptions 
% The experimental avionics bay design, known as the puck avionics bay,
% burger avbay, or sandwich avbay, raises concerns over whether venting is
% adequate for the flight computers. This script will assess the concern by
% simulating airflow throughout the puck avbay and comparing the pressure
% in the flight computer puck to a reference avionics bay simulated using
% the same approach. 
% 
% The puck avionics bay is modelled here as three discrete volumes:
% V1: Flight computer puck
% V2: Tracker puck
% V3: Intermediary volume (this arises because of the current puck design)
%
% The flight computer puck is connected to the intermediary volume by what
% is effectively a short tube; flow losses through this tube will be
% modeled as the sum of the inlet loss, exit loss, and viscous pressure
% loss. Losses at the orifice between V2 and V3 are always modeled as an
% exit; the dynamic pressure is discarded.
%
% Two cases are considered; in Case #1, the vent holes are positioned at
% the tracker disk. In Case #2, the vent holes are positioned at the flight
% computer disk. The difference can be seen via images in the figures
% folder within this repository.
%
% The puck avionics bay is compared to a conventional avbay, in which a
% single volume vents air to outside. Pressure losses through this hole are
% assumed to be equal to the dynamic pressure.
%% Setup
clc
clear;
close all;
% Run a test OR sim
filePath = "C:\ltestbed\TB-1 Airbrake Test Rocket\TB-1.ork"; 
if ~isfile(filePath)
    error("Error: not on path", filePath);
end
TB1 = openrocket(filePath);
sim = TB1.sims("15mph_URRG");
openrocket.simulate(sim);
% Get data and trim
simData = openrocket.get_data(sim);
data_range = timerange(eventfilter("LAUNCH"), eventfilter("MAIN"), "openleft");
data = simData(data_range, [("Air pressure"), ("Air temperature")]);
time = round(seconds(data.Time), 3);
Pdata = griddedInterpolant(time, data.("Air pressure")*10^-3); % interpolants for sim loop
Tdata = griddedInterpolant(time, data.("Air temperature"));
% Define airflow sim parameters
% initialize ambient conditions
air.R = 0.278;
air.P = data.("Air pressure")(1)*10^-3;
air.T = data.("Air temperature")(1);
air.rho = air.P/(air.R*air.T);
air.Cp = 1.005;
air.h = air.T*air.Cp;
% Control volumes
V1 = initControlVolume(1.004*10^-3, air);
V2 = initControlVolume(1.016*10^-3, air);
V3 = initControlVolume(8.554*10^-4, air);
conventional = initControlVolume(3.404*10^-3, air);
% Other parameters
A1 = 8.0623*10^-4;
A2 = 8.6*10^-4;
Avent = 90.008*10^-6;
L = 0.05;
f_L = 0.05;
D_L = sqrt((4*A1/pi)); % effective pipe diameter
k1 = 1; % pipe exit coefficient
k2 = 0.5; % pipe entrance coefficient
C_L = 1+k1+k2+f_L*(L/D_L);
% set timer parameters
t_end = time(end);
dt = 5*10^-3;

%% Simulate Puck Avbay; Vents in tracker puck
% setup timer and data collection
t = 0;
simTime = [];
pressures = [];
Verr = [];
% iterate
while (t <= t_end)
    % Update ambient conditions
    air.P = Pdata(t);
    air.T = Tdata(t);
    air.rho = air.P/(air.R*air.T);
    % Calculate internal mass and enthalpy flows
    [mdot, Hdot] = flowCalc(V1, V2, V3, A1, A2, C_L);
    % Vent hole flows
    [mdot(3), Hdot(3)] = ventCalc(V2, air, Avent);
    % Iterate masses and enthalpies
    % dM and dH change depending on vent hole locations
    dM = [-mdot(1)*dt, (mdot(1)-mdot(2))*dt, (mdot(2)-mdot(3))*dt];
    dH = [-Hdot(1)*dt, (Hdot(1)-Hdot(2))*dt, (Hdot(2)-Hdot(3))*dt];
    V1 = updateVolume(V1, dM(1), dH(1), air);
    V3 = updateVolume(V3, dM(2), dH(2), air);
    V2 = updateVolume(V2, dM(3), dH(3), air);
    pressures = [pressures; V1.P, V2.P, V3.P];
    Verr = [Verr, 100*(V1.P-air.P)/air.P];
    % Iterate timer
    simTime = [simTime; t];
    disp(t);
    t = t + dt;
end
% plot results
figure(1);
plotPuck(time, simTime, pressures, Verr, data, "Puck Avbay, Case 1");

%% Simulate Puck Avbay; Vents in flight computer puck
% Re-init volumes
air.P = data.("Air pressure")(1)*10^-3;
air.T = data.("Air temperature")(1);
air.rho = air.P/(air.R*air.T);
air.h = air.T*air.Cp;
V1 = initControlVolume(1.004*10^-3, air);
V2 = initControlVolume(1.016*10^-3, air);
V3 = initControlVolume(8.554*10^-4, air);
% setup timer and data collection
t = 0;
simTime = [];
pressures = [];
Verr = [];
% iterate
while (t <= t_end)
    % Update ambient conditions
    air.P = Pdata(t);
    air.T = Tdata(t);
    air.rho = air.P/(air.R*air.T);
    % Calculate internal mass and enthalpy flows
    [mdot, Hdot] = flowCalc(V1, V2, V3, A1, A2, C_L);
    % Vent hole flows
    [mdot(3), Hdot(3)] = ventCalc(V1, air, Avent);
    % Iterate masses and enthalpies
    % dM and dH change depending on vent hole locations
    dM = [-(mdot(1)+mdot(3))*dt, (mdot(1)-mdot(2))*dt, mdot(2)*dt];
    dH = [-(Hdot(1)+Hdot(3))*dt, (Hdot(1)-Hdot(2))*dt, Hdot(2)*dt];
    V1 = updateVolume(V1, dM(1), dH(1), air);
    V3 = updateVolume(V3, dM(2), dH(2), air);
    V2 = updateVolume(V2, dM(3), dH(3), air);
    pressures = [pressures; V1.P, V2.P, V3.P];
    Verr = [Verr, 100*(V1.P-air.P)/air.P];
    % Iterate timer
    simTime = [simTime; t];
    disp(t);
    t = t + dt;
end
% plot results
figure(2);
plotPuck(time, simTime, pressures, Verr, data, "Puck Avbay, Case 2");

%% Simulate Conventional Avbay
% setup timer and data collection
t = 0;
simTime = [];
convPressure = [];
Verr = [];
% iterate
while (t <= t_end)
    % Update ambient conditions
    air.P = Pdata(t);
    air.T = Tdata(t);
    air.rho = air.P/(air.R*air.T);
    % Vent hole flows
    [mdot, Hdot] = ventCalc(conventional, air, Avent);
    % Iterate mass and enthalpy
    dM = -mdot*dt;
    dH = -Hdot*dt;
    conventional = updateVolume(conventional, dM, dH, air);
    convPressure = [convPressure; conventional.P];
    Verr = [Verr, 100*(conventional.P-air.P)/air.P];
    % Iterate timer
    simTime = [simTime; t];
    disp(t);
    t = t + dt;
end
% plot results
figure(3);
plotConv(time, simTime, convPressure, Verr, data, "Conventional Avbay");

%% Functions
function controlVolume = initControlVolume(volume, air)
    controlVolume.V = volume; % it literally says what it is im not gonna tell u
    controlVolume.P = air.P; % air pressure
    controlVolume.T = air.T; % air temperature
    controlVolume.rho = air.rho; % air density
    controlVolume.M = volume*air.rho; % total air mass
    controlVolume.H = air.T*air.Cp*controlVolume.M; % total enthalpy
    controlVolume.h = controlVolume.H/controlVolume.M; % specific enthalpy
end
function [mdot, Hdot] = flowCalc(V1, V2, V3, A1, A2, C_L)
    rho1 = min(V1.rho, V3.rho);
    rho2 = min(V2.rho, V3.rho);
    v1 = sign(V1.P - V3.P)*sqrt(2*abs(V1.P-V3.P)/(rho1*C_L));
    v2 = sign(V3.P - V2.P)*sqrt(2*abs(V3.P-V2.P)/(rho2*C_L));
    mdot(1) = v1*A1*rho1;
    mdot(2) = v2*A2*rho2;
    if v1 >= 0
        Hdot(1) = mdot(1)*V1.h;
    else
        Hdot(1) = mdot(1)*V3.h;
    end
    if v2 >= 0
        Hdot(2) = mdot(2)*V3.h;
    else
        Hdot(2) = mdot(2)*V2.h;
    end
end
function [mdot, Hdot] = ventCalc(V1, air, Avent)
    rho = min(air.rho, V1.rho);
    v = sign(V1.P-air.P)*sqrt(2*(V1.P - air.P)/rho);
    mdot = v*rho*Avent;
    if v >= 0
        Hdot = mdot*V1.h;
    else
        Hdot = mdot*air.h;
    end
end

function Vol = updateVolume(Vol, dM, dH, air)
    Vol.M = Vol.M + dM;
    Vol.H = Vol.H + dH;
    Vol.h = Vol.H/Vol.M;
    Vol.T = Vol.h/air.Cp;
    Vol.rho = Vol.M/Vol.V;
    Vol.P = Vol.T*air.R*Vol.rho;
end
function plotPuck(time, simTime, pressures, Verr, data, titleStr)
    subplot(2, 1, 1);
    plot(time, data.("Air pressure")*10^-3, "b");
    hold on;
    plot(simTime, pressures(:,1), "r");
    plot(simTime, pressures(:,2), "m");
    plot(simTime, pressures(:,3), "k");
    legend("Ambient pressure", "CV 1", "CV 2", "CV 3");
    hold off;
    xlim([0, simTime(end)]);
    ylabel("Pressure [kPa]");
    title(titleStr);
    subplot(2, 1, 2)
    plot(simTime, Verr)
    xlim([0, simTime(end)]);
    xlabel("Time [s]");
    ylabel("P1 error [%]");
end
function plotConv(time, simTime, pressures, Verr, data, titleStr)
    subplot(2, 1, 1);
    plot(time, data.("Air pressure")*10^-3, "b");
    hold on;
    plot(simTime, pressures, "r");
    legend("Ambient pressure", "Avbay pressure");
    hold off;
    xlim([0, simTime(end)]);
    ylabel("Pressure [kPa]");
    title(titleStr);
    subplot(2, 1, 2)
    plot(simTime, Verr)
    xlim([0, simTime(end)]);
    xlabel("Time [s]");
    ylabel("P1 error [%]");
end