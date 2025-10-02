%% Monte Carlo Simulation
% Testbed TB-1
% Last updated 1 October 2025
%% Setup
clear; close all; clc;
% Run a test OR sim
filePath = "C:\ltestbed\TB-1 Airbrake Test Rocket\TB-1.ork"; 
if ~isfile(filePath)
    error("Error: not on path", filePath);
end
TB1 = openrocket(filePath);
% "15mph_URRG", "15mph_URRG_K455", and "15mph_URRG_L935" are valid currently
simName = "15mph_URRG";
sim = TB1.sims("15mph_URRG");
opts = sim.getOptions();
windBounds = [1.5 6.5]; % [min max] [m/s]
windRange = windBounds(2)-windBounds(1);
tempBounds = [260 300]; % [K]
tempRange = tempBounds(2)-tempBounds(1);
pressBounds = [90 100]*10^3; % [Pa]
pressRange = pressBounds(2) - pressBounds(1);
appReduction = 248;

%% Monte Carlo Loop
N = 200; % Number of samples
apogeeList = zeros([N,1]);
pressAppList = zeros([N,1]);
elapsed = tic;
for i = 1:N
    disp("Running simulation " + i + " of " + N)
    wind = windBounds(1) + rand()*windRange;
    temp = tempBounds(1) + rand()*tempRange;
    press = pressBounds(1) + rand()*pressRange;
    
    opts.setWindSpeedAverage(wind);
    opts.setLaunchTemperature(temp);
    opts.setLaunchPressure(press);
    TB1.simulate(sim);
    altData = openrocket.get_data(sim, [("Altitude"), ("Air pressure")]);
    apogeeList(i) = max(altData.("Altitude"));
    pressAppList(i) = pressalt("m", min(altData.("Air pressure")), "Pa")-pressalt("m", altData.("Air pressure")(1), "Pa");
end
%% Analysis
appErr = pressAppList - apogeeList; % Supposed measurement error
% Averages
avgAlt = mean(apogeeList); 
avgPressAlt = mean(pressAppList);
avgErr = mean(appErr);
% Standard deviations
sigAlt = std(apogeeList);
sigPressAlt = std(pressAppList);
sigErr = std(appErr);
% Conversions factors
C1 = 2.2369; % m/s to mph
tempsF = (tempBounds-273)*1.8 + 32;
pressBounds = pressBounds*10^-3;
%% Output
fprintf("\n%d Simulations run varying parameters in the listed ranges", N);
fprintf("\n   Wind Speed: %1.0f to %1.0f [m/s]; %2.0f to %2.0f [mph]",...
    windBounds(1), windBounds(2), windBounds(1)*C1, windBounds(2)*C1);
fprintf("\n   Launch Temperature: %2.1f to %2.1f [Celcius]; %2.1f to %2.1f [Fahrenheit]",...
    tempBounds(1)-273, tempBounds(2)-273, tempsF(1), tempsF(2));
fprintf("\n   Launch Pressure: %4.2f to %4.2f [kPa]",...
    pressBounds(1), pressBounds(2));
fprintf("\nISA used for atmosphere model");
fprintf("\n2 sigma bounds");

fprintf("\n\nApogee (geometric): %4.1f [m] +/- %3.1f [m]\n", avgAlt, 2*sigAlt);
fprintf("Apogee (indicated): %4.1f [m] +/- %3.1f [m]\n", avgPressAlt, 2*sigPressAlt);
fprintf("Apogee error: %2.1f [m] +/- %2.1f [m]\n", avgErr, 2*sigErr);

fprintf("\nAirbrake apogee reduction of %3.0f [m] corresponds to %2.1f sigma\n", appReduction, (appReduction/sigAlt));