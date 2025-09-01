%% Simulation Result Plotter
% Created by Juno Afko
% Last updated September 1, 2025
clear;
close all;
%% Setup
filePath = "C:\ltestbed\TB-1 Airbrake Test Rocket\TB-1.ork"; 
if ~isfile(filePath)
    error("Error: not on path", filePath);
end

TB1 = openrocket(filePath);
sim = TB1.sims("15mph_URRG");
openrocket.simulate(sim);
simData = openrocket.get_data(sim);
time = round(seconds(simData.Time), 3);

%% Plot performance data
figure;
plot_openrocket(simData, "Altitude", "Total velocity", ...
    start_ev = "LAUNCH", end_ev = "GROUND_HIT", labels = ["BURNOUT", "APOGEE", "MAIN"]);

figure;
plot_openrocket(simData, "Stability margin", ...
    start_ev = "LAUNCH", end_ev = "APOGEE", labels = ["LAUNCHROD"]);