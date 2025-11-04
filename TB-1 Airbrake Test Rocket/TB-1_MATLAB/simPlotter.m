%% Simulation Result Plotter
% Created by Juno Afko
% Last updated November 4, 2025
clear;
close all;
clc;
%% Setup
filePath = "C:\ltestbed\TB-1 Airbrake Test Rocket\TB-1.ork"; 
if ~isfile(filePath)
    error("Error: not on path", filePath);
end
% Get alternate drag curve
dragFilePath = "C:\ltestbed\TB-1 Airbrake Test Rocket\Data\TB1_RasDrag.CSV";
rasDrag = import_rasaero_aerodata(dragFilePath);
rasDrag = rasDrag.align("mach");
rasDrag = table(rasDrag.mach, rasDrag.pick{"aoa", 0, "field", "CD"}, ...
    VariableNames = ["MACH", "DRAG"]);
rasDrag.MACH(1) = 0;
% Other parameters
D = 6.2;
L = 79.7;

%% Simulate
TB1 = openrocket(filePath);
sim = TB1.sims("15mph_URRG");
openrocket.simulate(sim, outputs = "ALL", drag = rasDrag);
simData = openrocket.get_data(sim);
time = round(seconds(simData.Time), 3);
% Process data
% Adding stability percent, finding min and max
simData.("Stability percent") = 100*simData.("Stability margin")*D/L;
stabs = [min(simData.("Stability percent")), max(simData.("Stability percent"))];
% Flutter FoS
fins = TB1.component(class = "FinSet");
flutterMargin = FOS_finflutter(simData, fins);

%% Plot performance data
figure(name = "Flight Performance");
plot_openrocket(simData, "Altitude", "Total velocity", ...
    start_ev = "LAUNCH", end_ev = "GROUND_HIT", labels = ["BURNOUT", "APOGEE", "MAIN"]);

figure(name = "Ascent stability");
plot_openrocket(simData, "Stability percent", ...
    start_ev = "LAUNCH", end_ev = "APOGEE", labels = ("LAUNCHROD"));

%% Text output
fprintf("\nStablity at launch: %2.2f percent\nMaximum stability: %2.2f percent", stabs(1), stabs(2))
fprintf("\nFin flutter FoS: %1.2f\n", flutterMargin);