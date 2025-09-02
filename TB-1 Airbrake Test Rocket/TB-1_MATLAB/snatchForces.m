%% Snatch Force Calculator
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
main = TB1.component(name = "Parachute");
drogue = TB1.component(name = "Streamer");

conv = 0.2248089431; % Newton to pound-force

%% Find Main Snatch Force
D = main.getDiameter();
Cd = main.getCD();
A = (pi/4)*D^2;

% trim data
data_range = timerange(eventfilter("LAUNCHROD"), eventfilter("MAIN"), "openleft");
mainData = simData(data_range, :);

% Assign data points
% Kyle method
T = (mainData.("Air temperature")); 
T = T(end);
P = mainData.("Air pressure")/1000; 
P = P(end);
v = mainData.("Total velocity");
v = v(end);
F_snatch1 = snatchCalc(Cd, A, v, P, T);

%% Find Drogue Snatch Force
L = drogue.getStripLength();
W = drogue.getStripWidth();
A = L*W;
Cd = drogue.getCD();
filter = eventfilter("DROGUE");
drogueData = simData(filter, ["Total velocity", "Air temperature", "Air pressure"]);
T = drogueData.("Air temperature");
P = drogueData.("Air pressure")/1000;
v = drogueData.("Total velocity");
F_snatch2 = snatchCalc(Cd, A, v, P, T);

%% Display output
disp("Main snatch: " + F_snatch1/1000 + " [kN], " + F_snatch1*conv + " [lbf]"...
    + newline + "Drogue snatch: " + F_snatch2/1000 + " [kN], "...
    + F_snatch2*conv + " [lbf]")

%% Calculator Function
% Arguments are drag coefficient, reference area, airspeed, pressure,
% and temperature
function F = snatchCalc(Cd, A, v, P, T)
    R = 0.287;
    rho = P/(R*T);
    F = 0.5*Cd*A*rho*v^2;
end