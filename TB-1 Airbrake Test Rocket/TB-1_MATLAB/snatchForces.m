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

% run sim
mainData = openrocket.simulate(sim,outputs = ["Total velocity", "Air temperature", "Air pressure"]);

% trim data
data_range = timerange(eventfilter("LAUNCHROD"), eventfilter("APOGEE"), "openleft");
mainData = mainData(data_range, :);

% assign data points. now there is probably a better way to grab the end
% point but I forget its been a min
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

%% Calculator Function
% Arguments are drag coefficient, reference area, airspeed, pressure,
% and temperature
function F = snatchCalc(Cd, A, v, P, T)
    R = 0.287;
    rho = P/(R*T);
    F = 0.5*Cd*A*rho*v^2;
end