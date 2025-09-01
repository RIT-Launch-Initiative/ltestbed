%% Fin Optimizer for TB-1
% Created by Juno Afko
% Derived from finoptimization by Ares Bustinza-Nguyen
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
fins = TB1.component(class = "FinSet");
if ~isscalar(fins)
    error("Error: multiple fin sets found");
end