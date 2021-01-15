function [mm] = mmwat(gas1,nlevs);

% ASSUMES 98 levels!!

% Avagadro's number (molecules per mole)
navagadro=6.02214199E+23;

% Approximate mass of water per molecule (AMU)
mass=18.015;
cfac=10*mass/navagadro;

[~,nprof]=size(gas1');

mm = cfac*nansum(gas1');

