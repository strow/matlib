function y = ppmv2gg(p,t,q,mass_g)

%% ppmv to grams/g  ie units 10 to units 21

%% modelled on toppmv.m file, in this directory

%% modelled on Scott Hannon's toppmv.f file, in the KLAYERS dist
%% uses the same GUNIT values as would be in header.gunit (of a .rtp file)
%% input   p = pressure        (mb)
%%         t = temperature     (kelvin)
%%         q = gas amount      (whatever units in ppmv)
%%    mass_g = molar gas mass  (grams)

%%  iGasUnit = 21 %%%(from .rtp format) convert to gg

mdair = 28.966;  % molecular mass of dry air

% PPMV = MRgg*((MDAIR)/MASSF)*1E+6
rjunk =1E+6 * mdair/mass_g;
rjunk = 1/rjunk;
y = q * rjunk;

