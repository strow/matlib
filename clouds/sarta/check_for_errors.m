function prof = check_for_errors(pin)

prof = pin;

% Error checks:

ii = prof.cfrac12 > prof.cfrac | prof.cfrac12 > prof.cfrac2;
prof.cfrac12(ii) = min(prof.cfrac(ii), prof.cfrac2(ii));

prof.cprbot(prof.cprbot > prof.spres) = prof.spres(prof.cprbot > prof.spres);
prof.cprbot2(prof.cprbot2 > prof.spres) = prof.spres(prof.cprbot2 > prof.spres);
% Error: CPRTO2 > CPRBO2
prof.cprtop2(prof.cprtop2 > prof.cprbot2) = prof.cprbot2(prof.cprtop2 > prof.cprbot2);
% CPRTO1 outside allowed PLEVS1 to  SPRES range
prof.cprtop(prof.cprtop > prof.spres) = prof.spres(prof.cprtop > prof.spres);

% extra check
prof.cprtop(prof.cprtop > prof.cprbot) = prof.cprbot(prof.cprtop > prof.cprbot);

