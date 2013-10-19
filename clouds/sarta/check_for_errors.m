function [prof,iNotOK,ibadlist] = check_for_errors(pin,cfracSet,iFixLoop)

iNotOK = 0;
prof = pin;
ibadlist = [];

if cfracSet > 0
  ibad = find(prof.cngwat > 0 & prof.cfrac ~= cfracSet);
  ibadlist = [ibadlist ibad];
  fprintf(1,'  stage0 : fixing  %6i        ibad cfracs \n',length(ibad));
  prof.cfrac(ibad) = cfracSet;

  ibad = find(prof.cngwat2 > 0 & prof.cfrac2 ~= cfracSet);
  ibadlist = [ibadlist ibad];
  fprintf(1,'  stage0 : fixing  %6i        ibad cfracs \n',length(ibad));
  prof.cfrac2(ibad) = cfracSet;
end

%% cfrac checks
ibad = find(prof.cfrac12 > prof.cfrac | prof.cfrac12 > prof.cfrac2);
ibadlist = [ibadlist ibad];
iNotOK = iNotOK + length(ibad);
prof.cfrac12(ibad) = min(prof.cfrac(ibad), prof.cfrac2(ibad));
fprintf(1,'  stage1 : fixing  %6i        ibad cfracs \n',length(ibad));

%% need cprtop > 0
ibad = find(prof.cprtop <= 0 & prof.cprbot > 0);
ibadlist = [ibadlist ibad];
iNotOK = iNotOK + length(ibad);
prof.cprtop(ibad) = prof.cprtop(ibad) + 50;
prof.cprbot(ibad) = prof.cprbot(ibad) + 50;
ibad2 = find(prof.cprtop2 <= 0 & prof.cprbot2 > 0);
ibadlist = [ibadlist ibad2];
iNotOK = iNotOK + length(ibad2);
prof.cprtop2(ibad2) = prof.cprtop2(ibad2) + 50;
prof.cprbot2(ibad2) = prof.cprbot2(ibad2) + 50;
fprintf(1,'  stage2 : fixing  %6i %6i ibad cprtop < 0 \n',length(ibad),length(ibad2));

%% need cprbot < cprtop2 
ibad = find(prof.cprbot >= prof.cprtop2 & prof.cprbot > 0 & prof.cprtop2 > 0);
ibadlist = [ibadlist ibad];
iNotOK = iNotOK + length(ibad);
junk = prof.cprbot(ibad);
prof.cprbot(ibad) = prof.cprtop2(ibad)-10;
fprintf(1,'  stage3 : fixing  %6i        ibad cprbot2/cprtop1 \n',length(ibad));

% need cprtop < spres
ibad = find(prof.cprtop >= prof.spres & prof.cprtop > 0);
ibadlist = [ibadlist ibad];
prof.cprtop(ibad) = prof.spres(ibad)-50;
ibad2 = find(prof.cprtop2 >= prof.spres & prof.cprtop2 > 0);
ibadlist = [ibadlist ibad2];
prof.cprtop2(ibad2) = prof.spres(ibad2)-50;
fprintf(1,'  stage4 : fixing  %6i %6i ibad cprtop/spres \n',length(ibad),length(ibad2));
iNotOK = iNotOK + length(ibad) + length(ibad2);

%% need cprbot < spres
ibad = find(prof.cprbot >= prof.spres & prof.cprbot > 0);
ibadlist = [ibadlist ibad];
prof.cprbot(ibad) = prof.spres(ibad)-25;
ibad2 = find(prof.cprbot2 >= prof.spres & prof.cprbot2 > 0);
ibadlist = [ibadlist ibad2];
prof.cprbot2(ibad2) = prof.spres(ibad2)-25;
fprintf(1,'  stage5 : fixing  %6i %6i ibad cprbot/spres \n',length(ibad),length(ibad2));
iNotOK = iNotOK + length(ibad) + length(ibad2);

%% need cprtop < cprbot
ibad = find(prof.cprtop >= prof.cprbot & prof.cprtop > 0 & prof.cprbot > 0);
ibadlist = [ibadlist ibad];
prof.cprtop(ibad) = prof.cprbot(ibad)-25;
ibad2 = find(prof.cprtop2 >= prof.cprbot2 & prof.cprtop2 > 0 & prof.cprbot2 > 0);
ibadlist = [ibadlist ibad2];
prof.cprtop2(ibad2) = prof.cprbot2(ibad2)-25;
fprintf(1,'  stage6 : fixing  %6i %6i ibad cprtop/cprbot \n',length(ibad),length(ibad2));
iNotOK = iNotOK + length(ibad) + length(ibad2);
%[prof.cprtop(2832) prof.cprbot(2832) prof.cprtop2(2832) prof.cprbot2(2832)]

ibadlist = unique(ibadlist);

%%%%%%%%%%%%%%%%%%%%%%%%%
if iFixLoop >= 10 & iNotOK > 0
  %% if things aren't ok by now, do drastic measures!!!!
  fprintf(1,'iFixLoop >= 10, iNotOK == %4i getting desperate!!! \n',iNotOK')

  ibadX = find(prof.cprtop < 0);
  if length(ibadX) > 0
    fprintf(1,'desperate mode : found %4i prof.cprtop < 0 \n',length(ibadX));
    prof.cprtop(ibadX) = 10;
    prof.cprbot(ibadX) = 20;
    prof.cngwat(ibadX) = 0.0;
    prof.cfrac(ibadX) = 0.0;
    prof.cfrac12(ibadX) = 0.0;
 end

  ibadX = find(prof.cprtop > prof.spres);
  if length(ibadX) > 0
    fprintf(1,'desperate mode : found %4i prof.cprtop > prof.spres \n',length(ibadX));
    prof.cprtop(ibadX) = prof.spres(ibadX)-50;
    prof.cprbot(ibadX) = prof.spres(ibadX)-10;
    prof.cngwat(ibadX) = 0.0;
    prof.cfrac(ibadX) = 0.0;
    prof.cfrac12(ibadX) = 0.0;
 end

  ibadX = find(prof.cprtop2 < 0);
  if length(ibadX) > 0
    fprintf(1,'desperate mode : found %4i prof.cprtop2 < 0 \n',length(ibadX));
    prof.cprtop2(ibadX) = 10;
    prof.cprbot2(ibadX) = 20;
    prof.cngwat2(ibadX) = 0.0;
    prof.cfrac2(ibadX) = 0.0;
    prof.cfrac12(ibadX) = 0.0;
 end

  ibadX = find(prof.cprtop2 > prof.spres);
  if length(ibadX) > 0
    fprintf(1,'desperate mode : found %4i prof.cprtop2 > prof.spres \n',length(ibadX));
    prof.cprtop2(ibadX) = prof.spres(ibadX)-50;
    prof.cprbot2(ibadX) = prof.spres(ibadX)-10;
    prof.cngwat2(ibadX) = 0.0;
    prof.cfrac(ibadX) = 0.0;
    prof.cfrac12(ibadX) = 0.0;
 end

  ibadX = find(prof.cprbot >= prof.cprtop2 & prof.cprbot > 0 & prof.cprtop2 > 0);
  if length(ibadX) > 0
    fprintf(1,'desperate mode : found %4i prof.cprbot > prof.cprtop2, turn off lower cloud \n',length(ibadX));
    prof.cprtop2(ibadX) = -9999;
    prof.cprbot2(ibadX) = -9999;
    prof.cngwat2(ibadX) = 0.0;
    prof.cfrac2(ibadX) = 0.0;
    prof.cfrac12(ibadX) = 0.0;
 end

end