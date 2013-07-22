function [head, prof, profX] = readecmwf91_nearest_gasNcloudNEW(rtime, ...
                                                             lat, lon);

%%%%
% see /home/schou/Airs/mkcldfile_era_sergio.m
%%%%

% dumps in both gas profile info into prof and puts cloud info into 2 slabs,
% so that "ecmwf2cloud" and "old klayers" can be run!

% copied from /asl/matlab/gribtools/readecmwf91_nearest.m
% but now has the added functionality of including the cloud info, just like
% readecmwf91_grid_gasNcloud

%function [head, prof, profX] = readecmwf91_nearest_gasNcloud(mdate, ...
%                                                                 lat, lon);
% Routine to read in a 60 or 91 level ECMWF file and return a
% RTP-like structure of profiles that are the closest grid points
% to the specified (lat,lon) locations.
%
% Input:
%    rtime : Airs rtime (seconds since 1993)
%    lat : (1 x nprof) latitudes (degrees -90 to +90)
%    lon : (1 x nprof) longitude (degrees, either 0 to 360 or -180 to 180)
%
% Output:
%    head : (RTP "head" structure of header  info) 
%    prof : (RTP "prof" structure of gas + cloud info) 
%   profX : (RTP "prof" structure of cloud   info) 
%
% Note: uses external routines: p60_ecmwf.m, p91_ecmwf.m, readgrib_inv.m,
%    readgrib_rec.m, as well as the "wgrib" program.
%

% rtime = mattime2tai( datenum(2011,3,14,13,4,0) );
% readecmwf91_nearest_gasNcloud_slabprofNEW(rtime,0,0)

%
% Created: 19 Jan 2007  Sergio Machado
% Gas only 17 Mar 2006, Scott Hannon - re-write of old 60 level version
% Update:  23 October 2006, S.Hannon - changed from 60 level only to
%    60 or 91 level depending on file size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit this section as needed

%xstartup
addpath /home/sergio/MATLABCODE
addpath /asl/matlab/gribtools  % for readgrib_inv.m, readgrib_rec.m
addpath /asl/matlab/aslutil    % for mktemp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;

prof.rlat = lat;
prof.rlon = lon;
prof.rtime = rtime;
hattr = [];
pattr = set_attr([],'rtime','Seconds since 1993','profiles');




% Assign the output header structure
head.ptype = 0;
head.pfields = 1;
%head.pmin = min( prof.plevs(1,:) );
%head.pmax = max( prof.plevs(nlev,:) );
head.ngas = 2;
head.glist = [1; 3];
head.gunit = [21; 21];
head.nchan = 0;
head.mwnchan = 0;

fields = {'SP','SKT','10U','10V','TCC','CI','T','Q','O3','CC','CIWC','CLWC'};
[head, hattr, prof, pattr] = rtpadd_ecmwf_data(head, hattr, prof, pattr, fields);
keyboard

profX.plat     = prof.plat; 
profX.plon     = prof.plon; 
profX.stemp    = prof.stemp; 
profX.spres    = prof.spres; 
profX.plevs    = prof.plevs; 
%profX.landfrac = prof.landfrac;

%%% end of function %%%

tnow = toc; 
fprintf(1,' took %8.6f minutes to read in ECMWF file \n',tnow/60); 
 
disp('read the ECMWF file ... parsing cloud info ... '); 

tic;
nlev = nanmax(prof.nlevs);
profX = prof;
ecmwfcld2sartacld;        %%% set  the cloud profile info here 

%%%%
error('see /home/schou/Airs/mkcldfile_era_sergio.m')
%%%%

%%%%%%%%%%%%% this is to check versus klayers_trace %%%%%%%%%%%%%%%%%%%%%

iCheck = -1;
if iCheck > 0
  check_klayers_works
  end

tnow = toc; 
fprintf(1,'TOTAL : %8.6f more minutes to process slab info \n',tnow/60); 
