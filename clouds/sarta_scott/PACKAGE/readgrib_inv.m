function [rec,param,level] = readgrib_inv(file);

% function [rec,param,level] = readgrib_inv(file);
%
% Read an inventory of a GRIB file.  This is basically
% just a MATLAB wrapper for the "wgrib" program.
%
% Input:
%    file : {string} name of an NCEP/ECMWF model grib file
%
% Output:
%    rec   : [1 x Nrec] {integer} record number
%    param : [1 x Nrec] {string cell} parameter name
%    level : [1 x Nrec] {string cell} type of level
%          
% Requires:  wgrib in your path.  wgrib is a C program run via the
%	     shell that pull fields out of the grib file in a format
%            that can be read by Matlab.
%

% Created: 15 Mar 2006, Scott Hannon
% Update: 16 May 2006, Scott Hannon - change temporary inventory file
%    name from "dump" to "inv_<file>".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = dir(file);
invfile = ['inv_' d.name];

% Run wgrib and create temporary inventory file
eval(['! wgrib -s ' file ' | cut -f1,4,5 -d":" > ' invfile]);

% Read temporary text file
[rec,param,level] = textread(invfile,'%n%s%s\n','delimiter',':');

% Remove temporary inventory file
% eval(['! rm ' invfile])
%error('pah');

%%% end of function %%%
