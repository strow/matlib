function [emis efreq nemis] = interp_emis(land_efreq,land_emis,sea_efreq,sea_emis,landfrac,interp_efreq)
% [emis efreq nemis] = interp_emis(land_efreq,land_emis,sea_efreq,sea_emis,landfrac,interp_efreq)
%
% Generic function to interpolate land and sea emissivities.
%
% INPUTS
%   land_efreq   - land emissivity frequencies
%   land_emis    - land emissivities
%   sea_efreq    - sea emissivity frequencies
%   sea_emis     - sea emissivities
%   landfrac     - land fraction to linearly interpolate on
%   interp_efreq - final frequency space of the interpolated emissivities (optional)
%

% Written by Paul Schou  (27 July 2011)

if size(land_emis,2) ~= size(sea_emis,2)
  error('Land and sea emissivities must have the same number of profiles')
end
if size(land_emis,1) ~= size(land_efreq,1)
  error('Land_emis and land_efreq must be the same length')
end
land_efreq = land_efreq(:,1);
if size(sea_emis,1) ~= size(sea_efreq,1)
  error('Sea_emis and sea_efreq must be the same length')
end
sea_efreq = sea_efreq(:,1);

% interpolation frequency set
if nargin < 6
  interp_efreq = sort([land_efreq;sea_efreq]);
  interp_efreq = interp_efreq(interp_efreq >= max(min(land_efreq),min(sea_efreq)) & ...
                              interp_efreq <= min(max(land_efreq),max(sea_efreq)));
end
if size(interp_efreq,1) < 2
  error('Interp_freq not a valid size or dimension')
end

npro = size(land_emis,2);
land_nemis = size(land_emis,1);
sea_nemis = size(sea_emis,1);
interp_nemis = size(interp_efreq,1);

% predeclare the size of arrays
emis = zeros([max(land_nemis,sea_nemis), npro],'single');
efreq = ones([max(land_nemis,sea_nemis), npro],'single') * -9999;
nemis = zeros([1, npro],'uint8');

isea = landfrac < 1;
isea(max(sea_emis) < -9) = 0;
iland = landfrac > 0;
iland(max(land_emis) < -9) = 0;

% make sure we return emis data (sea emis), even if it is over land
isea(~isea & ~iland) = 1;

% replace the clean stuff, here we put all the sea fovs into the sea emis
emis(1:sea_nemis,isea & ~iland) = sea_emis(:,isea & ~iland);
efreq(1:sea_nemis,isea & ~iland) = repmat(sea_efreq(:),[1 sum(isea & ~iland)]);
nemis(1,isea & ~iland) = sea_nemis;
%if(any(emis(:) < -9)); error('emis < -9'); end

% replace the land fovs in the emissivity
emis(1:land_nemis,~isea & iland) = land_emis(:,~isea & iland);
efreq(1:land_nemis,~isea & iland) = repmat(land_efreq(:),[1 sum(~isea & iland)]);
nemis(1,~isea & iland) = land_nemis;
%if(any(emis(:) < -9)); error('emis < -9'); end

% interpolate the rest
lf = landfrac(isea & iland); lf(lf < 0) = 0; % landfrac
emis(1:interp_nemis,isea & iland) = ...
  bsxfun(@times,interp1([0;land_efreq;9000],land_emis([1 1:end end],isea & iland),interp_efreq,'linear'), lf) + ...
  bsxfun(@times,interp1([0;sea_efreq;9000],sea_emis([1 1:end end],isea & iland),interp_efreq,'linear'), (1-lf));
efreq(1:interp_nemis,isea & iland) = repmat(interp_efreq(:),[1 sum(isea & iland)]);
nemis(1,isea & iland) = interp_nemis;
%if(any(emis(:) < -9)); error('emis < -9'); end

