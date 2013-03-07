
% apdemo -- Demonstration of Apodization Functions
%
% Interferometric parameters are set in assignments, below.
% The most important parameter here is L1, the apodization
% limit.  The value chosen for L2 determines dv = 1/(2*L2),
% and so the amount of oversampling and the resolution of the
% plotted functions.  The value of v1 should be large enough
% to include significant wings of the spectral response
% functions of interest.
%
% See apod.m for the set of valid apodization types.

% H. Motteler, 8/24/98


L1    = input('path length > ');
atype = input('apodization type > ', 's');
aparg = input('optional parameter > ');

% interferometric parameters

L2 = 20;		% total path length (cm)
v1 = 50;		% max transform wavenumber (1/cm)

dv = 1/(2*L2);		% wavenumber increment (1/cm)

v1ind = round(v1/dv) + 1;

cospts = 2^nextpow2(v1ind) + 1;

vmax = dv * (cospts-1); 

dd = 1/(2*vmax);	% optical path increment (cm)

Lmax = dd * (cospts-1); 

L1ind = round(L1/dd) + 1;    
L2ind = round(L2/dd) + 1;    

L1a = dd * (L1ind-1); 
L2a = dd * (L2ind-1); 

L1pts = 0:dd:L1a; 
L2pts = 0:dd:L2;
v1pts = 0:dv:v1;

fprintf(1, 'v1=%g, vmax=%g, dv=%g\n', v1, vmax, dv);
fprintf(1, 'L1=%g, L2=%g, Lmax=%g, dd=%g\n', L1, L2, Lmax, dd);
fprintf(1, 'L1ind=%d, L2ind=%d, v1ind=%d, cospts=%d\n', ...
		L1ind, L2ind, v1ind, cospts);

% defined apodization to calculated response fnx

intf1 = zeros(cospts,1);
intf1(1:L1ind) = apod(L1pts', L1, atype, aparg);
spec2 = real(fft([intf1; flipud(intf1(2:cospts-1,1))]));
spec2 = spec2(1:cospts) / max(spec2(1:cospts));

% defined response fnx to calculated apodization

spec1 = zeros(cospts,1);
spec1(1:v1ind) = resp(v1pts', L1, atype, aparg);
intf2 = real(ifft([spec1; flipud(spec1(2:cospts-1,1))]));
intf2 = intf2(1:L2ind) / max(intf2(1:L2ind));


% apodization plot

dplot = 2.5;			% plot limit (cm)
dplotind = round(dplot/dd) + 1;

x  = L2pts(1:dplotind);
y1 = intf1(1:dplotind);
y2 = intf2(1:dplotind);

figure (1); 
h1 = subplot(3,1,1);
p = get(h1, 'Position'); p(2) = p(2)+.03; set(h1, 'Position', p);

plot(x, y1, x, y2);
title([apname(atype, aparg), ' apodization, L = ', num2str(L1), ' cm']);
xlabel('cm');
legend('defined apodization', 'apodization from response function', 1);
grid

% response function plot

vplot = 3;			% plot limit (1/cm)
vplotind = round(vplot/dv) + 1;

[t,i] =  min(abs(spec1(1:vplotind) - 0.5));
fwhm = 2*v1pts(i);
fprintf(1, 'FWHM = %g\n', fwhm);

x  = v1pts(1:vplotind);
y1 = spec1(1:vplotind);
y2 = spec2(1:vplotind);

figure (1)
h2 = subplot(3,1,2);
plot(x, y1, x, y2);
title([apname(atype, aparg), ' response function, L = ', num2str(L1), ' cm'])
xlabel('1/cm')
legend('defined response function', 'response function from apodization', 1);
text(fwhm/1.9, 0.5, [' FWHM = ',num2str(fwhm)])
grid

% response function zoom

zv0 = 0.5;
zv1 = 10;
zmin = -0.05;
zmax =  0.05;

totmin = sum(spec1(find(spec1 < 0)));
minspec = min(spec1);
zmin = min(minspec, zmin);
zmax = max(-minspec, zmax);

vplotind0 = round(zv0/dv) + 1;
vplotind1 = round(zv1/dv) + 1;

x  = v1pts(vplotind0:vplotind1);
y1 = spec1(vplotind0:vplotind1);

figure (1)
h3 = subplot(3,1,3);
p = get(h3, 'Position'); p(2) = p(2)-.03; set(h3, 'Position', p);

plot(x, y1);
title([apname(atype, aparg),' response (detail)']);
xlabel('1/cm')
axis([zv0, zv1, zmin, zmax]);
legend('defined response function', 1);
text(zv0+1, zmin/1.5, ...
	[' minima = ',num2str(fixpt(minspec,4))]);
text((zv0+zv1)/2, zmin/1.5, ...
	['total negative = ',num2str(fixpt(totmin,2))]);
grid

eval(['print -dpsc ', atype, num2str(aparg), 'L', ...
		num2str(fixpt(L1,1)*10), '.ps'])

