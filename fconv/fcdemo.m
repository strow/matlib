
% fcdemo - demonstrate fconv.m and fconvkc.m
%
% This demo needs the kcarta parameter file nastradB1.ip
% and the interferometric parameter file nastB1.m to work
% correctly.


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate some radiances
%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('nastradB1.f4') ~= 2

  if exist('nastradB1.ip') ~= 2
    error('can''t find kcarta driver file nastradB1.ip');
  end
  
  ! kcarta nastradB1.ip nastradB1.kc

  if exist('nastradB1.kc') ~= 2
    error('no kcarta output file nastradB1.kc')
  end

  readkc('nastradB1.kc', 'nastradB1.f4')

  if exist('nastradB1.f4') ~= 2
    error('no readkc output file nastradB1.f4');
  end
end
 
% Get interferometric parameters, including the number of
% monochromatic points in the relevant interval.
%
calcifp('nastB1');
global  v1 v2 dvk v1kind v2kind kcpts

% read the unchunked kcarta data file; in this example
% two columns of radiances
%
f = fopen('nastradB1.f4', 'r');
rkc = fread(f, [kcpts, 2], 'float');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo monochromatic convolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convolve the radiances with fconvkc(), with two apodizations
%
[rchbox, wch] = fconvkc(rkc, 'nastB1', 'box');
[rchkb6, wch] = fconvkc(rkc, 'nastB1', 'kb', 6);

% plot monochromatic and channel convolved data together
%
figure (1)
wkc = dvk*((v1kind:v2kind)-1);
plot(wkc, rad2bt(wkc, rkc(:,1)), wch, rad2bt(wch, rchkb6(:,1)))
legend('monochromatic', 'convolved, KB #6')
% print -dps fig1.ps

% plot the two different channel convolutions together
%
figure (2)
plot(wch, rchbox(:,1), wch, rchkb6(:,1))
legend('boxcar', 'KB #6')
% print -dps fig2.ps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo channel (re-)convolutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reconvolve the channel convolved data with fconv()
%
[rchkb6B, wnum2] = fconv(rchbox, wch, 'nastB1', 'kb', 6);
[rchboxB, wnum2] = fconv(rchbox, wch, 'nastB1', 'box');
[rchboxC, wnum2] = fconv(rchboxB, wch, 'nastB1', 'box');

% compare 2-step, box to KB#6, to direct KB#6
rms(rchkb6B - rchkb6) / rms(rchkb6)

% compare 2-step, boxcar to boxcar, to direct boxcar
rms(rchbox - rchboxB) / rms(rchbox)

% compare 2-step boxcar to 3-step boxcar
rms(rchboxB - rchboxC) / rms(rchboxB)

% fconv can work with subsets of the nominal channel set
% defined in the interferometric parameter file
%
[rchkb6D, wnumD] = fconv(rchbox(1:800,:), wch(1:800), 'nastB1', 'kb', 6);

% compare subset to convolution of entire interval
rms(rchkb6D - rchkb6B(1:800,:)) / rms(rchkb6B(1:800,:))

