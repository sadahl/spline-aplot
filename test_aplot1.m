%% Test script for spline_aplot_dev.m (manual testing)
% FIRST PHASE testing of general adaptive algorithm spline_aplot_dev
%
% Analyzes the effect of the parameters C, R and S in different cases.
% In particular, we tweak delta, skip, tbreak and tau.
% Visual control of the methods, comparing with standard matlab plotting.
% (see test_aplot2.m for a more systematic, reliable check).

% Loading testcurves.m first (alt: initialize other splines in B-form).
run testcurves.m;
% Test splines are named: sptest0, ..., sptest4.
% randomsplines is an array holding randomly generated splines.

% Parameters to experiment with.
% Load these first, then execute one section at a time, modifying these
% parameters to observe their effect on plotting.
C=3;
R=0.5;
N=200;
delta=2;
skip=1;
consec=5;
tbreak=1;
S=[N, delta, skip, consec, tbreak];
Sigma=0;

% Note: some of the curves plotted using MATLAB standard methods (red) are
% not smooth, generally because of high curvature and/or long arc length.
% It is not necessarily an aim to make the two curves coincide, as the blue
% curve (control polygon approximation) may more truly represent the curve
% in some cases.

%% sptest0:
% Simple quadratic Bezier function.
% tauconvert=1.543;
% angletau=.3; % degrees
% tau = tauconvert*(angletau*pi/180)^2;
% S(6) = tau;
[spr, tadded, flag] = spline_aplot_dev(sptest0, C, R, S, Sigma);
% [spr, tadded] = spline_aplot(sptest0, 2);
length(tadded)
splinePlot3(sptest0, spr, tadded);

%% sptest3:
%Curve, low degree, long arc, straight sections and curved sections.
S(1)=1000;
% tauconvert=1.543;
% angletau=1; % degrees
% tau = tauconvert*(angletau*pi/180)^2;
% S(6) = tau;
S(5) = 1;
[spr, tadded, flag] = spline_aplot_dev(sptest3, C, R, S, Sigma);
length(tadded)
splinePlot3(sptest3, spr, tadded);

%% sptest2:
%Function, high degree, many variations.
S(1)=1000;
% S(3)=5;
% S(4)=10;
% S(5)=1;
% S(6)=eps(1);
[spr, tadded, flag] = spline_aplot_dev(sptest2, C, R, S, Sigma);
length(tadded);
splinePlot3(sptest2, spr, tadded);
% splinePlot3(sptest2, spr);

%% sptest1:
% Function, quintic, many variations
S(1)=100;
[spr, tadded, flag] = spline_aplot_dev(sptest1, C, R, S, Sigma);
length(tadded)
splinePlot3(sptest1, spr, tadded);

%% sptest4:
% Curve, containing a small curl
S(1)=300; %C=4;
[spr, tadded, flag] = spline_aplot_dev(sptest4, C, R, S, Sigma);
length(tadded);
splinePlot3(sptest4, spr);
