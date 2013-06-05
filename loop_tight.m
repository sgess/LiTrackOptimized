clear all;

file_name = '~/Desktop/scan.mat';
savE = 1;
scan = 1;
test = 0;
global PARAM;
param_04_16_13;
PARAM.LI20.R16   = 85;
PARAM.LI20.BETA  = 3.0;
PARAM.LI20.T166  = 0;
PARAM.NRTL.R56 = 0.603;
PARAM.NRTL.T566 = 1.3;

global A;
A = load('slac.dat');

%sigz0 range
sigmaz_lo = 0.0055;
sigmaz_hi = 0.007;

% particle range
npart_lo = 2.0e10;
npart_hi = 2.2e10;

% comp phase range
nrtlPhas_lo = 89.7;
nrtlPhas_hi = 90.7;

% comp phase range
nrtlAmpl_lo = 0.0405;
nrtlAmpl_hi = 0.0415;

%2-10 range
phasOne_lo = -22;
phasOne_hi = -20.5;

%11-20 phase
phasTwo_lo = -8;
phasTwo_hi = -4;

% number of sample points
sigmaz_el   = 9;
npart_el    = 9;
nrtlPhas_el = 9;
nrtlAmpl_el = 9;
phasOne_el  = 9;
phasTwo_el  = 9;


