NMIN = 60;
NMAX = 1006;
highest = 1024;

savE = 0;
ploT = 1;

%fp = load('concat_full_pyro.mat');
fp = load('concat_1103.mat');
hp = load('concat_half_pyro.mat');
%load('more_pars_1103_59.mat');
load('more_pars_1103_389');

half = 0;
full = 1;

nOut = 3;

fp_lo_py = 389;
fp_hi_py = 59;

spec_axis = fp.cat_dat.yag_ax;
spec_thing = fp.cat_dat.YAG_SPEC(:,fp_hi_py);

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;
param_tcav;
nmax = 0;
%residual = 0;
[c,v] = min(residual);
for i = 1:1
%[a,b] = min(residual(1:nmax));
%pCurrent = params(:,b-1);
nmax = v+i-1;
pCurrent = params(:,nmax);

PARAM.INIT.SIGZ0 = pCurrent(1);   % Bunch Length
PARAM.INIT.SIGD0 = pCurrent(2);   % Initial Energy Spread
PARAM.INIT.NPART = pCurrent(3);   % Number of Particles
PARAM.INIT.ASYM  = pCurrent(4);   % Initial Gaussian Asymmetry
PARAM.NRTL.AMPL  = pCurrent(5);   % Amplitude of RF Compressor
PARAM.NRTL.PHAS  = pCurrent(6);   % RF Compressor Phase
PARAM.NRTL.R56   = pCurrent(7);   % RTL Compression
PARAM.NRTL.T566  = pCurrent(8);   % RTL Second order compression
decker           = pCurrent(9);   % 2-10 Phase
l_two            = pCurrent(10);  % 11-20 Phase
ramp             = pCurrent(11)+10.0*0.2;  % Ramp Phase
PARAM.LI20.BETA  = pCurrent(12);  % Beta Function
PARAM.LI20.R16   = pCurrent(13);  % Dispersion
PARAM.LI20.T166  = pCurrent(14);  % 2nd Order Dispersion

% Set dependent params
PARAM.LONE.PHAS = decker+ramp;  % Total PhasepCurrent(14)
PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

PARAM.LTWO.PHAS = l_two+ramp;  % Total Phase
PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain

OUT = LiTrackOpt('FACETpar');
OUT.I.PEAK(3)
xx = spec_axis';
Lineout = spec_thing;

Line_minBG = Lineout-Lineout(1);
line_x  = xx;
x_avg = mean(line_x);
[MaxLine,max_ind] = max(Line_minBG);
SumLine = sum(Line_minBG);
center = sum(line_x.*Line_minBG)/sum(Line_minBG);

SimDisp = interpSimX(OUT,line_x,PARAM.SIMU.BIN,center-x_avg);
SumX = sum(SimDisp);
normX = SumLine/SumX;
ProfXLi = normX*SimDisp;

%nmax
%residual(i) = sum(Line_minBG.*(ProfXLi - Line_minBG).^2);
%residual(i) = sum((ProfXLi - Line_minBG).^2);

if ploT
figure(1);
plot(line_x,Line_minBG,'b',line_x,ProfXLi,'g','linewidth',2);
xlabel('X (mm)');
legend('sYAG','LiTrack');
if savE; saveas(gca,'~/Desktop/matched_spec_full_pyro.png');end;

figure(2);
plot(OUT.Z.AXIS(:,nOut),OUT.Z.HIST(:,nOut),'r','linewidth',2);
xlabel('Z (mm)');
if savE; saveas(gca,'~/Desktop/bunch_prof_half_pyro.png');end;
end
end

figure(3);
plot(nmax,residual)