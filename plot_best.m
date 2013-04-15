load('concat_1103.mat');
load('retry_1103.mat');

global A;
A = load('slac.dat');

global PARAM;
param_tcav;

nmin = 60;
nmax = 2000;

nOut = 3;
savE = 0;

spec_axis = cat_dat.yag_ax;
spec_thing = cat_dat.YAG_SPEC(:,59);

[a,b] = min(residual(1:nmax));
pCurrent = params(:,b-1);

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
ramp             = pCurrent(11);  % Ramp Phase
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


figure(1);
plot(line_x,Line_minBG,'b',line_x,ProfXLi,'g','linewidth',2);
xlabel('X (mm)');
legend('sYAG','LiTrack');
if savE; saveas(gca,'~/Desktop/matched_spec_full_pyro.png');end;

figure(2);
plot(OUT.Z.AXIS(:,nOut),OUT.Z.HIST(:,nOut),'r','linewidth',2);
xlabel('Z (mm)');
if savE; saveas(gca,'~/Desktop/bunch_prof_full_pyro.png');end;