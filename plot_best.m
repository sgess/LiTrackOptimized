%load('concat_1103.mat');
%load('retry_1103.mat');
%load('more_pars_1103_59.mat');

global A;
A = load('slac.dat');

global PARAM;
%param_tcav;

nmin = 1;
nmax = 500;

nOut = 3;
savE = 0;
old = 0;

%spec_axis = cat_dat.yag_ax;
%spec_thing = cat_dat.YAG_SPEC(:,59);



if old
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

PARAM.INIT.SIGZ0 = 0.0075;   % Bunch Length
PARAM.INIT.SIGD0 = 0.0008;   % Initial Energy Spread
PARAM.INIT.NPART = 1.95e10;   % Number of Particles
PARAM.INIT.ASYM  = -0.15;   % Initial Gaussian Asymmetry
PARAM.NRTL.AMPL  = 0.04000;   % Amplitude of RF Compressor
PARAM.NRTL.PHAS  = 89.75;   % RF Compressor Phase
PARAM.NRTL.R56   = 0.602;   % RTL Compression
PARAM.NRTL.T566  = 1.3;   % RTL Second order compression
decker           = -19.9;   % 2-10 Phase
l_two            = 4;  % 11-20 Phase
ramp             = 0;  % Ramp Phase
PARAM.LI20.BETA  = 5;  % Beta Function
PARAM.LI20.R16   = 94;  % Dispersion
PARAM.LI20.T166  = 700;  % 2nd Order Dispersion
PARAM.LI20.EHI   = 0.010;
PARAM.LI20.ELO   = -0.020;

PARAM.LI20.R56   = 0.005;  
PARAM.LI20.T566  = 0.100;

% Set dependent params
PARAM.LONE.PHAS = decker+ramp;  % Total PhasepCurrent(14)
PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

PARAM.LTWO.PHAS = l_two+ramp;  % Total Phase
PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain
end


% param_04_16_13;
% 
% pars_init = [0.0066;        0.0008;         2.0e10;       -0.15];
% sens_init = [0.2;           0.2;            0.1;           0.3];
% name_init = {'INIT SIGZ0'; 'INIT SIGD0'; 'INIT NPART'; 'INIT ASYM'};
% 
% pars_nrtl = [0.0405;        90.00;         0.602;       1.3];
% sens_nrtl = [0.06;           0.01;            0.02;       0.1];
% name_nrtl = {'NRTL AMPL'; 'NRTL PHAS'; 'NRTL R56'; 'NRTL T566'};
% 
% pars_lone = [-23.0];
% sens_lone = [0.3];
% name_lone = {'LONE PHAS'};
% 
% pars_ltwo = [1];
% sens_ltwo = [1];
% name_ltwo = {'LTWO PHAS'};
% 
% pars_li20 = [1;             95;             50;        0.030;   -0.030;];
% sens_li20 = [0.5;           0.01;            1;         0.5;     0.5];
% name_li20 = {'LI20 BETA'; 'LI20 R16'; 'LI20 T166'; 'LI20 EHI'; 'LI20 ELO'};
% 
% pars = [pars_init; pars_nrtl; pars_lone; pars_ltwo; pars_li20];
% sens = [sens_init; sens_nrtl; sens_lone; sens_ltwo; sens_li20];
% name = [name_init; name_nrtl; name_lone; name_ltwo; name_li20];
% 
% nPar = length(pars);
% [Cent, Diff, lo_lims, hi_lims] = SetParLims(pars,sens);

[a,b] = min(residual(1:nmax));
%[a,b] = max(I_peak);
%b = 322;
pCurrent = params(:,b-1);
pars = pCurrent;
SetPars(pars, name, nPar);

PARAM.LI20.T166  = 0;
%PARAM.NRTL.AMPL = 0.0398;
%PARAM.LONE.PHAS = -24;
%PARAM.LI20.R56 = 0.005;
%PARAM.LI20.T566 = 0.1;

OUT = LiTrackOpt('FACETpar');
OUT.I.PEAK(3)
xx = spec_axis;
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


figure(2);
plot(line_x,Line_minBG,'b',line_x,ProfXLi,'g','linewidth',2);
xlabel('X (mm)');
legend('sYAG','LiTrack');
if savE; saveas(gca,'~/Desktop/matched_spec_fullest_pyro.png');end;

figure(3);
plot(OUT.Z.AXIS(:,nOut),OUT.Z.HIST(:,nOut),'r','linewidth',2);
xlabel('Z (mm)');
if savE; saveas(gca,'~/Desktop/bunch_prof_fullest_pyro.png');end;