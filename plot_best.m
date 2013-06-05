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

[a,b] = min(residual(1:500));
%[a,b] = max(I_peak);
%b = 322;
pCurrent = params(:,b-1);
pars = pCurrent;
SetPars(pars, name, nPar);

PARAM.SIMU.BIN = 256;
PARAM.INIT.NPART = 2.02e10;
PARAM.INIT.SIGZ0 = 0.00703;
PARAM.NRTL.AMPL = 0.03992;
PARAM.NRTL.PHAS = 90.32;
PARAM.NRTL.R56  = 0.6085;
PARAM.NRTL.T566  = 1.325;
PARAM.LONE.PHAS = -20.20;
PARAM.LTWO.PHAS = -3.5;
% PARAM.LI20.EHI = 0.0285;
% PARAM.LI20.R16 = 94;
PARAM.LI20.T166 = 0;
% PARAM.LI20.BETA = 4;

OUT = LiTrackOpt('FACETpar');
ip = OUT.I.PEAK(3)
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

resi = sum((ProfXLi - Line_minBG).^2);
parts = PARAM.INIT.NPART*OUT.I.PART(3)/PARAM.INIT.NESIM;
if resi < res_low; res_low = resi; end;

figure(2);
plot(line_x,Line_minBG,'b',line_x,ProfXLi,'g','linewidth',2);
xlabel('X (mm)');
legend('sYAG','LiTrack');
if savE; saveas(gca,'~/Desktop/matched_spec_fullest_pyro.png');end;

figure(3);
plot(OUT.Z.AXIS(:,nOut),OUT.Z.HIST(:,nOut),'r','linewidth',2);
xlabel('Z (mm)');
a = axis;
text(a(2)*.6,a(4)*.8,num2str(ip,'%0.2f'),'fontsize',16);
text(a(2)*.6,a(4)*.7,num2str(parts,'%0.2e'),'fontsize',16);
text(a(2)*.6,a(4)*.6,num2str(resi,'%0.1f'),'fontsize',16);
text(a(2)*.6,a(4)*.5,num2str(res_low,'%0.1f'),'fontsize',16);
if savE; saveas(gca,'~/Desktop/bunch_prof_fullest_pyro.png');end;