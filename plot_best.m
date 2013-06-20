%load('concat_1103.mat');
%load('retry_1103.mat');
%load('more_pars_1103_59.mat');

global A;
A = load('slac.dat');

global PARAM;

nmin = 1;
nmax = 200;

nOut = 3;
savE = 0;
old = 0;

[a,b] = min(residual(1:nmax));
pCurrent = params(:,b-1);
pars = pCurrent;
SetPars(pars, name, nPar);

load('~/Desktop/param_420.mat');
% PARAM.INIT.SIGZ0 = 0.0060;
% PARAM.INIT.SIGD0 = 6.5e-4;
% PARAM.INIT.NPART = 1.95e10;
% PARAM.INIT.ASYM  = -0.20;
% % 
% PARAM.NRTL.AMPL = 0.0403;
% PARAM.NRTL.PHAS = 89.9;
% 
% PARAM.LONE.RAMP = 0;
% PARAM.LONE.DECK = -20.5;
% PARAM.LONE.PHAS = PARAM.LONE.RAMP+PARAM.LONE.DECK;
% PARAM.LTWO.PHAS = PARAM.LONE.RAMP;
% PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
% PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain
% 
% 
% PARAM.LI20.R56  = 0.005;
% PARAM.LI20.T566 = 0.000;
% PARAM.LI20.R16  = 85;
% PARAM.LI20.T166 = 000;


OUT = LiTrackOpt('FACETpar');
ip = OUT.I.PEAK(3);

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

dZ = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
zzLi = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
ProfZLi = [0; OUT.Z.HIST(:,nOut); 0];

dE = OUT.E.AXIS(2,nOut) - OUT.E.AXIS(1,nOut);
eeLi = [OUT.E.AXIS(1,nOut)-dE; OUT.E.AXIS(:,nOut); OUT.E.AXIS(end,nOut)+dE];
ProfELi = [0; OUT.E.HIST(:,nOut); 0];
zd = hist2(OUT.Z.DIST(1:OUT.I.PART(nOut),nOut),OUT.E.DIST(1:OUT.I.PART(nOut),nOut),zzLi/1000,eeLi/100);

resi = sum((ProfXLi - Line_minBG).^2);
parts = PARAM.INIT.NPART*OUT.I.PART(3)/PARAM.INIT.NESIM;
if ~exist('res_low','var'); res_low = 1e10; end;
if resi < res_low; res_low = resi; end;

figure(2);
plot(line_x,Line_minBG,'b',line_x,ProfXLi,'g','linewidth',2);
xlabel('X (mm)');
legend('sYAG','LiTrack');
if savE; saveas(gca,'~/Desktop/matched_spec_fullest_pyro.png');end;

figure(3);
plot(zzLi,ProfZLi,'r','linewidth',2);
xlabel('Z (mm)');
a = axis;
text(a(2)*.6,a(4)*.8,num2str(ip,'%0.2f'),'fontsize',16);
text(a(2)*.6,a(4)*.7,num2str(parts,'%0.2e'),'fontsize',16);
text(a(2)*.6,a(4)*.6,num2str(resi,'%0.1f'),'fontsize',16);
text(a(2)*.6,a(4)*.5,num2str(res_low,'%0.1f'),'fontsize',16);
if savE; saveas(gca,'~/Desktop/bunch_prof_fullest_pyro.png');end;

figure(4);
imagesc(1000*zzLi,eeLi,flipud(rot90(zd,1)));
xlabel('Z (\mum)','fontsize',14);
ylabel('\delta (%)','fontsize',14);
set(gca,'YDir','normal');