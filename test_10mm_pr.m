%%
clear all;
show = 1;
nOut = 3;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

param_06_04_13;
PARAM.INIT.NPART = 1.90e10;
PARAM.LI20.R56 = 0.010;
PARAM.LI20.T566 = 0.100*0.4175*(PARAM.LI20.R56/0.005)^2;
PARAM.LI20.NLO = 0.005;
PARAM.LI20.NHI = 0.015;
PARAM.LI20.ELO = -0.015;

pars_nrtl = [0.0402];
sens_nrtl = [0.06];
name_nrtl = {'NRTL AMPL'};

pars_lone = [-19.45];
sens_lone = [0.3];
name_lone = {'LONE PHAS'};

pars_ltwo = [0];
sens_ltwo = [1];
name_ltwo = {'LTWO PHAS'};

pars = [pars_nrtl; pars_lone; pars_ltwo];
sens = [sens_nrtl; sens_lone; sens_ltwo];
name = [name_nrtl; name_lone; name_ltwo];

nPar = length(pars);
[Cent, Diff, lo_lims, hi_lims] = SetParLims(pars,sens);
SetPars(pars, name, nPar);

phase_ramp = -4:0.25:4;

I_peak18 = zeros(length(phase_ramp),1);
I_peak = zeros(length(phase_ramp),1);
eeLi18 = zeros(258,length(phase_ramp));
zzLi18 = zeros(258,length(phase_ramp));
ProfELi18 = zeros(258,length(phase_ramp));
ProfZLi18 = zeros(258,length(phase_ramp));
eeLi = zeros(258,length(phase_ramp));
zzLi = zeros(258,length(phase_ramp));
ProfELi = zeros(258,length(phase_ramp));
ProfZLi = zeros(258,length(phase_ramp));
zd18 = zeros(258,258,length(phase_ramp));
zd = zeros(258,258,length(phase_ramp));

for i=1:length(phase_ramp)
    
    pars(2) = pars_lone + phase_ramp(i);
    pars(3) = pars_ltwo + phase_ramp(i);
    SetPars(pars, name, nPar);
    
    %OUT = LiTrackOpt('FACETpar');
    OUT = LiTrackOpt('FACETconNOTCH');
    
    % Calculate particle fraction
    Part_frac18 = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;
    I_peak18(i)    = OUT.I.PEAK(2);
    Part_frac = 1 - OUT.I.PART(nOut)/PARAM.INIT.NESIM;
    I_peak(i)    = OUT.I.PEAK(nOut);
    
    dE18      = OUT.E.AXIS(2,2) - OUT.E.AXIS(1,2);
    eeLi18(:,i)    = [OUT.E.AXIS(1,2)-dE18; OUT.E.AXIS(:,2); OUT.E.AXIS(end,2)+dE18];
    ProfELi18(:,i) = [0; OUT.E.HIST(:,2); 0];
    
    dZ18      = OUT.Z.AXIS(2,2) - OUT.Z.AXIS(1,2);
    zzLi18(:,i)    = [OUT.Z.AXIS(1,2)-dZ18; OUT.Z.AXIS(:,2); OUT.Z.AXIS(end,2)+dZ18];
    ProfZLi18(:,i) = [0; OUT.Z.HIST(:,2); 0];
    
    zd18(:,:,i)      = hist2(OUT.Z.DIST(1:OUT.I.PART(2),2),OUT.E.DIST(1:OUT.I.PART(2),2),zzLi18(:,i)/1000,eeLi18(:,i)/100);
    
    dE      = OUT.E.AXIS(2,nOut) - OUT.E.AXIS(1,nOut);
    eeLi(:,i)    = [OUT.E.AXIS(1,nOut)-dE; OUT.E.AXIS(:,nOut); OUT.E.AXIS(end,nOut)+dE];
    ProfELi(:,i) = [0; OUT.E.HIST(:,nOut); 0];
    
    dZ      = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
    zzLi(:,i)    = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
    ProfZLi(:,i) = [0; OUT.Z.HIST(:,nOut); 0];
    
    zd(:,:,i)      = hist2(OUT.Z.DIST(1:OUT.I.PART(nOut),nOut),OUT.E.DIST(1:OUT.I.PART(nOut),nOut),zzLi(:,i)/1000,eeLi(:,i)/100);
    
%     figure(1);
%     clf;
%     subplot(2,2,1);
%     plot(eeLi18,ProfELi18,'g','linewidth',2);
%     xlabel('\delta (%)','fontsize',14);
%     subplot(2,2,2);
%     imagesc(1000*zzLi18,eeLi18,rot90(zd18,1));
%     xlabel('Z (\mum)','fontsize',14);
%     ylabel('\delta (%)','fontsize',14);
%     subplot(2,2,3);
%     text(0.1,0.8,num2str(I_peak18(i),'%0.2f'),'fontsize',16);
%     subplot(2,2,4);
%     plot(1000*zzLi18,ProfZLi18,'r','linewidth',2);
%     xlabel('Z (\mum)','fontsize',14);
%     
%     figure(2);
%     clf;
%     subplot(2,2,1);
%     plot(eeLi,ProfELi,'g','linewidth',2);
%     xlabel('\delta (%)','fontsize',14);
%     subplot(2,2,2);
%     imagesc(1000*zzLi,eeLi,rot90(zd,1));
%     xlabel('Z (\mum)','fontsize',14);
%     ylabel('\delta (%)','fontsize',14);
%     subplot(2,2,3);
%     text(0.1,0.8,num2str(I_peak(i),'%0.2f'),'fontsize',16);
%     subplot(2,2,4);
%     plot(1000*zzLi,ProfZLi,'r','linewidth',2);
%     xlabel('Z (\mum)','fontsize',14);

end

%%

plot(phase_ramp,I_peak18,'b--*',phase_ramp,I_peak,'r--*');
xlabel('Phase Ramp (deg)','fontsize',14);
ylabel('Peak Current (kA)','fontsize',14);
l = legend('S18','S20','location','northeast');
set(l,'fontsize',14);
saveas(gca,'/Users/sgess/Desktop/presentations/figs/phase_ramp_10mm.pdf');
%%

i = 11;

figure(1);
plot(eeLi18(:,i),ProfELi18(:,i),'g','linewidth',2);
xlabel('\delta (%)','fontsize',14);
title('LI18 Energy Spectrum','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li18_spectrum_' num2str(i) '.pdf']);

figure(2);
imagesc(1000*zzLi18(:,i),eeLi18(:,i),flipud(rot90(zd18(:,:,i),1)));
xlabel('Z (\mum)','fontsize',14);
ylabel('\delta (%)','fontsize',14);
set(gca,'YDir','normal')
title('LI18 z-d Phase Space','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li18_zd_' num2str(i) '.png']);

figure(3);
plot(1000*zzLi18(:,i),ProfZLi18(:,i),'r','linewidth',2);
xlabel('Z (\mum)','fontsize',14);
title('LI18 Bunch Profile','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li18_profile_' num2str(i) '.pdf']);


figure(4);
plot(eeLi(:,i),ProfELi(:,i),'g','linewidth',2);
xlabel('\delta (%)','fontsize',14);
title('LI20 Energy Spectrum','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li20_spectrum_' num2str(i) '.pdf']);

figure(5);
imagesc(1000*zzLi(:,i),eeLi(:,i),flipud(rot90(zd(:,:,i),1)));
xlabel('Z (\mum)','fontsize',14);
ylabel('\delta (%)','fontsize',14);
set(gca,'YDir','normal')
title('LI20 z-d Phase Space','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li20_zd_' num2str(i) '.png']);

figure(6);
plot(1000*zzLi(:,i),ProfZLi(:,i),'r','linewidth',2);
xlabel('Z (\mum)','fontsize',14);
title('LI20 Bunch Profile','fontsize',14);
saveas(gca,['/Users/sgess/Desktop/presentations/figs/Li20_profile_' num2str(i) '.pdf']);