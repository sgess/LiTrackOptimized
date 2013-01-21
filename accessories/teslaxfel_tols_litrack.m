from_sig = 0;					% ==1: calculate I_pk from sigma, otherwise from peak polyfit 
%fac = 2*pi;					% calculate I_pk from gaussian estimate
fac = 12;						% calculate I_pk from uniform estimate

dE_E_tol     = 1E-3;			% calculate jitter tolerance based on (e.g.) 1E-3 dE/E jitter
dIpk_Ipk_tol = 0.12;			% calculate jitter tolerance based on (e.g.) 10% dIpk/Ipk jitter
dphi         = 0.2;				% phase error used to test sensitivity [deg, S or X]
dV_V         = 0.002;			% relative voltage error used to test sensitivity [ ]

seed  = 1;						% random seed used
c     = 2.99792458E8;

fn    = input('Litrack file name (e.g. teslaxfel): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end

eval(fnf);
beamline(:,1) = abs(beamline(:,1))     % turn off plots and sample markers
beamline0 = beamline;

Ef    = 20.5;
dt    = 2;								% timing test range [psec]
dq    = 0.02;							% relative charge test [ ]

dz    = dt*1E-9*c;						% [psec -> mm]

[zpos,dE_E,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0,sigzG0,sigEG0,Ipk_fit0] = litrack(fn,seed,0,0,0,beamline0);
sigE0 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk0  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
else
  Ipk0  = Ipk_fit0;
end
dd0 = 100*(E_barcuts0/Ef-1);

% Do gun-timing jitter test:
% =========================
figure(1)
[zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipk_fit1] = litrack(fn,seed,-dz,0,0,beamline0);
sigE1 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk1  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
else
  Ipk1  = Ipk_fit1;
end
dd1 = 100*(E_barcuts1/Ef-1);
[zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipk_fit2] = litrack(fn,seed, dz,0,0,beamline0);
sigE2 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk2  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
else
  Ipk2  = Ipk_fit2;
end
dd2 = 100*(E_barcuts2/Ef-1);
[Etol_dt,mean_E] = jitter_tol([-dt 0 dt],[dd1 dd0 dd2],dE_E_tol*100,1,1);
subplot(221)
plot_polyfit([-dt 0 dt],[dd1 dd0 dd2],1,2,'\Delta{\itt}','\Delta{\itE}/{\itE}','psec','%')
title(sprintf('%5.2f psec',Etol_dt))
enhance_plot;
[Ztol_dt,mean_I] = jitter_tol([-dt 0 dt],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
subplot(223)
plot_polyfit([-dt 0 dt],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itt}','{\itI_{pk}}','psec','kA')
title(sprintf('%5.2f psec',Ztol_dt))
enhance_plot;

% Do charge jitter test:
% =====================
[zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipk_fit1] = litrack(fn,seed,0,-dq,0,beamline0);
sigE1 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk1  = 1e-9*(1-dq)*1E-3*c/(sqrt(fac)*sigz);
else
  Ipk1  = Ipk_fit1;
end
dd1 = 100*(E_barcuts1/Ef-1);
[zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipk_fit2] = litrack(fn,seed,0, dq,0,beamline0);
sigE2 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk2  = 1e-9*(1+dq)*1E-3*c/(sqrt(fac)*sigz);
else
  Ipk2  = Ipk_fit2;
end
dd2 = 100*(E_barcuts2/Ef-1);
[Etol_dq,mean_E] = jitter_tol(100*[-dq 0 dq],[dd1 dd0 dd2],dE_E_tol*100,1,1);
subplot(222)
plot_polyfit(100*[-dq 0 dq],[dd1 dd0 dd2],1,2,'\Delta{\itQ}/{\itQ}','\Delta{\itE}/{\itE}','%','%')
title(sprintf('%5.2f %%',Etol_dq))
enhance_plot;
[Ztol_dq,mean_I] = jitter_tol(100*[-dq 0 dq],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
subplot(224)
plot_polyfit(100*[-dq 0 dq],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itQ}/{\itQ}','{\itI_{pk}}','%','kA')
title(sprintf('%5.2f %%',Ztol_dq))
enhance_plot;

% Do RF phase and voltage jitter tests:
% ====================================
i = find(beamline0(:,1)==11 & beamline0(:,2)~=0 & beamline0(:,5)~=0);	% Do only RF with V~=0 and with wake ON
n = length(i);
Etol_dphi = zeros(n,1);
Ztol_dphi = zeros(n,1);
Etol_dV   = zeros(n,1);
Ztol_dV   = zeros(n,1);
for j = 1:n
  figure(1+j);
  beamline         = beamline0;
% Do all RF phase jitter tests:
% ============================
  beamline(i(j),3) = beamline0(i(j),3) - dphi;
  [zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipk_fit1] = litrack(fn,seed,0,0,0,beamline);
  sigE1 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk1  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
  else
    Ipk1  = Ipk_fit1;
  end
  dd1 = 100*(E_barcuts1/Ef-1);
  beamline(i(j),3) = beamline0(i(j),3) + dphi;
  [zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipk_fit2] = litrack(fn,seed,0,0,0,beamline);
  sigE2 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk2  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
  else
    Ipk2  = Ipk_fit2;
  end
  dd2 = 100*(E_barcuts2/Ef-1);
  [Etol_dphi(j),mean_E] = jitter_tol([-dphi 0 dphi],[dd1 dd0 dd2],dE_E_tol*100,1,1);
  subplot(221)
  plot_polyfit([-dphi 0 dphi],[dd1 dd0 dd2],1,2,'\Delta{\it\phi}','\Delta{\itE}/{\itE}','deg','%')
  title(sprintf('%5.2f deg',Etol_dphi(j)))
  enhance_plot;
  [Ztol_dphi(j),mean_I] = jitter_tol([-dphi 0 dphi],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
  subplot(223)
  plot_polyfit([-dphi 0 dphi],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\it\phi}','{\itI_{pk}}','deg','kA')
  title(sprintf('%5.2f deg',Ztol_dphi(j)))
  enhance_plot;
% Do all RF voltage jitter tests:
% ==============================
  beamline         = beamline0;
  beamline(i(j),2) = beamline0(i(j),2)*(1 - dV_V);
  [zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipk_fit1] = litrack(fn,seed,0,0,0,beamline);
  sigE1 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk1  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
  else
    Ipk1  = Ipk_fit1;
  end
  dd1 = 100*(E_barcuts1/Ef-1);
  beamline(i(j),2) = beamline0(i(j),2)*(1 + dV_V);
  [zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipk_fit2] = litrack(fn,seed,0,0,0,beamline);
  sigE2 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk2  = 1e-9*1E-3*c/(sqrt(fac)*sigz);
  else
    Ipk2  = Ipk_fit2;
  end
  dd2 = 100*(E_barcuts2/Ef-1);
  [Etol_dV(j),mean_E] = jitter_tol(100*[-dV_V 0 dV_V],[dd1 dd0 dd2],dE_E_tol*100,1,1);
  subplot(222)
  plot_polyfit(100*[-dV_V 0 dV_V],[dd1 dd0 dd2],1,2,'\Delta{\itV}/{\itV}','\Delta{\itE}/{\itE}','%','%')
  title(sprintf('%5.2f %%',Etol_dV(j)))
  enhance_plot;
  [Ztol_dV(j),mean_I] = jitter_tol(100*[-dV_V 0 dV_V],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
  subplot(224)
  plot_polyfit(100*[-dV_V 0 dV_V],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itV}/{\itV}','{\itI_{pk}}','%','kA')
  title(sprintf('%5.2f %%',Ztol_dV(j)))
  enhance_plot;
end
