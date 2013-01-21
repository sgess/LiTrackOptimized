from_sig = 0;					% ==1: calculate I_pk from sigma (uses "fac"), ==2: from Gauss-fit?, otherwise from peak of histogram 
fac = 2*pi;					% calculate I_pk from gaussian estimate
%fac = 12;						% calculate I_pk from uniform estimate

fn    = input('Litrack file name (e.g. lcls0): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end
dE_E_tol     = nprompt('Tolerable rms relative energy jitter (unitless)',1E-3);			% calculate jitter tolerance based on (e.g.) 1E-3 dE/E jitter
dIpk_Ipk_tol = nprompt('Tolerable rms relative peak-current jitter (unitless)',0.12);	% calculate jitter tolerance based on (e.g.) 12% dIpk/Ipk jitter
dT_tol       = nprompt('Tolerable rms final timing jitter (fsec)',100);					% calculate jitter tolerance based on (e.g.) 100 fs timing jitter
Ef           = nprompt('Final nominal beam energy (GeV)',14.1);                         % final e- energy [GeV]
lcls_nom     = prompt('Show LCLS nominal sensitivities','yn','y');

if lcls_nom=='y'
  load lcls_jitter_sens0.mat;
else
  Etol_dt0 = 0;
  Ztol_dt0 = 0;
  Ttol_dt0 = 0;
  Etol_dq0 = 0;
  Ztol_dq0 = 0;
  Ttol_dq0 = 0;
  Etol_dphi0 = zeros(5,1);
  Ztol_dphi0 = zeros(5,1);
  Ttol_dphi0 = zeros(5,1);
  Etol_dV0 = zeros(5,1);
  Ztol_dV0 = zeros(5,1);
  Ttol_dV0=  zeros(5,1);
end

dphi         = 0.2;						% phase error used to test sensitivity [deg, S or X]
dV_V         = 0.002;					% relative voltage error used to test sensitivity [ ]
dt    = 1;								% timing test range [psec]
dq    = 0.02;							% relative charge test [ ]

seed  = 1;								% random seed used
c     = 2.99792458E8;

eval(fnf);
beamline(:,1) = abs(beamline(:,1));     % turn off plots and sample markers
beamline0 = beamline;
Q0    = Ne*1.602e-19;

dz    = dt*1E-9*c;						% [psec -> mm]

[zpos,dE_E,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0,sigzG0,sigEG0,Ipkm0,Ipk_fit0] = litrack(fn,seed,0,0,0,beamline0);
sigE0 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk0  = Q0*1E-3*c/(sqrt(fac)*sigz);
elseif from_sig==2
  Ipk0  = Ipk_fit0;  
else
  Ipk0  = Ipkm0;
end
dd0 = 100*(E_barcuts0/Ef-1);
dt0 = z_bar0*1E-3/c*1E15;

% Do gun-timing jitter test:
% =========================
figure(1)
[zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipkm1,Ipk_fit1] = litrack(fn,seed,-dz,0,0,beamline0);
sigE1 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk1  = Q0*1E-3*c/(sqrt(fac)*sigz);
elseif from_sig==2
  Ipk1  = Ipk_fit1;  
else
  Ipk1  = Ipkm1;
end
dd1 = 100*(E_barcuts1/Ef-1);
dt1 = z_bar1*1E-3/c*1E15 - dt0;
[zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipkm2,Ipk_fit2] = litrack(fn,seed, dz,0,0,beamline0);
sigE2 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk2  = Q0*1E-3*c/(sqrt(fac)*sigz);
elseif from_sig==2
  Ipk2  = Ipk_fit2;  
else
  Ipk2  = Ipkm2;
end
dd2 = 100*(E_barcuts2/Ef-1);
dt2 = z_bar2*1E-3/c*1E15 - dt0;
[Etol_dt,mean_E] = jitter_tol([-dt 0 dt],[dd1 dd0 dd2],dE_E_tol*100,1,1);
subplot(231)
plot_polyfit([-dt 0 dt],[dd1 dd0 dd2],1,2,'\Delta{\itt}','\Delta{\itE}/{\itE}','psec','%',0,0,2)
title(sprintf('%6.3f psec (%6.3f)',Etol_dt,Etol_dt0))
enhance_plot('times',14,2,5);
[Ztol_dt,mean_I] = jitter_tol([-dt 0 dt],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
subplot(232)
plot_polyfit([-dt 0 dt],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itt}','{\itI_{pk}}','psec','kA',0,0,2)
title(sprintf('%6.3f psec (%6.3f)',Ztol_dt,Ztol_dt0))
enhance_plot('times',14,2,5);
[Ttol_dt,mean_T] = jitter_tol([-dt 0 dt],[dt1 0 dt2],dT_tol,1,1);
subplot(233)
plot_polyfit([-dt 0 dt],[dt1 0 dt2],1,2,'\Delta{\itt}','\Delta{\itt_f}','psec','fsec',0,0,2)
title(sprintf('%6.3f psec (%6.3f)',Ttol_dt,Ttol_dt0))
enhance_plot('times',14,2,5);

% Do charge jitter test:
% =====================
[zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipkm1,Ipk_fit1] = litrack(fn,seed,0,-dq,0,beamline0);
sigE1 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk1  = Q0*(1-dq)*1E-3*c/(sqrt(fac)*sigz);
elseif from_sig==2
  Ipk1  = Ipk_fit1;  
else
  Ipk1  = Ipkm1;
end
dd1 = 100*(E_barcuts1/Ef-1);
dt1 = z_bar1*1E-3/c*1E15 - dt0;
[zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipkm2,Ipk_fit2] = litrack(fn,seed,0, dq,0,beamline0);
sigE2 = 100*std(dE_E);
sigz  = std(zpos);
if from_sig==1
  Ipk2  = Q0*(1+dq)*1E-3*c/(sqrt(fac)*sigz);
elseif from_sig==2
  Ipk2  = Ipk_fit2;  
else
  Ipk2  = Ipkm2;
end
dd2 = 100*(E_barcuts2/Ef-1);
dt2 = z_bar2*1E-3/c*1E15 - dt0;
[Etol_dq,mean_E] = jitter_tol(100*[-dq 0 dq],[dd1 dd0 dd2],dE_E_tol*100,1,1);
subplot(234)
plot_polyfit(100*[-dq 0 dq],[dd1 dd0 dd2],1,2,'\Delta{\itQ}/{\itQ}','\Delta{\itE}/{\itE}','%','%',0,0,2)
title(sprintf('%6.3f %% (%6.3f)',Etol_dq,Etol_dq0))
enhance_plot('times',14,2,5);
[Ztol_dq,mean_I] = jitter_tol(100*[-dq 0 dq],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
subplot(235)
plot_polyfit(100*[-dq 0 dq],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itQ}/{\itQ}','{\itI_{pk}}','%','kA',0,0,2)
title(sprintf('%6.3f %% (%6.3f)',Ztol_dq,Ztol_dq0))
enhance_plot('times',14,2,5);
[Ttol_dq,mean_T] = jitter_tol(100*[-dq 0 dq],[dt1 0 dt2],dT_tol,1,1);
subplot(236)
plot_polyfit([-dq 0 dq],[dt1 0 dt2],1,2,'\Delta{\itQ}/{\itQ}','\Delta{\itt_f}','%','fsec',0,0,2)
title(sprintf('%6.3f %% (%6.3f)',Ttol_dq,Ttol_dq0))
enhance_plot('times',14,2,5);

% Do RF phase and voltage jitter tests:
% ====================================
i = find(beamline0(:,1)==11 & beamline0(:,2)~=0 & beamline0(:,5)~=0);	% Do only RF with V~=0 and with wake ON
n = length(i);
Etol_dphi = zeros(n,1);
Ztol_dphi = zeros(n,1);
Ttol_dphi = zeros(n,1);
Etol_dV   = zeros(n,1);
Ztol_dV   = zeros(n,1);
Ttol_dV   = zeros(n,1);
for j = 1:n
  figure(1+j);
  beamline         = beamline0;
% Do all RF phase jitter tests:
% ============================
  beamline(i(j),3) = beamline0(i(j),3) - dphi;
  [zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipkm1,Ipk_fit1] = litrack(fn,seed,0,0,0,beamline);
  sigE1 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk1  = Q0*1E-3*c/(sqrt(fac)*sigz);
  elseif from_sig==2
    Ipk1  = Ipk_fit1;  
  else
    Ipk1  = Ipkm1;
  end
  dd1 = 100*(E_barcuts1/Ef-1);
  dt1 = z_bar1*1E-3/c*1E15 - dt0;
  beamline(i(j),3) = beamline0(i(j),3) + dphi;
  [zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipkm2,Ipk_fit2] = litrack(fn,seed,0,0,0,beamline);
  sigE2 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk2  = Q0*1E-3*c/(sqrt(fac)*sigz);
  elseif from_sig==2
    Ipk2  = Ipk_fit2;  
  else
    Ipk2  = Ipkm2;
  end
  dd2 = 100*(E_barcuts2/Ef-1);
  dt2 = z_bar2*1E-3/c*1E15 - dt0;
  [Etol_dphi(j),mean_E] = jitter_tol([-dphi 0 dphi],[dd1 dd0 dd2],dE_E_tol*100,1,1);
  subplot(231)
  plot_polyfit([-dphi 0 dphi],[dd1 dd0 dd2],1,2,'\Delta{\it\phi}','\Delta{\itE}/{\itE}','deg','%',0,0,2)
  title(sprintf('%6.3f deg (%6.3f)',Etol_dphi(j),Etol_dphi0(j)))
  enhance_plot('times',14,2,5);
  [Ztol_dphi(j),mean_I] = jitter_tol([-dphi 0 dphi],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
  subplot(232)
  plot_polyfit([-dphi 0 dphi],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\it\phi}','{\itI_{pk}}','deg','kA',0,0,2)
  title(sprintf('%6.3f deg (%6.3f)',Ztol_dphi(j),Ztol_dphi0(j)))
  enhance_plot('times',14,2,5);
  [Ttol_dphi(j),mean_T] = jitter_tol([-dphi 0 dphi],[dt1 0 dt2],dT_tol,1,1);
  subplot(233)
  plot_polyfit([-dphi 0 dphi],[dt1 0 dt2],1,2,'\Delta{\it\phi}','\Delta{\itt_f}','deg','fsec',0,0,2)
  title(sprintf('%6.3f deg (%6.3f)',Ttol_dphi(j),Ttol_dphi0(j)))
  enhance_plot('times',14,2,5);

% Do all RF voltage jitter tests:
% ==============================
  beamline         = beamline0;
  beamline(i(j),2) = beamline0(i(j),2)*(1 - dV_V);
  [zpos,dE_E,E_bar1,dFWpct1,ZFWmm1,z_bar1,E_barcuts1,fcut1,sigzG1,sigEG1,Ipkm1,Ipk_fit1] = litrack(fn,seed,0,0,0,beamline);
  sigE1 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk1  = Q0*1E-3*c/(sqrt(fac)*sigz);
  elseif from_sig==2
    Ipk1  = Ipk_fit1;  
  else
    Ipk1  = Ipkm1;
  end
  dd1 = 100*(E_barcuts1/Ef-1);
  dt1 = z_bar1*1E-3/c*1E15 - dt0;
  beamline(i(j),2) = beamline0(i(j),2)*(1 + dV_V);
  [zpos,dE_E,E_bar2,dFWpct2,ZFWmm2,z_bar2,E_barcuts2,fcut2,sigzG2,sigEG2,Ipkm2,Ipk_fit2] = litrack(fn,seed,0,0,0,beamline);
  sigE2 = 100*std(dE_E);
  sigz  = std(zpos);
  if from_sig==1
    Ipk2  = Q0*1E-3*c/(sqrt(fac)*sigz);
  elseif from_sig==2
    Ipk2  = Ipk_fit2;  
  else
    Ipk2  = Ipkm2;
  end
  dd2 = 100*(E_barcuts2/Ef-1);
  dt2 = z_bar2*1E-3/c*1E15 - dt0;
  [Etol_dV(j),mean_E] = jitter_tol(100*[-dV_V 0 dV_V],[dd1 dd0 dd2],dE_E_tol*100,1,1);
  subplot(234)
  plot_polyfit(100*[-dV_V 0 dV_V],[dd1 dd0 dd2],1,2,'\Delta{\itV}/{\itV}','\Delta{\itE}/{\itE}','%','%',0,0,2)
  title(sprintf('%6.3f %% (%6.3f)',Etol_dV(j),Etol_dV0(j)))
  enhance_plot('times',14,2,5);
  [Ztol_dV(j),mean_I] = jitter_tol(100*[-dV_V 0 dV_V],[Ipk1 Ipk0 Ipk2],Ipk0*dIpk_Ipk_tol,1,1);
  subplot(235)
  plot_polyfit(100*[-dV_V 0 dV_V],[Ipk1 Ipk0 Ipk2],1,2,'\Delta{\itV}/{\itV}','{\itI_{pk}}','%','kA',0,0,2)
  title(sprintf('%6.3f %% (%6.3f)',Ztol_dV(j),Ztol_dV0(j)))
  enhance_plot('times',14,2,5);
  [Ttol_dV(j),mean_T] = jitter_tol(100*[-dV_V 0 dV_V],[dt1 0 dt2],dT_tol,1,1);
  subplot(236)
  plot_polyfit(100*[-dV_V 0 dV_V],[dt1 0 dt2],1,2,'\Delta{\itV}/{\itV}','\Delta{\itt_f}','%','fsec',0,0,2)
  title(sprintf('%6.3f %% (%6.3f)',Ttol_dV(j),Ttol_dV0(j)))
  enhance_plot('times',14,2,5);
end

if lcls_nom=='y'
	p_z = [	% 12% peak current jitter: new for 135 MeV inj., 26JUN03
		Ztol_dphi0(1) Ztol_dphi(1) 0.10	% L0 RF phase [deg]
		Ztol_dphi0(2) Ztol_dphi(2) 0.10	% L1 RF phase [deg]
		Ztol_dphi0(3) Ztol_dphi(3) 0.50	% Lx RF phase [deg]
		Ztol_dphi0(4) Ztol_dphi(4) 0.07	% L2 RF phase [deg]
		Ztol_dphi0(5) Ztol_dphi(5) 0.50	% L3 RF phase [deg]
		Ztol_dV0(1)   Ztol_dV(1) 0.10	% L0 RF voltage [%]
		Ztol_dV0(2)   Ztol_dV(2) 0.10	% L1 RF voltage [%]
		Ztol_dV0(3)   Ztol_dV(3) 0.25	% Lx RF voltage [%]
		Ztol_dV0(4)   Ztol_dV(4) 0.10	% L2 RF voltage [%]
		Ztol_dV0(5)   Ztol_dV(5) 0.50	% L3 RF voltage [%]
		Ztol_dt0(1)   Ztol_dt(1) 0.80	% gun timing [psec]
		Ztol_dq0(1)   Ztol_dq(1) 2.00	% bunch charge [%]
     								];
	
	p_E = [	% 0.10 % relative e- beam energy jitter: new for 135 MeV inj., 26JUN03
		Etol_dphi0(1) Etol_dphi(1) 0.10	% L0 RF phase [deg]
		Etol_dphi0(2) Etol_dphi(2) 0.10	% L1 RF phase [deg]
		Etol_dphi0(3) Etol_dphi(3) 0.50	% Lx RF phase [deg]
		Etol_dphi0(4) Etol_dphi(4) 0.07	% L2 RF phase [deg]
		Etol_dphi0(5) Etol_dphi(5) 0.15	% L3 RF phase [deg]
		Etol_dV0(1)   Etol_dV(1)   0.10	% L0 RF voltage [%]
		Etol_dV0(2)   Etol_dV(2)   0.10	% L1 RF voltage [%]
		Etol_dV0(3)   Etol_dV(3)   0.25	% Lx RF voltage [%]
		Etol_dV0(4)   Etol_dV(4)   0.10	% L2 RF voltage [%]
		Etol_dV0(5)   Etol_dV(5)   0.08	% L3 RF voltage [%]
		Etol_dt0(1)   Etol_dt(1)   0.80	% gun timing [psec]
		Etol_dq0(1)   Etol_dq(1)   4.00	% bunch charge [%]
										];
	
	r_z0 = p_z(:,3)./p_z(:,1);
	r_E0 = p_E(:,3)./p_E(:,1);
	r_z  = p_z(:,3)./p_z(:,2);
	r_E  = p_E(:,3)./p_E(:,2);
	
	n_z = length(r_z);
	n_E = length(r_E);
	
	disp(' ')
	disp(sprintf('Nominal sigZ budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_z0'*r_z0),n_z))
	disp(sprintf('New     sigZ budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_z'*r_z),n_z))
	disp(sprintf('Nominal dE/E budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_E0'*r_E0),n_E))
	disp(sprintf('New     dE/E budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_E'*r_E),n_E))
	disp(' ')
end