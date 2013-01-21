dE_E_tol     = 1E-3;			% calculate jitter tolerance based on (e.g.) 1E-3 dE/E jitter
dIpk_Ipk_tol = 0.10;			% calculate jitter tolerance based on (e.g.) 10% dIpk/Ipk jitter

seed  = 1;						% random seed used
c     = 2.99792458E8;

fn    = input('Litrack file name (e.g. lcls): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end

rg = prompt('Calculate peak current from RMS or from Gaussian-fitted bunch length?','rg','r');

eval(fnf);
beamline(:,1) = abs(beamline(:,1))     % turn off plots and sample markers

Ef    = nprompt('Final nominal energy in GeV',14.346,0,1000);

ynt   = prompt('Do timing jitter test?','yn','y');
ynq   = prompt('Do charge jitter test?','yn','y');

if ynt=='y'
  dt    = nprompt('Maximum abs. of timing scan (in psec)',2.0,0,1000);
  n     = nprompt('Number of ODD timing points to take',3,3,99);
  if rem(n,2)==0
    error('Must enter ODD number of points to scan - quitting.')
  end
  dz    = abs(dt)*1E-9*c;			% [psec -> mm]
  z0    = -dz:(2*dz/(n-1)):dz;			% varied offset [mm]
  t0    = z0*1E9/c;
end

if ynq == 'y'
  dq    = nprompt('Maximum abs. of relative charge scan (in %)',2.0,0,10);
  nq    = nprompt('Number of ODD charge points to take',3,3,99);
  if rem(nq,2)==0
    error('Must enter ODD number of points to scan - quitting.')
  end
  dQ_Q  = (-dq:(2*dq/(nq-1)):dq)/100;			% varied charge [ ]
end

E_bar  = zeros(1,n);
E_barcuts  = zeros(1,n);
dFWpct = zeros(1,n);
ZFWmm  = zeros(1,n);
z_bar  = zeros(1,n);
fcut   = zeros(1,n);
sigE   = zeros(1,n);
sigz   = zeros(1,n);

if ynq == 'y'
  E_barq  = zeros(1,nq);
  E_barcutsq  = zeros(1,nq);
  dFWpctq = zeros(1,nq);
  ZFWmmq  = zeros(1,nq);
  z_barq  = zeros(1,nq);
  fcutq   = zeros(1,nq);
  sigEq   = zeros(1,nq);
  sigzq   = zeros(1,nq);
  Ipkq    = zeros(1,nq);
end

j = (n+1)/2;
[zpos,dE_E,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0,sigzG0] = ...
  LiTrack(fn,seed,0,0,0,beamline);
%zposj = zeros(size(zpos));
%dE_Ej = zeros(size(dE_E));
%zposj(:,j) = zpos;
%dE_Ej(:,j) = dE_E;
sigE(j) = 100*std(dE_E);
if rg == 'r'
  sigz(j) = 1E3*std(zpos);	% for LCLS
else
  sigz(j) = 1E3*sigzG0;		% for NLC
end
dFWpct(j) = dFWpct0;
ZFWmm(j) = ZFWmm0;
z_bar(j) = z_bar0;
E_barcuts(j) = E_barcuts0;
fcut(j) = fcut0;

%zposjq  = zeros(size(zpos));
%dE_Ejq  = zeros(size(dE_E));
%zposjq(:,j) = zpos;
%dE_Ejq(:,j) = dE_E;
sigEq(j) = 100*std(dE_E);
sigzq(j) = 1E3*std(zpos);
Ipkq(j)  = 1/sigzq(j);
dFWpctq(j) = dFWpct0;
ZFWmmq(j) = ZFWmm0;
z_barq(j) = z_bar0;
E_barcutsq(j) = E_barcuts0;
fcutq(j) = fcut0;

if ynt=='y'
  for j = 1:n
    if j~=(n+1)/2
      [zpos,dE_E,E_bar(j),dFWpct(j),ZFWmm(j),z_bar(j),E_barcuts(j),fcut(j),sigzG(j)] = ...
        LiTrack(fn,seed,z0(j),0,0,beamline);
%      zposj(:,j) = zpos;
%      dE_Ej(:,j) = dE_E;
      sigE(j) = 100*std(dE_E);
      if rg == 'r'
        sigz(j) = 1E3*std(zpos);	% for LCLS
      else
        sigz(j) = 1E3*sigzG(j);		% for NLC
      end
    end
  end
  sigz0  = sigz(n-(n-1)/2);
  dd = 100*(E_barcuts/Ef-1);
  [Etol,mean_E] = jitter_tol(t0,dd,dE_E_tol*100,1,1);
  [Ztol,mean_I] = jitter_tol(t0,sigz,sigz0*dIpk_Ipk_tol,1,1);

  figure(1)

  clf;
  subplot;
  subplot(221)
  plot_spline(t0,dd)
  hold on
  plot(t0,dd,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dt dt]);
%  xlabel('\Delta{\itt}_0 /psec')
  ylabel('\langle{\it\DeltaE}/{\itE}_0\rangle /%')
  title(['{\itE}_{0} = ' sprintf('%7.3f GeV, tol=%4.2f psec',Ef,Etol)])
  enhance_plot

  dEw = sigE;
  subplot(222)
  plot_spline(t0,dEw)
  hold on
  plot(t0,dEw,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dt dt]);
%  xlabel('\Delta{\itt}_0 /psec')
  ylabel('({\it\DeltaE}/{\itE}_0)_{rms} /%')
  enhance_plot

  zw = sigz;
  subplot(223)
  plot_spline(t0,zw)
  hold on
  plot(t0,zw,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dt dt]);
  xlabel('\Delta{\itt}_0 /psec')
  ylabel('{\it\sigma_z} /mm')
  title(sprintf('tol=%4.2f psec',Ztol))
  enhance_plot

  z_nom = z_bar((n + 1)/2);
  subplot(224)
  plot_spline(t0,1E9*(z_bar - z_nom)/c)
  hold on
  plot(t0,1E9*(z_bar - z_nom)/c,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dt dt]);
  xlabel('\Delta{\itt}_0 /psec')
  ylabel('\langle\Delta{\itt_f}\rangle /psec')
  enhance_plot

  Ha = gca;
  Hy = get(Ha,'YLabel');
  set(Hy,'VerticalAlignment','baseline');
end

if ynq=='y'
  for j = 1:nq
    if j~=(nq+1)/2
      [zposq,dE_Eq,E_barq(j),dFWpctq(j),ZFWmmq(j),z_barq(j),E_barcutsq(j),fcutq(j)] = ...
        LiTrack(fn,seed,0,dQ_Q(j),0,beamline);
%      zposjq(:,j) = zposq;
%      dE_Ejq(:,j) = dE_Eq;
      sigEq(j) = 100*std(dE_Eq);
      sigzq(j) = 1E3*std(zposq);
      Ipkq(j)  = (1+dQ_Q(j))/sigzq(j);
    end
  end
  sigz0q = sigzq(nq-(nq-1)/2);
  Ipk0q  = 1/sigz0q;
  figure(2)
  ddq = 100*(E_barcutsq/Ef-1);

  
  [Etolq,mean_Eq] = jitter_tol(100*dQ_Q,ddq,dE_E_tol*100,1,1);
  [Ztolq,mean_Iq] = jitter_tol(100*dQ_Q,Ipkq,Ipk0q*dIpk_Ipk_tol,1,1);

  figure(2);
  subplot;
  subplot(221)
  plot_spline(100*dQ_Q,ddq)
  hold on
  plot(100*dQ_Q,ddq,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dq dq]);
%  xlabel('\Delta{\itQ}/{\itQ}_0 /%')
  ylabel('\langle{\it\DeltaE}/{\itE}_0\rangle /%')
  title(['{\itE}_{0} = ' sprintf('%7.3f GeV, tol=%4.2f %%',Ef,Etolq)])
  enhance_plot

  dEw = sigEq;
  subplot(222)
  plot_spline(100*dQ_Q,dEw)
  hold on
  plot(100*dQ_Q,dEw,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dq dq]);
%  xlabel('\Delta{\itQ}/{\itQ}_0 /%')
  ylabel('({\it\DeltaE}/{\itE}_0)_{rms} /%')
  enhance_plot

  zw = Ipkq/Ipk0q;
  subplot(223)
  plot_spline(100*dQ_Q,zw)
  hold on
  plot(100*dQ_Q,zw,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dq dq]);
  xlabel('\Delta{\itQ}/{\itQ}_0 /%')
  ylabel('{\itI_{pk}}/{\itI_{pk}}_0')
  title(sprintf('tol=%4.2f %%',Ztolq))
  enhance_plot

  z_nomq = z_barq((nq + 1)/2);
  subplot(224)
  plot_spline(100*dQ_Q,1E9*(z_barq - z_nomq)/c)
  hold on
  plot(100*dQ_Q,1E9*(z_barq - z_nomq)/c,'or')
  hold off
  h = gca;
  set(h,'XLim',[-dq dq]);
  xlabel('\Delta{\itQ}/{\itQ}_0 /%')
  ylabel('\langle\Delta{\itt_f}\rangle /psec')
  enhance_plot

  Ha = gca;
  Hy = get(Ha,'YLabel');
  set(Hy,'VerticalAlignment','baseline');
end