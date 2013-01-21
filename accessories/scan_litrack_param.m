seed  = 1;

fn     = input('Litrack file name (e.g. lcls): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end

eval(fnf)

Ef     = nprompt('Final nominal energy in GeV',28.5,0.01,1000);
pname  = input('Name of scanned parameter (e.g. BC1 Voltage): ','s');
punit  = input('Unit of scanned parameter: ','s');
beamline(:,1) = abs(beamline(:,1))     % turn off plots and sample markers
rblim = length(beamline(:,1));
cblim = length(beamline(1,:));
rb     = nprompt('Scan parameter in what beamline row',5,1,rblim);
cb     = nprompt('Scan parameter in what beamline column',3,1,cblim);
disp(' ')
disp(sprintf(['Value of parameter selected = %7.3e ' punit],beamline(rb,cb)))
disp(' ')
dp     = nprompt(['Absolute value of maximum CHANGE to ' pname ' (in ' punit ')']);
n      = nprompt('Number of ODD points to take',11);
if rem(n,2)==0
  error('Must enter ODD number of points to scan - quitting.')
end
if n <= 1
  error('Need at least two points')
end
param  = beamline(rb,cb) + (-dp:(2*dp/(n-1)):dp);	% varied parameter

disp(' ')
disp([pname ' will be scanned as:'])
disp(sprintf(['%7.3g ' punit],param))
disp(' ')

c      = 2.99792458E8;
E_bar  = zeros(1,n);
E_barcuts  = zeros(1,n);
dFWpct = zeros(1,n);
ZFWmm  = zeros(1,n);
z_bar  = zeros(1,n);
fcut   = zeros(1,n);
sigE   = zeros(1,n);
sigz   = zeros(1,n);
sigzg  = zeros(1,n);

for j = 1:n
  disp(['==> Run #' int2str(j)])
  [zpos,dE_E,E_bar(j),dFWpct(j),ZFWmm(j),z_bar(j),E_barcuts(j),fcut(j),sigzg(j)] = ...
                       LiTrack(fn,seed,0,0,[param(j) rb cb],beamline);
  sigE(j) = 100*std(dE_E);
  sigz(j) = 1E3*std(zpos);
end
sigz0  = sigz(n-(n-1)/2);
dd = 100*(E_barcuts/Ef-1);
[Etol,mean_E] = jitter_tol((param-param((n+1)/2)),dd,0.10,1,1);
[Ztol,mean_I] = jitter_tol((param-param((n+1)/2)),sigz,0.12*sigz0,1,1);
param_min = min(param);
param_max = max(param);

clf;
subplot;
subplot(221)
%plot_spline(param,dd)
hold on
plot(param,dd,'or')
hold off
h = gca;
set(h,'XLim',[param_min param_max]);
ylabel('\langle{\it\DeltaE}/{\itE}_0\rangle /%')
title(['{\itE}_{0} = ' sprintf(['%7.3f GeV, tol=%7.3g ' punit],Ef,Etol)])
enhance_plot

dEw = sigE;
subplot(222)
%plot_spline(param,dEw)
hold on
plot(param,dEw,'or')
hold off
h = gca;
set(h,'XLim',[param_min param_max]);
ylabel('({\it\DeltaE}/{\itE}_0)_{rms} /%')
enhance_plot

zw = sigz;
subplot(223)
%plot_spline(param,zw)
hold on
plot(param,zw,'or')
plot(param,sigzg,'gd')
plot(param,ZFWmm/2.355,'ms')
hold off
h = gca;
set(h,'XLim',[param_min param_max]);
xlabel([pname ' /' punit])
ylabel('{\it\sigma_z} /mm')
title(['s=FW/2.355, d=gauss-fit, ' sprintf(['tol=%7.3g ' punit],Ztol)])
enhance_plot

z_nom = z_bar((n + 1)/2);
dt = 1E9*(z_bar - z_nom)/c;
subplot(224)
%plot_spline(param,dt)
hold on
plot(param,dt,'or')
hold off
h = gca;
set(h,'XLim',[param_min param_max]);
xlabel([pname ' /' punit])
ylabel('\langle\Delta{\itt_f}\rangle /psec')
enhance_plot

Ha = gca;
Hy = get(Ha,'YLabel');
set(Hy,'VerticalAlignment','baseline');
