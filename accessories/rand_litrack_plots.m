nhist = 50;
wd  = 0.36;
ht  = 0.34;
lf1 = 0.08;
lf2 = 0.60;
bt1 = 0.60;
bt2 = 0.10;
lim = 250;
figure(1)
clf;
subplot;
subplot('position',[lf1 bt1 wd ht])
[N,X]=hist(dd,nhist);
bar(X,N,'c')
%axis([-0.4 0.4 0 lim])
%xlim([-0.3 0.3])
ver_line(0,'r--')
ver_line(-0.1,'g:')
ver_line(+0.1,'g:')
xlabel('\langle\Delta{\itE}/{\itE}_0\rangle /%')
ylabel('\itN')
title([sprintf('MEAN=%7.3f %%, ',mean(dd)) sprintf('RMS=%7.3f %%',std(dd))])
enhance_plot

dEw = sigE;
subplot('position',[lf2 bt1 wd ht])
[N,X]=hist(dEw,nhist);
bar(X,N,'c')
%axis([0 0.04 0 lim])
%xlim([0 2.0])
ver_line(sigE0,'r--')
xlabel('(\Delta{\itE}/{\itE}_0)_{rms} /%')
ylabel('\itN')
title([sprintf('MEAN=%7.4f %%, ',mean(dEw)) sprintf('RMS=%7.4f %%',std(dEw))])
enhance_plot

%zw = sigz;
zw = sigzG(:,end);
[mx,imx] = max(zw);
%zw(imx) = [];
subplot('position',[lf1 bt2 wd ht])
[N,X]=hist(zw*1E3,nhist);
bar(X,N,'c')
%axis([0 35 0 lim])
%xlim([0 35])
ver_line(sigz0(end)*1E3,'r--')
ver_line(sigz0(end)*1E3*0.9,'g:')
ver_line(sigz0(end)*1E3*1.1,'g:')
xlabel('{\it\sigma_z} /\mum')
ylabel('\itN')
title([sprintf('MEAN=%5.2f ',mean(zw*1E3)) '\mum, ' sprintf('RMS=%5.2f ',std(zw*1E3)) '\mum'])
enhance_plot

ts = 1E9*(z_bar(:,end)-mean(z_bar(:,end)))/c;     % timing error [ps]
subplot('position',[lf2 bt2 wd ht])
[N,X]=hist(ts,nhist);
bar(X,N,'c')
%[q,dq] = gauss_plot(X,N,1,0);
%axis([-0.5 0.5 0 lim])
%xlim([-1 1])
ver_line(0,'r--')
xlabel('\langle\Delta{\itt_f}\rangle /ps')
ylabel('\itN')
title([sprintf('MEAN=%5.3f ps, ',mean(ts)) sprintf('RMS=%5.3f ps',std(ts))])
enhance_plot

figure(2)
clf;
subplot;
subplot('position',[lf1 bt1 wd ht])
[N,X]=hist(Ipk/1e3,nhist);
bar(X,N,'c')
%axis([0 5 0 lim])
%xlim([0 35])
ver_line(Ipk0/1e3,'r--')
ver_line(Ipk0*0.9/1e3,'g:')
ver_line(Ipk0*1.1/1e3,'g:')
xlabel('{\itI_{pk}} /kA')
ylabel('\itN')
title([sprintf('MEAN=%7.3f kA, ',mean(Ipk)/1e3) sprintf('RMS=%7.3f kA',std(Ipk)/1e3)])
enhance_plot

subplot('position',[lf2 bt1 wd ht])
[mx,imx] = max(ZFWmm);
ZFWmmc = ZFWmm(:,end);
%ZFWmmc(imx) = [];
[N,X]=hist(1E12*ZFWmmc/c,nhist);
bar(X,N,'c')
%axis([0 350 0 lim])
%xlim([0 250])
ver_line(1E12*ZFWmm0(end)/c,'r--')
ver_line(0.9*1E12*ZFWmm0(end)/c,'g:')
ver_line(1.1*1E12*ZFWmm0(end)/c,'g:')
xlabel('\Delta{\itt}_{FWHM} /fs')
ylabel('\itN')
title([sprintf('MEAN=%5.1f fs, ',mean(1E12*ZFWmmc/c)) sprintf('RMS=%5.1f fs',std(1E12*ZFWmmc/c))])
enhance_plot
