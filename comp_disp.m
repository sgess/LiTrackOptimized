clear all;
%load('../DATA/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/slim.mat');
load('/Users/sgess/Desktop/data/2013/slims/slim_10794.mat');
%load('../DATA/scans/this_scan.mat');
load('~/Desktop/this_scan.mat');


% dx = xx(2)-xx(1);
% xs = repmat(xx,length(sy),1);
% cent = sum(xs.*sy,2)./sum(sy,2);
% off = round(cent/dx);
% new = zeros(size(sy));
% 
% for i = 1:length(sy)
%     
%     new(i,:) = interp1(xx-cent(i),sy(i,:),xx,'linear',0);
%     new(i,:) = new(i,:)/max(new(i,:));
%     
% end


savE = 1;



LineMinBG = data.YAG.spectra(:,144) - data.YAG.spectra(1,144);
cent = sum(data.YAG.axis.*LineMinBG)/sum(LineMinBG);
cs = interp1(data.YAG.axis-cent,LineMinBG,xx,'linear',0);
cs = cs/max(cs);

cm = repmat(cs,3125,1);

i80=r_16==80;
i85=r_16==85;
i90=r_16==90;
i95=r_16==95;
i100=r_16==100;

res80 = sum((new(i80,:) - cm).^2,2);
[a80,b80] = min(res80);
bi80 = b80+0*3125;

res85 = sum((new(i85,:) - cm).^2,2);
[a85,b85] = min(res85);
bi85 = b85+1*3125;

res90 = sum((new(i90,:) - cm).^2,2);
[a90,b90] = min(res90);
bi90 = b90+2*3125;

res95 = sum((new(i95,:) - cm).^2,2);
[a95,b95] = min(res95);
bi95 = b95+3*3125;

res100 = sum((new(i100,:) - cm).^2,2);
[a100,b100] = min(res100);
bi100 = b100+4*3125;

bInds = [bi80 bi85 bi90 bi95 bi100];

figure(1);
plot(xx,cs,'k','linewidth',3);
hold on;
plot(xx,new(bi80,:),'bs',xx,new(bi85,:),'gv',xx,new(bi90,:),'rd',xx,new(bi95,:),'co',xx,new(bi100,:),'m^');
xlabel('X (mm)','fontsize',14);
legend('YAG','\eta = 80 mm','\eta = 85 mm','\eta = 90 mm','\eta = 95 mm','\eta = 100 mm');
if savE; saveas(gca,'~/Desktop/plots/disp_test/spectra.pdf'); end;


figure(2);
plot(unique(r_16),1000*ampl([bi80 bi85 bi90 bi95 bi100]),'*--','linewidth',2);
xlabel('\eta (mm)','fontsize',14);
ylabel('Compressor Amplitude (MV)','fontsize',14);
if savE; saveas(gca,'~/Desktop/plots/disp_test/nrtl_ampl.pdf'); end;

figure(3);
plot(unique(r_16),deck([bi80 bi85 bi90 bi95 bi100]),'*--','linewidth',2);
xlabel('\eta (mm)','fontsize',14);
ylabel('2-10 Phase (deg)','fontsize',14);
if savE; saveas(gca,'~/Desktop/plots/disp_test/deck_phas.pdf'); end;

figure(4);
plot(unique(r_16),phas([bi80 bi85 bi90 bi95 bi100]),'*--','linewidth',2);
xlabel('\eta (mm)','fontsize',14);
ylabel('Compressor Phase (deg)','fontsize',14);
if savE; saveas(gca,'~/Desktop/plots/disp_test/comp_phas.pdf'); end;

figure(5);
plot(unique(r_16),1000*sigz([bi80 bi85 bi90 bi95 bi100]),'*--','linewidth',2);
xlabel('\eta (mm)','fontsize',14);
ylabel('\sigma_z (mm)','fontsize',14);
if savE; saveas(gca,'~/Desktop/plots/disp_test/sigmaz.pdf'); end;

figure(6);
plot(unique(r_16),part([bi80 bi85 bi90 bi95 bi100]),'*--','linewidth',2);
xlabel('\eta (mm)','fontsize',14);
ylabel('Particles','fontsize',14);
if savE; saveas(gca,'~/Desktop/plots/disp_test/npart.pdf'); end;

figure(7);
plot(1000*zz,bl(bi80,:),'bs',1000*zz,bl(bi85,:),'gv',1000*zz,bl(bi90,:),'rd',1000*zz,bl(bi95,:),'co',1000*zz,bl(bi100,:),'m^');
xlabel('Z (\mum)','fontsize',14);
legend(['\eta = 80 mm, I_{max} = ' num2str(I_max(bi80),'%0.1f') ' (kA)'],...
    ['\eta = 85 mm, I_{max} = ' num2str(I_max(bi85),'%0.1f') ' (kA)'],...
    ['\eta = 90 mm, I_{max} = ' num2str(I_max(bi90),'%0.1f') ' (kA)'],...
    ['\eta = 95 mm, I_{max} = ' num2str(I_max(bi95),'%0.1f') ' (kA)'],...
    ['\eta = 100 mm, I_{max} = ' num2str(I_max(bi100),'%0.1f') ' (kA)']);
if savE; saveas(gca,'~/Desktop/plots/disp_test/profs.pdf'); end;