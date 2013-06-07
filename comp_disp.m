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

figure(1);
plot(xx,cs,'k','linewidth',3);
hold on;
plot(xx,new(bi80,:),'bs',xx,new(bi85,:),'gv',xx,new(bi90,:),'rd',xx,new(bi95,:),'co',xx,new(bi100,:),'m^');
xlabel('X (mm)','fontsize',14);
legend('YAG','\eta = 80 mm','\eta = 85 mm','\eta = 90 mm','\eta = 95 mm','\eta = 100 mm');


