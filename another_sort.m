clear all;

load('~/Desktop/data/2013/scans/enrg_scan.mat');
load('/Users/sgess/Desktop/data/2013/slims/slim_10794.mat');
specs = zeros(250,200);

for i = 1:length(data.YAG.good_shot)
    LineMinBG = data.YAG.spectra(:,i) - data.YAG.spectra(1,i);
    cent = sum(data.YAG.axis.*LineMinBG)/sum(LineMinBG);
    cs = interp1(data.YAG.axis-cent,LineMinBG,xx,'linear',0);
    specs(i,:) = cs/max(cs);
    cm = repmat(specs(i,:),nsim,1);
    
    res = sum((sy - cm).^2,2);
    [a(i),b(i)] = min(res);
end


%%

% Machine parameter histograms
figure(1);
subplot(2,3,1);
hist(part(b(data.YAG.good_shot)),10);
xlabel('Particles','fontsize',14);

subplot(2,3,2);
hist(sigz(b(data.YAG.good_shot)),10);
xlabel('\sigma_z','fontsize',14);

subplot(2,3,3);
hist(ampl(b(data.YAG.good_shot)),10);
xlabel('Compressor Amplitude','fontsize',14);

subplot(2,3,4);
hist(phas(b(data.YAG.good_shot)),10);
xlabel('Compressor Phase','fontsize',14);

subplot(2,3,5);
hist(deck(b(data.YAG.good_shot)),10);
xlabel('2-10 Phase','fontsize',14);

subplot(2,3,6);
hist(enrg(b(data.YAG.good_shot)),10);
xlabel('2-10 Phase','fontsize',14);

% Find best fit
[min_res,min_ind] = min(a+1000*(~data.YAG.good_shot'));
[max_res,max_ind] = max(a.*(1-(~data.YAG.good_shot))');
figure(2);
plot(xx,specs(max_ind,:),'g',xx,sy(b(max_ind),:),'b','linewidth',2);
xlabel('X (mm)','fontsize',14);
legend('Data','Sim');
title('Worst fit');
figure(3);
hist(I_max(b(data.YAG.good_shot)),30);
title('Peak Current Histogram');
figure(4);
plot(data.YAG.pyro(data.YAG.good_shot),I_max(b(data.YAG.good_shot)),'*');
title('Pyro vs. I');
figure(5);
plot(xx,specs(min_ind,:),'g',xx,sy(b(min_ind),:),'b','linewidth',2);
xlabel('X (mm)','fontsize',14);
legend('Data','Sim');
title('Best fit');
figure(6);
plot(a(data.YAG.good_shot),I_max(b(data.YAG.good_shot)),'*');
title('residual vs. I');