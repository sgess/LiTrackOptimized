clear all;

load('~/Desktop/data/2013/scans/new_scan.mat');
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