clear all;
load('../DATA/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/slim.mat');
load('../DATA/scans/this_scan.mat');

max_sy = max(sy,[],2);
maxes = repmat(max_sy,1,513);
specs = sy./maxes;


LineMinBG = data.YAG.spectra(:,144) - data.YAG.spectra(1,144);
cent = sum(data.YAG.axis.*LineMinBG)/sum(LineMinBG);
cs = interp1(data.YAG.axis-cent,LineMinBG,xx,'linear',0);
cs = cs/max(cs);

cm = repmat(cs,nsim,1);

all_res = sum((specs - cm).^2,2);

[a,b] = min(all_res);

plot(xx,cs,xx,specs(b,:));



