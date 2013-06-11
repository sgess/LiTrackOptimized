load('../DATA/scans/scan.mat');

n_interp = 129;
nsim = numel(I_max);
new  = zeros(nsim,n_interp);

cent = squeeze(sum(xx.*sy)./sum(sy));

x2 = linspace(-4,4,n_interp)';
%dx = xx(2)-xx(1);
% xs = repmat(xx,1,size(I_max));
% cent = squeeze(sum(xs.*sy,2)./sum(sy,2));
% off = round(cent/dx);

% 
% for i = 1:length(sy)
%     
%     new(i,:) = interp1(xx-cent(i),sy(i,:),xx,'linear',0);
%     new(i,:) = new(i,:)/max(new(i,:));
%     
% end

for i = 1:11
    disp(i);
    for j = 1:11
        disp(j);
        for k = 1:11
            for l = 1:11
                for m = 1:11
                    
                    ind = sub2ind(size(I_max),l,m,k,j,i);
                    new(ind,:) = interp1(xx(:,i,j,k,l,m)-cent(i,j,k,l,m),sy(:,i,j,k,l,m),x2,'linear',0);
                end
            end
        end
    end
end

%%
load('../DATA/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/slim.mat');

blarge = data.EPICS.BLEN_LI20_3014_BRAW(data.YAG.epics_index) > 1e4;
ress = 1000*ones(length(blarge),1);
inds = zeros(length(blarge),1);
specs = zeros(200,250);
%mn = max(new,[],2);
%ms = repmat(mn,1,129);
%newer = new./ms;

for i=1:length(blarge)
    disp(i);
    if blarge(i)

        LineMinBG = data.YAG.spectra(:,i) - data.YAG.spectra(1,i);
        cent2 = sum(data.YAG.axis.*LineMinBG)/sum(LineMinBG);
        cs = interp1(data.YAG.axis-cent2,LineMinBG,xx,'linear',0);
        specs(:,i) = cs/max(cs);
        cm = repmat(cs,nsim,1);

        res = sum((sy - cm).^2,2);
        [ress(i),inds(i)] = min(res);
        
    end
end


%%

specs = zeros(129,250);

for i=1:length(blarge)
    disp(i);
    if blarge(i)

        LineMinBG = data.YAG.spectra(:,i) - data.YAG.spectra(1,i);
        cent2 = sum(data.YAG.axis.*LineMinBG)/sum(LineMinBG);
        cs = interp1(data.YAG.axis-cent2,LineMinBG,x2,'linear',0);
        specs(:,i) = cs/max(cs);
        %cm = repmat(cs,1,nsim)';

        %res = sum((newer - cm).^2,2);
        %[ress(i),inds(i)] = min(res);
        
    end
end

%%

%for i=1:length(blarge); if blarge(i); plot(x2,specs(:,i),x2,newer(inds(i),:)); pause; end; end
for i=1:length(blarge); if blarge(i); plot(xx,specs(:,i),xx,sy(inds(i),:)); pause; end; end

%%
Is = zeros(length(blarge),1);

for i=1:length(blarge); if blarge(i); Is(i) = I_max(ind2sub(size(I_max),inds(i))); end; end