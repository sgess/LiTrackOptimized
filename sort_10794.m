clear all;
load('~/Desktop/scan.mat');
load('/Users/sgess/Desktop/data/2013/slims/slim_10794.mat');

ax = data.YAG.axis(121:end);
pix = length(ax);

YAGS = zeros(pix,11^5);

SIGZ_VEC = zeros(11^5,1);
PART_VEC = zeros(11^5,1);
PHAS_VEC = zeros(11^5,1);
AMPL_VEC = zeros(11^5,1);
LIEL_VEC = zeros(11^5,1);
I_VEC    = zeros(11^5,1);

for i = 1:11
    i
    for j = 1:11
        j
        for k = 1:11
            for l = 1:11
                for m = 1:11
                    
                    ind = sub2ind(size(I_max),i,j,k,l,m);
                    
                    SIGZ_VEC(ind) = SIGZ(i);
                    PART_VEC(ind) = PART(j);
                    PHAS_VEC(ind) = PHAS(k);
                    AMPL_VEC(ind) = AMPL(l);
                    LIEL_VEC(ind) = LIEL(m);
                    
                    I_VEC(ind) = I_max(i,j,k,l,m);                    
                    
                    YAGS(:,ind) = interp1(xx(:,i,j,k,l,m),sy(:,i,j,k,l,m),ax,'linear','extrap');
                    
                end
            end
        end
    end
end

%%
dax = ax(2)-ax(1);
off = round(data.YAG.x_cent(1)/dax);
spec = repmat(data.YAG.spectra(121:end,1)/max(data.YAG.spectra(121:end,1)),[1,11^5]);
cent = zeros(pix,11^5);
cent(1:(pix-off+1),:) = spec(off:end,:);
maxYAGS = repmat(max(YAGS),[pix,1]);
nYAGS = YAGS./maxYAGS;

diff = nYAGS-cent;
res = sum(diff.^2);
min(res)

%%

data.YAG.pyro = data.EPICS.BLEN_LI20_3014_BRAW(data.YAG.epics_index);
data.YAG.good_shot = data.YAG.pyro > 5000;

specs = data.YAG.spectra(121:end,data.YAG.good_shot);
maxs = max(specs);
offs = round(data.YAG.x_cent(data.YAG.good_shot)/dax);

ress = zeros(sum(data.YAG.good_shot),1);
inds = zeros(sum(data.YAG.good_shot),1);
for i=1:sum(data.YAG.good_shot)
    i
    spec = repmat(specs(:,i)/maxs(i),[1,11^5]);
    cent = zeros(pix,11^5);
    cent(1:(pix-offs(i)+1),:) = spec(offs(i):end,:);
    
    diff = nYAGS-cent;
    res = sum(diff.^2);
    [a,b] = min(res);
    ress(i) = a;
    inds(i) = b;
    
end

data.YAG.scan_ind = zeros(40,1);
data.YAG.scan_res = zeros(40,1);
data.YAG.scan_ind(data.YAG.good_shot)=inds;
data.YAG.scan_res(data.YAG.good_shot)=ress;