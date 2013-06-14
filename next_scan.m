clear all;

file_name = '~/Desktop/wide_scan.mat';
savE = 1;
n_out = 3;
n_interp = 200;

global A;
A = load('slac.dat');

global PARAM;
param_04_16_13;
PARAM.LI20.R16 = 90;

deck_lo = -23.5;
part_lo = 1.90E10;
sigz_lo = 0.0060;
ampl_lo = 0.0395;
phas_lo = 89.0;

deck_hi = -20.5;
part_hi = 2.20E10;
sigz_hi = 0.0080;
ampl_hi = 0.0405;
phas_hi = 91.0;

deck_el = 11;
part_el = 11;
sigz_el = 11;
ampl_el = 11;
phas_el = 11;

decks = linspace(deck_lo,deck_hi,deck_el);
parts = linspace(part_lo,part_hi,part_el);
sigzs = linspace(sigz_lo,sigz_hi,sigz_el);
ampls = linspace(ampl_lo,ampl_hi,ampl_el);
phass = linspace(phas_lo,phas_hi,phas_el);

inds = zeros(deck_el,part_el,sigz_el,ampl_el,phas_el);
nsim = numel(inds);
deck = zeros(nsim,1);
part = zeros(nsim,1);
sigz = zeros(nsim,1);
ampl = zeros(nsim,1);
phas = zeros(nsim,1);

I_max = zeros(nsim,1);
N_par = zeros(nsim,1);

zz = linspace(-0.4,0.4,n_interp);
ee = linspace(-4,4,n_interp);
xx = linspace(-4,4,n_interp);

bl = zeros(nsim,n_interp);
es = zeros(nsim,n_interp);
sy = zeros(nsim,n_interp);


for i = 1:deck_el
    for j = 1:part_el
        for k = 1:sigz_el
            for l = 1:ampl_el
                for m = 1:phas_el
                    
                    ind = sub2ind(size(inds),m,l,k,j,i);
                    disp(['Percent Complete: ' num2str(100*ind/nsim,'%0.2f')]);
                    
                    PARAM.LONE.DECK  = decks(i);
                    PARAM.INIT.NPART = parts(j);
                    PARAM.INIT.SIGZ0 = sigzs(k);
                    PARAM.NRTL.AMPL  = ampls(l);
                    PARAM.NRTL.PHAS  = phass(m);
                    
                    deck(ind) = decks(i);
                    part(ind) = parts(j);
                    sigz(ind) = sigzs(k);
                    ampl(ind) = ampls(l);
                    phas(ind) = phass(m);
                    
                    PARAM.LONE.PHAS = PARAM.LONE.RAMP+PARAM.LONE.DECK;
                    PARAM.LTWO.PHAS = PARAM.LONE.RAMP;
                    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS);
                    PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS);
                    
                    OUT = LiTrackOpt('FACETpar');
                    
                    I_max(ind) = OUT.I.PEAK(n_out);
                    N_par(ind) = OUT.I.PART(n_out);
                    
                    zcent = sum(OUT.Z.AXIS(:,n_out).*OUT.Z.HIST(:,n_out))/sum(OUT.Z.HIST(:,n_out));
                    ecent = sum(OUT.E.AXIS(:,n_out).*OUT.E.HIST(:,n_out))/sum(OUT.E.HIST(:,n_out));
                    xcent = sum(OUT.X.AXIS.*OUT.X.HIST)/sum(OUT.X.HIST);
                    
                    bl(ind,:) = interp1(OUT.Z.AXIS(:,n_out)-zcent,OUT.Z.HIST(:,n_out),zz,'linear',0);
                    es(ind,:) = interp1(OUT.E.AXIS(:,n_out)-ecent,OUT.E.HIST(:,n_out),ee,'linear',0);
                    sy(ind,:) = interp1(OUT.X.AXIS-xcent,OUT.X.HIST,xx,'linear',0);
                    sy(ind,:) = sy(ind,:)/max(sy(ind,:));
                    
                    
                end
            end
        end
    end
end

if savE; save(file_name); end