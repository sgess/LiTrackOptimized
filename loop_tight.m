clear all;

file_name = '~/Desktop/scan.mat';
savE = 1;
n_out = 3;
n_interp = 513;

global A;
A = load('slac.dat');

global PARAM;
param_04_16_13;

disp_lo = 80;
deck_lo = -21.5;
part_lo = 1.90E10;
sigz_lo = 0.0058;
ampl_lo = 0.0395;
phas_lo = 89.9;

disp_hi = 100;
deck_hi = -19.5;
part_hi = 2.30E10;
sigz_hi = 0.0074;
ampl_hi = 0.0405;
phas_hi = 90.5;

disp_el = 5;
deck_el = 5;
part_el = 5;
sigz_el = 5;
ampl_el = 5;
phas_el = 5;

disps = linspace(disp_lo,disp_hi,disp_el);
decks = linspace(deck_lo,deck_hi,deck_el);
parts = linspace(part_lo,part_hi,part_el);
sigzs = linspace(sigz_lo,sigz_hi,sigz_el);
ampls = linspace(ampl_lo,ampl_hi,ampl_el);
phass = linspace(phas_lo,phas_hi,phas_el);

inds = zeros(disp_el,deck_el,part_el,sigz_el,ampl_el,phas_el);
nsim = numel(inds);
r_16 = zeros(nsim,1);
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

bl      = zeros(nsim,n_interp);
es      = zeros(nsim,n_interp);
sy      = zeros(nsim,n_interp);

for i = 1:disp_el
    for j = 1:deck_el
        for k = 1:part_el
            for l = 1:sigz_el
                for m = 1:ampl_el
                    for n = 1:phas_el
                        
                        ind = sub2ind(size(inds),n,l,m,k,j,i);
                        disp(['Percent Complete: ' num2str(ind/nsim,'%0.2f')]);
                        
                        PARAM.LI20.R16   = disps(i);
                        PARAM.LONE.DECK  = decks(j);
                        PARAM.INIT.NPART = parts(k);
                        PARAM.INIT.SIGZ0 = sigzs(l);
                        PARAM.NRTL.AMPL  = ampls(m);
                        PARAM.NRTL.PHAS  = phass(n);
                        
                        r_16(ind) = disps(i);
                        deck(ind) = decks(j);
                        part(ind) = parts(k);
                        sigz(ind) = sigzs(l);
                        ampl(ind) = ampls(m);
                        phas(ind) = phass(n);
                        
                        PARAM.LONE.PHAS = PARAM.LONE.RAMP+PARAM.LONE.DECK;
                        PARAM.LTWO.PHAS = PARAM.LONE.RAMP;
                        PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS);
                        PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS);
                        
                        OUT = LiTrackOpt('FACETpar');
                        
                        I_max(ind) = OUT.I.PEAK(n_out);
                        N_par(ind) = OUT.I.PART(n_out);
                        
                        bl(ind,:) = interp1(OUT.Z.AXIS(:,n_out)-mean(OUT.Z.AXIS(:,n_out)),OUT.Z.HIST(:,n_out),zz,'linear',0);
                        es(ind,:) = interp1(OUT.E.AXIS(:,n_out)-mean(OUT.E.AXIS(:,n_out)),OUT.E.HIST(:,n_out),ee,'linear',0);
                        sy(ind,:) = interp1(OUT.X.AXIS-mean(OUT.X.AXIS),OUT.X.HIST,xx,'linear',0);
                        
                        
                    end
                end
            end
        end
    end
end


