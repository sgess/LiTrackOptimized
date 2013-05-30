clear all;

file_name = '~/Desktop/scan.mat';
savE = 1;
scan = 1;
test = 0;
global PARAM;
param_04_16_13;
PARAM.LI20.R16   = 85;
PARAM.LI20.BETA  = 3.0;
PARAM.LI20.T166  = 0;
    
global A;
A = load('slac.dat');

if scan == 1
    
    n_out = 3;
    
    %sigz0 range
    sz_lo = 0.006;
    sz_hi = 0.007;
    
    % particle range
    n_lo = 1.9e10;
    n_hi = 2.1e10;
    
    % comp phase range
    c_lo = 89.5;
    c_hi = 90.5;
    
    % comp phase range
    a_lo = 0.0400;
    a_hi = 0.0410;
    
    %2-10 range
    li_lo = -23;
    li_hi = -21;
    
    
    
    
    
    % number of sample points
    sz_el = 11;
    n_el  = 11;
    c_el  = 11;
    a_el  = 11;
    li_el = 11;
    
    if test
        sz_el = 1;
        n_el  = 1;
        c_el  = 1;
        a_el  = 1;
        li_el = 1;
    end
    
    % phase and ampl vec
    SIGZ = linspace(sz_lo,sz_hi,sz_el);
    PART = linspace(n_lo,n_hi,n_el);
    PHAS = linspace(c_lo,c_hi,c_el);
    AMPL = linspace(a_lo,a_hi,a_el);
    LIEL = linspace(li_lo,li_hi,li_el);
    
    I_max   = zeros(sz_el,n_el,c_el,a_el,li_el);
    
    N_par   = zeros(sz_el,n_el,c_el,a_el,li_el);
    
    bl      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
    es      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
    sy      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
    
    zz      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
    ee      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
    xx      = zeros(PARAM.SIMU.BIN,sz_el,n_el,c_el,a_el,li_el);
     
    
    for i = 1:sz_el
        for j = 1:n_el
            for k = 1:c_el
                for l = 1:a_el
                    for m = 1:li_el
            

                        PARAM.INIT.SIGZ0 = SIGZ(i);
                        PARAM.INIT.NPART = PART(j);
                        PARAM.NRTL.PHAS  = PHAS(k);
                        PARAM.NRTL.AMPL  = AMPL(l);
                        PARAM.LONE.PHAS  = LIEL(m);
                        
                        
                        PARAM.LONE.GAIN  = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS);
                        
                        OUT = LiTrackOpt('FACETpar');
                        
                        I_max(i,j,k,l,m) = OUT.I.PEAK(n_out);
                        N_par(i,j,k,l,m) = OUT.I.PART(n_out);
                        
                        bl(:,i,j,k,l,m) = OUT.Z.HIST(:,n_out);
                        es(:,i,j,k,l,m) = OUT.E.HIST(:,n_out);
                        zz(:,i,j,k,l,m) = OUT.Z.AXIS(:,n_out);
                        ee(:,i,j,k,l,m) = OUT.E.AXIS(:,n_out);
                        
                        sy(:,i,j,k,l,m) = OUT.X.HIST;
                        xx(:,i,j,k,l,m) = OUT.X.AXIS;
                    end
                end
            end
        end
    end
    
    if savE
        save(file_name,'PARAM','SIGZ','PART','PHAS','AMPL','LIEL',...
            'I_max','N_par','bl','es','zz','ee','sy','xx');
    end
end