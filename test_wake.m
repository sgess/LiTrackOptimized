function time = test_wake(nbin,npart,nOut,it,wake_fun)

global PARAM;
PARAM.SIMU.BIN = nbin;
PARAM.INIT.NESIM = npart;
Lacc = 1;

OUT = LiTrackOpt('FACETpar');
Z = OUT.Z.DIST(:,nOut);
P = PARAM.INIT.NPART*OUT.I.PART(nOut)/PARAM.INIT.NESIM;

time = zeros(1,it);

if strcmp(wake_fun,'fast_wake')
    for i=1:it
        
        tic;
        fast_wake(Z,Lacc,P,PARAM.SIMU.BIN);
        time(i) = toc;
        
    end
end
    
if strcmp(wake_fun,'long_wake')
    for i=1:it
        
        tic;
        long_wake(Z,Lacc,P,PARAM.SIMU.BIN);
        time(i) = toc;
        
    end
end