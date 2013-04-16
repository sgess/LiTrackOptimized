figure(7);

twop = 56;

for i = (twop+1):length(params)
    
    par_avg(:,i) = mean(params(:,(i-twop):i),2);
    
end

nmin = 60;
nmax = 2000;

subplot(3,3,1);
plot(par_avg(1,nmin:nmax));
title('NDR bunch length');

subplot(3,3,2);
plot(par_avg(2,nmin:nmax));
title('NDR energy spread');

subplot(3,3,3);
plot(par_avg(3,nmin:nmax));
title('Number of particles');

subplot(3,3,4);
plot(par_avg(4,nmin:nmax));
title('NDR bunch asymmetry');

subplot(3,3,5);
plot(par_avg(5,nmin:nmax));
title('NRTL compressor amplitude');

subplot(3,3,6);
plot(par_avg(6,nmin:nmax));
title('NRTL compressor phase');

subplot(3,3,7);
plot(par_avg(7,nmin:nmax));
title('2-10 Phase');


subplot(3,3,8);
plot(par_avg(8,nmin:nmax));
title('11-20 Phase');

subplot(3,3,9);
plot(par_avg(9,nmin:nmax));
title('Phase Ramp');