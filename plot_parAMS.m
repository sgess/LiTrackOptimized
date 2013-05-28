figure(6);

twop = 1;

for i = 1:length(params)
    
    par_avg(:,i) = params(:,i);
    
end

subplot(7,2,1);
plot(sim_time,par_avg(1,:));
title('NDR bunch length');
datetick;

subplot(7,2,2);
plot(par_avg(2,nmin:nmax));
title('NDR energy spread');

subplot(7,2,3);
plot(par_avg(3,nmin:nmax));
title('Number of particles');

subplot(7,2,4);
plot(par_avg(4,nmin:nmax));
title('NDR bunch asymmetry');

subplot(7,2,5);
plot(par_avg(5,nmin:nmax));
title('NRTL compressor amplitude');

subplot(7,2,6);
plot(par_avg(6,nmin:nmax));
title('NRTL compressor phase');

subplot(7,2,7);
plot(par_avg(7,nmin:nmax));
title('NRTL R56');

subplot(7,2,8);
plot(par_avg(8,nmin:nmax));
title('NRTL T566');

subplot(7,2,9);
plot(par_avg(9,nmin:nmax));
title('2-10 Phase');


subplot(7,2,10);
plot(par_avg(10,nmin:nmax));
title('11-20 Phase');

subplot(7,2,11);
plot(par_avg(11,nmin:nmax));
title('Phase Ramp');

subplot(7,2,12);
plot(par_avg(12,nmin:nmax));
title('Beta at sYAG');

subplot(7,2,13);
plot(par_avg(13,nmin:nmax));
title('Dispersion at sYAG');

subplot(7,2,14);
plot(par_avg(14,nmin:nmax));
title('T166 at sYAG');