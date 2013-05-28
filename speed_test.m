nOut = 3;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;
param_tcav;

npart = 2e5;
nbin = 128;
it = 100;


nbin = 64;
f_time_64 = test_wake(nbin,npart,nOut,it,'fast_wake');
nbin = 128;
f_time_128 = test_wake(nbin,npart,nOut,it,'fast_wake');
nbin = 256;
f_time_256 = test_wake(nbin,npart,nOut,it,'fast_wake');
nbin = 512;
f_time_512 = test_wake(nbin,npart,nOut,it,'fast_wake');

nbins = [64 128 256 512];
l_time = [median(l_time_64) median(l_time_128) median(l_time_256) median(l_time_512)];
f_time = [median(f_time_64) median(f_time_128) median(f_time_256) median(f_time_512)];

loglog(nbins, l_time,'b:s',nbins, f_time,'r:s','linewidth',3);
axis([50 1000 0.001 10]);
set(gca,'fontsize',16);
xlabel('Wakefield Bins','fontsize',18);
ylabel('Time (s)','fontsize',18);
title('Wake calculation with 200K Particles','fontsize',18);
l1 = legend('long\_wake()','fast\_wake()','location','northwest');




npart = 2e5;
nbin = 128;
it = 100;


npart = 0.5e5;
f_time_50 = test_wake(nbin,npart,nOut,it,'fast_wake');
npart = 1e5;
f_time_100 = test_wake(nbin,npart,nOut,it,'fast_wake');
npart = 2e5;
f_time_200 = test_wake(nbin,npart,nOut,it,'fast_wake');
npart = 4e5;
f_time_400 = test_wake(nbin,npart,nOut,it,'fast_wake');

npart = 0.5e5;
l_time_50 = test_wake(nbin,npart,nOut,it,'long_wake');
npart = 1e5;
l_time_100 = test_wake(nbin,npart,nOut,it,'long_wake');
npart = 2e5;
l_time_200 = test_wake(nbin,npart,nOut,it,'long_wake');
npart = 4e5;
l_time_400 = test_wake(nbin,npart,nOut,it,'long_wake');

nparts = [0.5e5 1e5 2e5 4e5];
l_time = [median(l_time_50) median(l_time_100) median(l_time_200) median(l_time_400)];
f_time = [median(f_time_50) median(f_time_100) median(f_time_200) median(f_time_400)];

loglog(nparts, l_time,'b:s',nparts, f_time,'r:s','linewidth',3);
axis([4e4 1e6 0.001 1]);
set(gca,'fontsize',16);
xlabel('Number of Particles','fontsize',18);
ylabel('Time (s)','fontsize',18);
title('Wake calculation with 128 Bins','fontsize',18);
l1 = legend('long\_wake()','fast\_wake()','location','northwest');