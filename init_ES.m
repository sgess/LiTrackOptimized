function [w, dt] = init_ES(nPar)

w0=1000;
w00=5000;

w=zeros(1,nPar);

pr = primes(w0);
lpr = length(pr);

% The w(i) values are w(i) = sqrt(pi), where p1,...,p17 are the 17 primes less than w0. 
% This is a easy quick way to get pretty good "independence" between the
% varius frequencies.
for j2=nPar:-1:1;
    w(j2)=w00*(pr(lpr-5*j2).^0.5);
end

% ES Time Step Size, choose dt small enough so that it takes 20 steps for
% the highest frequency cos(w(17)n dt) to complete one full oscillation
dt=(2*pi)/(8*w(nPar));