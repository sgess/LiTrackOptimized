function x = gen_dist(N,tail,xsig,xbar,slew);

%	x = gen_dist(N,tail,xsig,xbar,slew);
%
%	Function to generate a generalized distribution by adding
%	multiple offset gaussians together.
%
%    INPUTS:	N:	The total number of particles
%		tail:	The width of the rising and trailing pulse
%			edges (0 < tail < 1) - [e.g.
%			get a steep rise for tail=0.1, 
%			slow for tail=0.5]
%		xsig:	rms of net distribution
%		xbar:	mean of net distribution
%		slew:	slew on pulse (-1 < f < 1) - [get
%			positive slope for f>0]
%    OUTPUTS:	x:	N by 1 random variable
%
%===================================================================

if tail > 1 | tail <= 0
  error('Tail must be <1 and >0')
end

if slew >= 1 | slew <= -1
  error('slew must be < 1 and > -1')
end

n  = 2/tail;
nn = round(n);
x  = [];
fj = 2*((0:(n/(nn-1)):n) - n/2)/n;
Nn = round((N/n)*(1 + slew*fj));
dN = N - sum(Nn);
Nn = Nn + round(dN/nn);
dN = N - sum(Nn);
Nn(round(nn/2)) = Nn(round(nn/2)) + dN;
for j = 1:nn
  xj = tail*randn(Nn(j),1) + 1.5*j*tail;
  x  = [x; xj];
end

x = xsig*x/std(x);
x = x - mean(x) + xbar;
