function x = parabolic_dist(N,sig);

%	x = parabolic_dist(N,sig);
%
%	Function to generate a very smooth parabolic distribution:
%
%		f(x) = 3/4/sqrt(5)/sig*(1-1/5*(x/sig).^2),   -sqrt(5)*sig < x < sqrt(5)*sig
%
%	INPUTS:		N:		Number of particles
%				sig:	The standard deviation around the mean (rms)
%	OUTPUTS:	x:		The N particles distributed in a parabola with rms=sig & mean=0

%======================================================================================

Nx = round(N/2);				% find integer value of N/2 (symmetric dist.)
x0 = 4*sqrt(5)*sig/3/N;			% normalization
x1  = zeros(Nx,1);				% form half-sized array initialized to 0
j  = 0;
while 1							% loop until all bins filled (from center out to one edge)
  j = j + 1;
  f  = 1 - 1/5*(x1(j)/sig).^2;
  dx = x0/f;
  if x1(j)+dx > sqrt(5)*sig;	% if beyond edge (= sqrt(5)*sig), stop..
    n = j;						% ..mark total number of particles to here
    break;						% ..jump out of loop
  end
  x1(j+1) = x1(j) + dx;			% add another particle to full array
end
x1 = x1(1:n);					% truncate array to where we stopped
x = [flipud(-x1(2:n)); x1];		% copy left-side to form the right-side of this symm. dist.
nn = length(x);					% total particles in full dist.
if nn==N						% if we came out just right, we're done
  return						% back to main program
end
if nn>N							% if we ended up with a few too many
  x((N+1):end) = [];			% throw away last few
else							% if we ended up with too few
  x = [x' zeros(1,(N-nn))]';	% add a few at x = 0 (not so nice & maybe unnecessary)
end
