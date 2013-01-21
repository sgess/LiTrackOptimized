        function [sigx_tol,mean_y] = jitter_tol(x,y,sigy_tol,gu,noplot);

%       [sigx_tol,mean_y] = jitter_tol(x,y,sigy_tol,gu,noplot);
%
%       Function to calculate the rms jitter tolerance on parameter "x" for
%       a desired rms tolerance "sigy_tol" on parameter "y", when "y"
%       is quadratically dependent on "x" as:
%
%               y = a + bx + cx^2 ,
%
%       for zero-mean "x" jitter which is Gaussian or Uniform
%       ("gu" = 1 or 2).
%
%   INPUTS:     x:              The raw x data (at least 3 points needed)
%               y:              The raw y data ( "    "  "   "      "   )
%               sigy_tol:       The rms tolerance on y (in y's input units)
%               gu:             (Optional, DEF=1) gu=1 is for Gaussian
%                               x-jitter, gu=2 is for uniform x-jitter
%               noplot:         (Optional, DEF=0) If noplot==1, get no plot
%
%   OUTPUTS:    sigx_tol:       The tolerance on parameter x (in x's input
%                               units)
%               mean_y:         The new mean value of y given rms x-jitter
%                               of sigx_tol

%========================================================================

if ~exist('gu')
  gu = 1;
end

if ~exist('noplot')
  noplot = 0;
end

if length(x) < 3
  error('Need at least 3 x-points')
end

if length(y) < 3
  error('Need at least 3 y-points')
end

if length(y) ~= length(x)
  error('Need equal number of x and y points')
end

[q,dq] = plot_polyfit(x,y,1,2,'X','Y',' ',' ',noplot);

a = q(1);
b = q(2);
c = q(3);

if gu == 1
  tt = -(b^2) + sqrt(b^4 + 8*c^2*sigy_tol^2);
  if abs(tt)>1E-10
    sigx_tolp = sqrt(  tt/(4*c^2)  );
    sigx_tolm = sqrt( -tt/(4*c^2)  );
  else
    if b==0
      sigx_tolp = 999;
      sigx_tolm = 0;
    else
      sigx_tolp = sigy_tol/b;
      sigx_tolm = 0;
	end
  end
else
  sigx_tolp = sqrt(  (-(b^2) + sqrt(b^4 + (16/5)*c^2*sigy_tol^2))/((8/5)*c^2)  );
  if sigx_tolp==0
    sigx_tolp = sigy_tol/b;
  end
  sigx_tolm = sqrt(  (-(b^2) - sqrt(b^4 + (16/5)*c^2*sigy_tol^2))/((8/5)*c^2)  );
  if sigx_tolm==0
    sigx_tolm = sigy_tol/b;
  end
end

if isreal(sigx_tolp)
  sigx_tol = sigx_tolp;
else
  sigx_tol = sigx_tolm;
end

mean_y = a + c*sigx_tol^2;
if noplot==0
  hor_line(mean_y+sigy_tol)
  hor_line(mean_y-sigy_tol)
  ver_line( sigx_tol)
  ver_line(-sigx_tol)
  enhance_plot
end