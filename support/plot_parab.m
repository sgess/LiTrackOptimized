function [q,dq,xf,yf] = plot_parab(x,y,dy,xtxt,ytxt,xunt,yunt,noplot);

%PLOT_PARAB
%
%               [q,dq,xf,yf] = plot_parab(x,y,dy,xtxt,ytxt,xunt,yunt,noplot);
%                                   
%               Function to plot X data vs Y data, and fit to a parabola
%               of the form:
%
%                       Y = A*(X-B)**2 + C
%
%               (with error bars if given), Will plot AND return fit
%               coefficients and their errors if output variables are
%               provided.  Otherwise only the plot is generated.
%
%                                
%     INPUTS:   x:      The X-axis data vector (column or row vector)
%               y:      The Y-axis data vector (column or row vector)
%               dy:     (Optional) The Y-axis data vector error 
%                       (or =1 for no errors given; in this case the fit
%                       coefficient errors are rescaled such that 
%                       CHISQ/NDF = 1)
%               xtxt:   (Optional) The text string describing the X-axis data
%               ytxt:   (Optional) The text string describing the Y-axis data
%               xunt:   (Optional) The text string of X-axis units
%               yunt:   (Optional) The text string of Y-axis units
%               noplot: (Optional,DEF=0) If =1, get no graphics plot
%
%     OUTPUTS:  q:      (Optional) If function is written as 
%                       "q = plot_parab(..." then "q" will be a column
%                       vector of fit ciefficients with 
%
%                         q(1) = A,   q(2) = B,   q(3) = C
%
%                       If function is written as "plot_parab(...", then no
%                       output is echoed (plot only)
%               dq:     (Optional) Error on "p" if function is written as 
%                       "[q,dq] = plot_parab(..."
%               xf:     Smooth fitted curve (hor axis) to overlay on data
%               yf:     Smooth fitted curve (ver axis) to overlay on data
%
% P.Emma 02/27/92   Original
% P.Emma 04/10/01   (no-plot option and xf, yf output)

%==========================================================================
                                                       
x  = x(:);                     
y  = y(:);

if ~exist('dy')
  dy = 1;
end

dy = dy(:);

nx = length(x);
ny = length(y);

if length(dy)==1
  bars = 0;
else
  bars = 1;
  if ny~=length(dy)
    error('Y-data and Y-error-data are unequal length vectors')
  end                    
end

if nx~=ny
  error('X-data and Y-data are unequal length vectors')
end

if ny < 3
  error('Not enough data points to fit to this order')
end

if exist('xtxt')==0,
  xtxt= ' ';
end
if exist('ytxt')==0,
  ytxt= ' ';
end

if exist('xunt')==0,
  xunt= ' ';
end
if exist('yunt')==0,
  yunt= ' ';
end

if exist('noplot')==0,
  noplot= 0;
end

Q = [];
for j = 0:2
  Q = [Q x.^j];
end

if bars
  [yfit,dyfit,p,dp,chisq,Cv] = fit(Q,y,dy);
else
  [yfit,dyfit,p,dp,chisq,Cv] = fit(Q,y);
end

if p(3) ~= 0
  A = p(3);
  B = -p(2)/(2*p(3));
  C = p(1) - p(2)^2/(4*p(3));
else
  error('The quadratic coefficient = 0 ...  cannot calc A,B,C')
end

grad_A = [0 0 1]';
grad_B = [0 -1/(2*p(3)) p(2)/(2*p(3)^2)]';
grad_C = [1 -p(2)/(2*p(3)) p(2)^2/(4*p(3)^2)]';

dA = sqrt(grad_A'*Cv*grad_A);
dB = sqrt(grad_B'*Cv*grad_B);
dC = sqrt(grad_C'*Cv*grad_C);

[Qr,Qc] = size(Q);
NDF = Qr-Qc;

difs   = (yfit-y).^2;
difbar = sqrt(mean(difs));
xsig   = std(x);
ysig   = std(y);

lims(1) = min(x);
lims(2) = max(x);
xf = lims(1):(lims(2)-lims(1))/250:lims(2);
Q = [];
for j = 0:2
  Q = [Q xf(:).^j];
end
p = p(:);
yf = Q*p;
xf = xf';

if noplot==0
  if bars
    plot_bars(x,y,dy,'o')
  else
    plot(x,y,'o')
  end
  hold on
  plot(xf(:),yf(:),'-r')
  hold off
  
  title([ytxt ' VS ' xtxt ' Y=A*(X-B)**2+C'])
  xlabel([xtxt ' /' xunt])              
  ylabel([ytxt ' /' yunt])
  
  text(scx(0.60),scy(0.89),[sprintf('A = %8.5g+-%5.3g',A,dA)],'sc')  
  text(scx(0.60),scy(0.86),[sprintf('B = %8.5g+-%5.3g',B,dB)],'sc')  
  text(scx(0.60),scy(0.83),[sprintf('C = %8.5g+-%5.3g',C,dC)],'sc')  
  
  text(scx(0.70),scy(0.02),[sprintf('RMS=%8.5g ',difbar) yunt],'sc')
  text(scx(0.70),scy(0.05),[sprintf('NDF=%5.0f ',NDF)],'sc')
  if bars
    text(scx(0.02),scy(0.02),[sprintf('chisq/N = %5.3f',chisq)],'sc')
  end
  if exist('r')
    text(scx(0.02),scy(0.05),[sprintf('rho     = %6.3f',r)],'sc')
  end
end

if nargout == 1
  q  =  [ A  B  C];
end
if nargout > 1
  q  =  [ A  B  C];
  dq =  [dA dB dC];
end
