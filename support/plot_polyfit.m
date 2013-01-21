function [q,dq,xf,yf] = plot_polyfit(x,y,dy,n,xtxt,ytxt,xunt,yunt,noplot,axisv,ll,res);
                                         
%       [q,dq,xf,yf] = plot_polyfit(x,y,dy,n,xtxt,ytxt[,xunt,yunt,noplot,axisv,ll,res]);
%                                   
%       Function to plot x data vs y data, and fit to a polynomial
%       of order "n" (with error bars if given), and order skipping
%       available (see below).  Will plot AND return fit coefficients
%       and their errors if output variables are provided.  Otherwise
%       only the plot is generated.
%
%                                
%     INPUTS:   x:      The X-axis data vector (column or row vector)
%               y:      The Y-axis data vector (column or row vector)
%               dy:     The Y-axis data vector error (or =1 for no errors
%                       given; in this case the fit coefficient errors
%                       are rescaled such that CHISQ/NDF = 1)
%               n:      Polynomial order to fit:
%                       e.g.  n=4 implies 0th, 1st, 2nd, 3rd, & 4th order.
%                       while n=[1 0 1 1 0] implies only 0th, 2nd, & 3rd order
%                       fit
%               xtxt:   (Optional) The text string describing the X-axis data
%               ytxt:   (Optional) The text string describing the Y-axis data
%               xunt:   (Optional) The text string of X-axis units
%               yunt:   (Optional) The text string of Y-axis units
%		        noplot:	(Optional,DEF=0) If noplot=1, get no graphical plot
%               axisv:  (Optional) 4-element vector of forced plot scale
%                       [xmin xmax ymin ymax], or not used if 1-element
%               ll:     (Optional) If "ll"==1, fit coefficients will be shown
%                       at lower left of plot, rather than at upper right, if
%						ll=2 get no coefficients on plot
%               res:    (Optional) If "r"==1, plot residuals of fit
%     OUTPUTS:  q:      (Optional) If function is written as "q = polyfit(..."
%                       then "q" will be a column vector of fit coefficients
%                       with row-1 as lowest order...
%                       If function is written as "plot_polyfit(...", then no
%                       output is echoed (plot only)
%               dq:     (Optional) Error on "p" if function is written as 
%                       "[q,dq] = polyfit(..."
%               xf:     (Optional) More densely sampled abscissa if written as
%                       "[q,dq,xf,yf] = polyfit(..."
%               yf:     (Optional) Fitted function if function is written as 
%                       "[q,dq,xf,yf] = polyfit(..."
%
% P.Emma 11/13/91

%==========================================================================

x  = x(:);                     
y  = y(:);
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

if ny < n+1
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
  noplot = 0;
end
if exist('ll')==0,
  ll = 0;
end
if exist('res')==0,
  res = 0;
end

nn = length(n);
if nn == 1
  mm = n+1;
  m  = 0:1:mm;
else
  i = find(n);
  m = i - 1;
  mm = length(m);
end

Q = [];

f = mean(abs(x));
x1 = x/f;

for j = 1:mm
  Q = [Q x1.^m(j)];
end

if bars
  [yfit,dyfit,p1,dp1,chisq] = fit(Q,y,dy);
else                        
  [yfit,dyfit,p1,dp1] = fit(Q,y);
end

for j = 1:mm
  p(j)  = p1(j)/(f^m(j));
  dp(j) = dp1(j)/(f^m(j));
end

if mm == 2
  if (m(1)==0) & (m(2)==1)
    r = p(2)*std(x)/std(y);
  end
end

[Qr,Qc] = size(Q);
NDF = Qr-Qc;

difs   = (yfit-y).^2;
difbar = sqrt(mean(difs));
xsig   = std(x);
ysig   = std(y);

if res == 1
  y = y-yfit;
  ytxt = [ytxt ' (FIT RESIDUALS)'];
end

if noplot==0
  if bars
    plot_bars(x,y,dy,'or')
  else
    plot(x,y,'or')
  end
  if exist('axisv')==1,
    if length(axisv)==4
      axis(axisv)
    end
  end
end  

if res == 0
  if noplot==0
    lims = get(gca,'XLim');
  else
    lims = [min(x) max(x)];
  end
  xp = lims(1):(lims(2)-lims(1))/250:lims(2);
  Q = [];
  for j = 1:mm
    Q = [Q xp(:).^m(j)];
  end
  p = p(:);
  yp = Q*p;
  if noplot==0
    hold on
    plot(xp(:),yp(:),'-b')
    hold off
  end  
end

if noplot==0
  title([ytxt ' VS ' xtxt])
  xlabel([xtxt ' /' xunt])              
  ylabel([ytxt ' /' yunt])

  if ll == 1
    xpn = 0.20;
    ypn = 0.41;
  else
    xpn = 0.10;
    ypn = 1.00;
  end

  if ll<2
    for j = 1:mm
      text(scx(xpn),scy(ypn-j*0.050), ...
        [sprintf('p%1.0f=%7.4g+-%6.4g',m(j),p(j),dp(j))])  
    end
    text(scx(0.10),scy(0.030),[sprintf('RMS=%5.3g ',difbar) yunt])
    text(scx(0.10),scy(0.080),[sprintf('NDF=%5.0f ',NDF)])
    if bars
      text(scx(0.02),scy(0.030),[sprintf('chisq/N = %5.3f',chisq)])
    end
    if exist('r')
      text(scx(0.02),scy(0.080),[sprintf('rho     = %6.3f',r)])
    end
  end
end

if nargout == 1
  q  =  p;
end
if nargout == 2
  q  =  p;
  dq = dp;
end
if nargout >= 3
  q  =  p;
  dq = dp;
  xf = xp;
  yf = yp;
end