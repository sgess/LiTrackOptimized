function [X,Y,Z,dx,dy] = contour_plot(x,y,Nx,Ny,fig);

%	function [X,Y,Z,dx,dy] = contour_plot(x,y,Nx,Ny,fig);
%
%	Function to generate 2D color image from 2D particle ditribution.
%
%	INPUTS:		x:		One of the two coordinates in the array of particles
%	     		y:		The other of the two coordinates in the array of particles
%				Nx:		The number of hor. bins
%				Ny:		The number of ver. bins
%				fig:	[Optional,DEF=1] if fig==1, get plot, else no plot,
%						just output arguments
%	OUTPUTS:	X:		The array of hor. bin-centers
%				Y:		The array of ver. bin-centers
%				Z:		The height of each point (pixel strength or color)
%				dx:		The calibration of x-input-units per bin (e.g., mm/bin)
%				dy:		The calibration of y-input-units per bin (e.g., mm/bin)

%==============================================================================

if ~exist('fig')
  fig = 1;
end

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

dx   = (xmax - xmin)/Nx;
dy   = (ymax - ymin)/Ny;
X    = ((xmin+dx/2):dx:(xmax-dx/2));
Y    = ((ymin+dy/2):dy:(ymax-dy/2));

Z = zeros(Ny,Nx);
for j = 1:Ny
  i  = find(abs(y-Y(j))<dy/2);
  xi = x(i);
  yi = y(i);
  for k = 1:Nx
    ii = find(abs(xi-X(k))<dx/2);
    Z(j,k) = length(ii);
  end 
end

if fig==1
  subplot(221)
  imagesc(X*1E3,Y*1E3,Z);
  axis xy
%  colorbar
  v = axis;
  xlabel('{\itx} (mm)')
  ylabel('{\ity} (mm)')
  enhance_plot('times',20,1);
  subplot(223)
  [Xrms,Xmean,XArea] = rms_calc(X*1E3,sum(Z),0.02);
  stairs(X*1E3,sum(Z),'r-')
  xlim([v(1) v(2)])
  xlabel('{\itx} (mm)')
  ylabel('{\itN}')
  title(['{\it\sigma_x}' sprintf(' = %6.4f mm',Xrms)])
  enhance_plot('times',20,3);
  subplot(222)
  [Yrms,Ymean,YArea] = rms_calc(Y*1E3,sum(Z'),0.02);
  stairs(sum(Z'),Y*1E3,'b-')
  ylim([v(3) v(4)])
  xlabel('{\itN}')
  ylabel('{\ity} (mm)')
  title(['{\it\sigma_y}' sprintf(' = %6.4f mm',Yrms)])
  enhance_plot('times',20,3);
end