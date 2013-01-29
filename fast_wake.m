function [dE,zc,sigz] = fast_wake(z,L,Ne,Nbin,fn)
%===============================================================================
%        [dE,zc,sigz] = fast_wake(z[,L,Ne,Nbin,fn]);
%
%	Function to return the wakefield induced energy profile vs. z for
%	a set of given axial coordinates "z".
%
%  INPUTS:	z:		The internal axial coordinates, within the bunch, of
%					each electron with respect to any fixed point [m]
%			L:		(Optional, DEF=1 m) The length of the linac [m]
%			Ne:		(Optional, DEF=1  ) The number of electrons in the bunch
%			Nbin:   (Optional, DEF=100) The number of bins to use
%			fn:		(Optional, DEF='slac.dat') File name containing longitudinal
%					point wake (DEF='slac.dat')
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			sigz:	rms bunch length (standard deviation) [m]
%===============================================================================

% Number of simulated particles
nn = length(z);

% RMS longitudinal spread
sigz = std(z);

% Wakefile: A(:,1) = Z (m), A(:,2) = W(Z) (V/C/m)
%fn = strtrim(fn);
%A  = load(fn);
global A;
nA = length(A);

% Histogram particles into Nbins with zeros at ends for interpolation
[N,zc] = hist(z,Nbin-2);
% Bin spacing (m)
dzc = zc(2)-zc(1);
% Zero padded axis
zc = [zc(1)-dzc zc zc(Nbin-2)+dzc];
zc = zc';
% Zero padded histogram
N = [0 N 0];
N = N';

% max Z in wake file (last point usually projected estimate)
maxz_fn = A(nA-1);
if (zc(Nbin)-zc(1)) > maxz_fn
    disp(' ')
    disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and ' fn ' is inaccurate there'])
    disp(' ')
end

% delta vector
dE = zeros(Nbin,1);
% SI electric charge
e = 1.6022E-19;
% scale variable (mC/V ?)
scl = -L*(Ne/nn)*e*1E-6;

% Bin separation vector
dzi = dzc*((1:Nbin)' - 1);

% Interpolated wake vector
Wf = interp1(A(:,1),A(:,2),dzi);

% Self wake bin
Wf(1) = Wf(1)/2;

% Sum delta due to wake from leading bins
for j =1:Nbin
    dE(j) = sum(scl*N(j:-1:1).*Wf(1:j));
end
