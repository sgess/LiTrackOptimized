function [dE,zc,sigz] = long_wake(z,L,Ne,Nbin,fn,pcwake);

%        [dE,zc,sigz] = long_wake(z[,L,Ne,Nbin,fn,pcwake]);
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
%			pcwake:	(Optional, DEF=none) Point-charge wake used instead of file
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			sigz:	rms bunch length (standard deviation) [m]

%=============================================================================

nn   = length(z);
sigz = std(z);

if nn < 100
  disp(' ')
  disp('Probably too few particles in your input "z" array')            
  disp(' ')
end
if nn > 5E6
  disp(' ')
  disp('>5E6 particles in your "z" array - let''s not get carried away now')
  disp(' ')
end
  
if sigz < 1E-6
  disp(' ')
  disp('Bunch length of <1 micron may be too short for this Green''s function')
  disp(' ')
end

if ~exist('L')
  L = 1; 					% default length of S-band linac [m]
end
if ~exist('Ne')
  Ne = 1; 					% default number of e- in bunch      
end
if ~exist('Nbin')
  Nbin = 100;  				% default number simulation particles
end
if ~exist('fn')
  fn = 'slac.dat';  		% default number simulation particles
end  

if ~exist('pcwake')			% if a point-charge wake is not passed in, use file
  cmnd = ['load ' fn];
  eval(cmnd);
  idot = find(fn=='.');
  if isempty(idot)
    error('Point wake file name needs "." in name string')
  end
  cmnd   = ['A = ' fn(1:(idot(1)-1)) ';'];
  eval(cmnd);
else
  [rpc,cpc] = size(pcwake);
  if cpc ~=2
    error('Point-charge wake function needs two columns')
  end
  A = pcwake;				% pcwake is 2-columns, 1st is z [m], 2nd is Wake [V/C/m]
end

zfnvar = A(:,1);			% m
Wfnvar = A(:,2);			% V/C/m
nA     = length(zfnvar);

[N,zc] = hist(z,Nbin-2);		% add zero height bins to both ends...
dzc = mean(diff(zc));			% ...so particles in 1/2-end bins get...
zc = [zc(1)-dzc zc zc(Nbin-2)+dzc];	% ...linearly interpolated...
N = [0 N 0];				% ...properly

maxz_fn = zfnvar(nA-1);			% max Z in wake file (last point usually projected estimate)
if (max(z)-min(z)) > maxz_fn
  disp(' ')
  if ~exist('pcwake')
    disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and ' fn ' is inaccurate there'])
  else
    disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and the RW-wake is inaccurate there'])
  end
  disp(' ')
end

dE  = zeros(Nbin,1);
e   = 1.6022E-19;
scl = -L*(Ne/nn)*e*1E-6;
for j = 1:Nbin				% jth bin is test-bin
  zcj = zc(j);				% save test bin z
  for k = 1:Nbin			% kth bin is field-bin
    zck = zc(k);
    if zck > zcj			% no wake when field-bin behind test-bin
      break
    end
    dz = zcj - zck;			% separation of field & test-bins (>=0)
    [ddz,ii] = min(abs(zfnvar-dz));
    ddz = sign(zfnvar(ii)-dz)*ddz;
    if ddz > 0
      i1 = ii - 1;
      i2 = ii;
      if i1 == 0
        error('Index into zeroth entry of wakefield array - try finer steps in wake file');
      end  
      dz1 = zfnvar(i1);
      dz2 = zfnvar(i2);
      W1  = Wfnvar(i1);
      W2  = Wfnvar(i2);
      W   = W2 - (W2-W1)*ddz/(dz2-dz1);
    elseif ddz < 0
      i1 = ii;
      i2 = ii + 1;
      if i2 > length(zfnvar)
        error('WARN: Index to a point beyond wakefield array - try extending the wake file');
      end
      dz1 = zfnvar(i1);
      dz2 = zfnvar(i2);
      W1  = Wfnvar(i1);
      W2  = Wfnvar(i2);
      W   = W1 - (W2-W1)*ddz/(dz2-dz1);
    else
      W   = Wfnvar(ii)/2;		% use 1/2 wake for self loading of bin
    end
    dE(j) = dE(j) + scl*N(k)*W;		% add each field-bin's wake to each test-bin's energy loss
  end
end