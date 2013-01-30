function LT_OUTPUT = LiTrackOpt(fn)
%function LT_OUTPUT = LiTrackOpt(fn)
%
%	Function to do 2-D longitudinal phase space tracking of many
%	electrons through several accelerator sections including bunch
%	compression and the longitudinal wakefield of the accelerating
%  	structures.
%
%	The beam and the various sections are described in your M-file
%	called "*_lit.m" where "*" is a string contained in the above
%	input argument: fn.
%
%  	The initial zpos and dE/E particle coordinates can be 'LiTrack'-
%	generated Matlab gaussian distributions, uniform distributions, or
%	they can be input from a user's ASCII, 2-column file (Z/mm, dE/E/%).
%
%	If no output arguments are provided, this function generates plots
%	(see below).  If at least one output argument is provided, no plots
%	are generated and the longitudinal coordinates of the N particles
%	are returned as described below.
%
%  INPUTS:
%	fn:		    A string which describes the leading characters of the beam-
%			    line M-file with name "*_lit.m" where the fn string is represented
%			    here as "*" (SEE FOR EXAMPLE the file LCLS_LIT.M which is
%			    internall documented and would be run by setting fn='lcls').
%
%  OUTPUTS:     
%	LT_OUPUT:	Struct containg beam distributions and moments
%
%==============================================================================================================

fontsize = 14;


start_time = get_time;
disp(' ')
disp(['LiTrack started:' start_time])
disp(' ')

%load beamline
fnf  = [fn '_lit'];				% build string as file name of beamline-file (BL_file)
fnfm = [fnf '.m'];				% build string as file name of beamline-file (BL_file) including ".m"

Ne = 0;							% must initialize "Ne" before "eval" to make vs. 7.0.1 work (Nov. 2, 2004)
eval(fnf);					    % run BL-file which is just an M-file (*.m)


% default to using all input file points
if ~exist('nsamp','var');       
    nsamp = 1;
end

% random seed used (defaults to 1)
if ~exist('seed','var')         
    seed = 1;
end

% default to no dist. 1%-halo
if ~exist('halo','var')				
    halo = 0;
    halo_pop = 0;					
end

% add LiTrack input 'phase' offset
if exist('z0in','var')
    z0_bar = z0_bar + z0in*1E-3;	
end

% add LiTrack input 'charge' error
if exist('dQ_Q','var')
    if abs(dQ_Q) > 0.8
        error('dQ_Q is > 80%, this seems unreasonable - quitting.')
    end
    Ne1 = Ne*(1+dQ_Q);				
else
    Ne1 = Ne;
end

% number of plot columns (3 for publish-size plots)
if splots == 1
    nn  = 3;					    
    pz  = 3;
    
% number of plot columns (2 for big display-only plots)
else
    nn  = 2;					    
    pz  = 4;
end

% Constants
elec   = 1.6022E-9;				    % Coulombs/(1E10 e-)
cspeed = 2.99792458E8;				% light speed [m/sec]
jf     = 0;					        % start with figure window-1
jc     = 0;                         % start with zero output locations (where cod < 0 | cod == 99)

% Generate Gaussian Dist:
% =======================

inpf = 'gaussian random';			% text to put on plots
%rGen = rng(seed,'twister');                   % set the gaussian random generator seed
d0 = sigd0*randn(Nesim,1) + d0_bar;	% always gaussian dE/E
z0 = asym_gaussian(Nesim,sigz0,z0_bar,...
    asym,cut,tail,halo,halo_pop);	% generate asymmetric gaussian with cuts & tails


% Start doing real calculations:
% ==============================

Ne0    = Ne1;					% save N electrons to get fraction lost
Nesim0 = Nesim;					% save N macro-particles to get fraction lost
iswake = 0;					    % default until a wake is calculated (=0 turns off wake-induced voltage plot)
ecuts  = 0;					    % default to no dE/E cuts shown on plots
zcuts  = 0;					    % default to no Z cuts shown on plots
z      = z0;        			% axial position within bunch w.r.t. mean position [m]
d      = d0;      				% relative energy deviation w.r.t. mean energy [ ]
E      = E0*(1 + d0);			% absolute energy of particles [GeV]
Ebar   = mean(E);				% mean particle energy [GeV]
Ebarcuts = mean(E);				% mean particle energy after later cuts [GeV]
z_bar  = 1E3*mean(z);			% mean z-position [mm]
sigz   = sigz0;      			% rms bunch length (w.r.t. mean) [m]
sigd   = sigd0;      			% rms relative energy spread [ ]

% nb=number of beamline sections in BL-file (e.g. accelerator, compressor, ...)
[nb,ncol] = size(beamline);		            


% Get orders from beamline
for j  = 1:nb					% loop over all beamline sections of BL-file
  cod  = beamline(j,1);			% beamline section code (e.g. 11=accelerator section - from K. Bane convention)
  if (abs(cod)==1)				% do-nothing code (e.g., used to plot input beam)
       %	do nothing - plot only if cod=-1 (see below)
  end
  if (abs(cod)==2)				% dump z,dE/E ASCII file
    fnout = 'LiTrack_zd_output.dat';
    dump_LiTrack_output(1E3*z,100*d,fnout);
    disp(['2-column ASCII output file written: [z(mm) dE/E(%)]: ' fnout])
  end
  if (abs(cod)==11)	|| (abs(cod)==10)	% ACCELERATION SECTION (11 and 10)
    ecuts = 0;					% no dE/E cuts shown on plots
    zcuts = 0;					% no Z cuts shown on plots
    Eacc   = beamline(j,2);		% nominal acc (w/o wake and for phi=crest(=0)) [GeV]
    phi    = beamline(j,3);		% acc phase (crest=0, low-E head at phi < 0) [deg]
    lam    = beamline(j,4);		% RF wavelength [m]
    if lam<=0
      error('RF wavelength cannot be <= 0')
    end
    wakeon = beamline(j,5);		% wakeON=1,2, wakeOFF=0
    Lacc   = beamline(j,6);		% length of acc section (scales wake) [m]
    phir   = phi*pi/180;		% RF phase in radians
    if wakeon					% if wakes calc switched ON...
      iswake = 1;				% turns wake plot on
      nwake_fn = length(wake_fn(:,1));		        % count number of files provided
      if wakeon > nwake_fn
        error('Need multiple wake function file names when "wakeON/OFF" > 1')
      end
      wake_fn1 = strtrim(wake_fn(wakeon,:));			        % select proper wake function depending on wakeon (=1,2,...)
      disp(['Using wake function: ' wake_fn1])		% echo wake file being used
       %[dE_wake,zc_wake] = long_wake(z,Lacc,Ne1,...
       %                    Nbin,wake_fn1);	        % calculate dE_wake in MeV from Z-coordinates, acc length, N-particles, etc.
      [dE_wake,zc_wake] = fast_wake(z,Lacc,Ne1,Nbin,wake_fn1);
      dE_wakes = interp1(zc_wake,dE_wake,z,'*linear');	% inerpolate between wake calc points in bunch to evaluate dE for each e-
      dE_loss = 1E-3*mean(dE_wakes);	% wake loss [GeV]
    else                    			% if wake calc switched OFF...
      iswake = 0;				    	% no wake plot
      dE_wake = zeros(Nbin,1); 			% dE_wake is all zeros
      dE_loss = 0;						% wake loss = 0 without wakes ON
    end
    if abs(cod) == 10					% special case where Eacc is final energy, rather than acc-Voltage
      Eacc = -dE_loss + (Eacc - mean(E))/cos(phir);	% Eacc was final energy, now is acc-volts again [GeV]
    end
    Erf = E + Eacc*cos(phir + 2*pi*z/lam);	% energy of each particle from RF shape alone (no wake yet)
    if wakeon
      E = Erf + dE_wakes*1E-3;		% energy from RF phase and wake added [GeV]
    else
      E = Erf;						% energy from RF phase and NO wake added [GeV]
    end
    Ebar = mean(E);				    % mean particle energy [GeV]
    Ebarcuts = Ebar;				% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;		% relative energy deviation w.r.t. mean energy [ ]
  end
  if abs(cod) == 12                	% zero-out all energy deviations for diagnostics ONLY
    d = zeros(size(d));
    E = mean(E)*ones(size(E));
  end
  if abs(cod) == 13                	% energy feedback with two phase-opposed sections, each of eV0 volts @ crest
    ecuts = 0;					    % no dE/E cuts shown on plots
    zcuts = 0;					    % no Z cuts shown on plots
    Efin   = beamline(j,2);		    % Energy setpoint (goal) [GeV]
    eV0    = beamline(j,3);		    % acc. voltage available at crest for each of two fdbk sections [GeV]
    if eV0==0
      error('Feedback voltage of zero will not correct energy')
    end
    phi1r  = beamline(j,4)*pi/180;  % acc phase of 1st section (crest=0, low-E head at phi < 0) [deg]
    phi2r  = beamline(j,5)*pi/180;	% acc phase of 2nd section (crest=0, low-E head at phi < 0) [deg]
    lam    = beamline(j,6);		    % RF wavelength [m]
    if lam<=0
      error('RF wavelength cannot be <= 0')
    end
    iswake  = 0;				    	% no wake plot
    options = optimset;
    dphi    = fminsearch('fdbk_fun',0,options,phi1r,phi2r,(Efin-mean(E))/eV0);
    fprintf('Energy feedback phase set to %8.3f deg\n',dphi*180/pi)
    En      = mean(E) + eV0*cos(phi1r+dphi) + eV0*cos(phi2r-dphi);
    if abs(En-Efin)/Efin > 1E-4
      fprintf('Energy feedback phase maxed out at %8.3f deg\n',dphi*180/pi)
    end
    E = E  + eV0*cos(phi1r+dphi+2*pi*z/lam) + eV0*cos(phi2r-dphi+2*pi*z/lam);
    Ebar = mean(E);				    % mean particle energy [GeV]
    Ebarcuts = Ebar;				% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;		% relative energy deviation w.r.t. mean energy [ ]
  end                               % end code==13, energy feedback card
  if abs(cod)==15				% resistive-wall wakefield (15)
    r0   = beamline(j,2);		% Beam-pipe radius [m]
    Lng  = beamline(j,3);		% Beam-pipe radius [m]
    sigc = beamline(j,4);		% Surface conductivity [(Ohm-m)^-1]
	tau  = beamline(j,5);		% relaxation time (sec) - if =zero, use DC wake
	rf   = beamline(j,6);		% rf=1: cylindrical chamber,  rf=2: parallel plates chamber
    if r0<=0
      error('Resistive-wall wake cannot be calculated for negative or zero radius')
    end
    if sigc<=0
      error('Resistive-wall wake cannot be calculated for negative or zero conductivity')
    end
    if tau<0
      error('Resistive-wall wake cannot be calculated for negative relaxation time')
    end
    if rf<0 || rf>2
      error('Last paremeter in resistive-wall wake must be 1 or 2, for cylindrical or rectangular chambers')
    end
    iswake = 1;				% turns wake plot on
	Z0 = 120*pi;			% free-space impedance [Ohm]
	s0 = (2*r0^2/(Z0*sigc))^(1/3);
	zmax = 2.01*(max(z) - min(z));
    zpc = 0:(zmax/1000):zmax;
    if tau==0
  	  pcwakeW = rw_wakefield(zpc,r0,s0);		% DC-wake
	else
  	  pcwakeW = rw_wakefield(zpc,r0,s0,tau,rf);	% AC-wake
	end
	pcwake = [zpc(:) -pcwakeW(:)];
    %[dE_wake,zc_wake] = long_wake(z,Lng,Ne1,...
	%							  Nbin,0,pcwake);	        % calculate dE_wake in MeV from Z-coordinates, acc length, N-particles, etc.
    [dE_wake,zc_wake] = fast_wake(z,Lacc,Ne1,Nbin,wake_fn1);                          
    dE_wakes = interp1(zc_wake,dE_wake,z,'*linear');		% inerpolate between wake calc points in bunch to evaluate dE for each e-
    dE_loss = 1E-3*mean(dE_wakes);	% wake loss [GeV]
    E = E + dE_wakes*1E-3;			% energy from RF phase and wake added [GeV]
    Ebar = mean(E);				    % mean particle energy [GeV]
    Ebarcuts = Ebar;				% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;		% relative energy deviation w.r.t. mean energy [ ]
  end
  if abs(cod) == 26                	% USER'S ENERGY CUTS (26) - doesn't change Ebar
    ecuts = 1;					    % show dE/E cuts on plots
    d1 = beamline(j,2);				% minimum dE/E to allow through [ ]
    d2 = beamline(j,3);				% maximum dE/E "   "     "      [ ]
    if d1 >= d2					    % bomb out if max<min (BT-file error)
      error(['Energy cuts (26) must have dE/E_min (col 2) < dE/E_max (col3) in ' fnfm])
    end
    i = find(d>d1 & d<d2);			% bomb out if cuts too tight
    if length(i) < 1
      error('Energy cuts (26) Emin=%7.4f %% and Emax=%7.4f %% threw out all particles',d1*100,d2*100)
    end
    Ni = length(i);				    % count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
    d = d(i);					    % reduce dE/E array inpose cuts
    z = z(i);					    % reduce Z array inpose cuts
    E = E(i);					    % reduce energy array inpose cuts
    Ebarcuts = mean(E);				% mean energy after cuts [GeV]
    disp([sprintf('E-cut (26): %6.3e',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					    % reduce number of simulation particles
  end
  if abs(cod) == 28                	% Notch collimator for M. Hogan
%    ecuts = 1;					    % show dE/E cuts on plots
    d1 = beamline(j,2);				% minimum dE/E for notch-collimator edge [ ]
    d2 = beamline(j,3);				% maximum dE/E for notch-collimator edge [ ]
    if d1 >= d2					    % bomb out if max<min (BT-file error)
        disp('No notch cut');
        %error(['Notch-collimator (28) must have dE/E_min (col 2) < dE/E_max (col3) in ' fnfm])
    else
        i = find(d<d1 | d>d2);			% bomb out if notch too wide
        if length(i) < 1
            error('Notch-collimator (28) Emin=%7.4f %% and Emax=%7.4f %% threw out all particles',d1*100,d2*100)
        end
        Ni = length(i);				    % count particles left after cuts
        Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
        d = d(i);					    % reduce dE/E array inpose cuts
        z = z(i);					    % reduce Z array inpose cuts
        E = E(i);					    % reduce energy array inpose cuts
        Ebarcuts = mean(E);				% mean energy after cuts [GeV]
        disp([sprintf('Notch-collimator (28) cut: %6.3e',100*(1-Ni/Nesim)) '% of bunch'])  
        Nesim = Ni;					    % reduce number of simulation particles
    end
  end
  if abs(cod) == 29                	% USER'S ABSOLUTE ENERGY CUTS (29)
    ecuts = 1;					    % show dE/E cuts on plots
    E1 = beamline(j,2);				% minimum E to allow through [GeV]
    E2 = beamline(j,3);				% maximum E "   "     "      [GeV]
	d1 = E1/Ebar - 1;
	d2 = E2/Ebar - 1;
    if E1 >= E2					    % bomb out if max<min (BT-file error)
      error(['Absolute energy cuts (29) must have E_min (col 2) < E_max (col 3) in ' fnfm])
    end
    i = find(E>E1 & E<E2);			% bomb out if cuts too tight
    if length(i) < 1
      error('Absolute energy cuts (29) Emin=%7.4f GeV and Emax=%7.4f GeV threw out all particles',E1,E2)
    end
    Ni = length(i);				    % count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
    d = d(i);					    % reduce dE/E array inpose cuts
    z = z(i);					    % reduce Z array inpose cuts
    E = E(i);					    % reduce energy array inpose cuts
    Ebarcuts = mean(E);				% mean energy after cuts [GeV]
    Ebar = mean(E);					% mean energy after cuts [GeV]
    disp([sprintf('E-cut (29): %6.3e',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					    % reduce number of simulation particles
  end
  if abs(cod) == 25				    % AUTO-APERTURTE ENERGY WINDOW CUTS (25) - doesn't change Ebar
    ecuts = 1;					    % show dE/E cuts on plots
    iswake = 0;					    % turn off induced voltage plot
    dw = beamline(j,2);				% energy width to allow (max and min set to maximize transmission) [ ]
    dspan = max(d) - min(d);  		% get full span of dE/E
    if dw >= dspan/2				% if E-window is >= 1/2 of full dE/E span...
      Nbin0 = 1;				    % 1-bin needed
    else					        % if E-window < 1/2 of full dE/E span...
      Nbin0 = round(dspan/dw);		% bin dE/E initially (>1)
    end
    Nbin0_max = 250;				% reasonable lower limit on 25-code E-width (dE/E full span/Nbin0_max)
    if Nbin0 > Nbin0_max			% limit N-bins to reasonable scale
      fprintf(['Auto aperture (25) of Ewid=%7.4f%% is just too narrow in ' fnfm '\n'],dw*100)
      fprintf('...estimate %7.4f%% at minimum acceptable for this beam\n',100*dspan/Nbin0_max)
      error('QUITIING')
    end
    if Nbin0 > 1				    % if E-window is < 1/2 dE/E full span...
      [Nd0,D0] = hist(d,Nbin0);		% bin all dE/E to find rough location of max density
      [Nd0max,iNd0max] = max(Nd0);	% find max bin as rough location of most dense population
      dD0 = mean(diff(D0));			% dE/E bin size [ ]
      d10 = D0(iNd0max) - dD0;		% rough minimum dE/E of beam core with dw width
      d20 = D0(iNd0max) + dD0;		% rough maximum dE/E of beam core with dw width
      icore = find(d>d10 & d<d20);	% pointers to core particles
      [Nd,D] = hist(d(icore),Nbin);	% re-bin only core dE/E to find more precise integration limits (d1 and d1+dw)
    else
      [Nd,D] = hist(d,Nbin);		% re-bin all dE/E to find more precise integration limits (d1 and d1+dw)
    end
    dD = mean(diff(D));				% dE/E bin size [ ]
    nj = min([round(dw/dD) Nbin]);	% number of bins which approx. add up to the dw width wanted (max allowed = Nbin)
    if nj < 2					    % bomb if window to narrow
      error(['Auto aperture (25) of Ewid=%7.4f%% is too narrow in ' fnfm],dw*100)
    end
    A = zeros(Nbin-nj+1,1);			% initialize area array
    for jj = 1:(Nbin-nj+1)
      A(jj) = sum(Nd(jj:(jj+nj-1)));% find transmission for each possible integration 1st-limit (const width)
    end
    [Amax,iAmax] = max(A);			% find set of bins with most beam
    if iAmax == Nbin-nj+1			% if most beam is near the high energy edge of dE/E distribution...
      derr = 0;					    % bias window 1/2-bin low so that the high-E dense are is not cut off
    elseif iAmax == 1				% if most beam is near the low energy edge of dE/E distribution...
      derr = dD;				    % bias window 1/2-bin high so that the low-E dense are is not cut off
    else					        % if dense portion of beam is not near edge...
      derr = dD/2;				    % no window bias (1/2-bin shift is necessary for accuracy)
    end 
    d1 = D(iAmax) - derr;			% find first integration limit (low energy cut) for max transmission
    d2 = d1 + dw;				    % find high energy cut
    i = find(d>d1 & d<d2);         	% bomb if window too tight...
    if length(i) < 2
      error('Auto aperture (25), Ewid=%7.4f%% sets Emin=%7.4f%%, Emax=%7.4f%% and throws out all particles',dw*100,d1*100,d2*100)
    end
    Ni = length(i);				    % count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;				% rescale N-particles to reflect cuts
    d = d(i);					    % reduce dE/E array inpose cuts
    z = z(i);					    % reduce Z array inpose cuts
    E = E(i);					    % reduce energy array inpose cuts
    Ebarcuts = mean(E);				% mean energy after cuts [GeV]
    disp([sprintf('E-window-cut (25): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					    % reduce number of simulation particles
  end
  if abs(cod) == 27                         % USER'S constant-dN/N dE/E cuts (27)
%    zcuts = 1;					            % show dE/E-cuts on plots
    dN_N = beamline(j,2);			        % fraction of max-dE/E-amplitude particles to cut [ ]
    no_charge_loss = beamline(j,3);			% if==1, no real charge cut intended, just better binning
    [dsort,idsort] = sort(abs(d-mean(d)));	% sort the absolute value of dE/E values (min to max)
    N1 = round(Nesim*dN_N);			        % throw out last N1 particles
    z(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of zpos
    d(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of dE/E
    E(idsort((Nesim-N1):Nesim)) = [];		% now throw them out of E
    Ni = length(z);				            % count particles left after cuts
    if no_charge_loss==0                    % if real charge cut intended...
      Ne1 = Ne1*Ni/Nesim;				        % ...rescale N-particles to reflect cuts
    end
    Ebarcuts = mean(E);				        % mean energy after cuts [GeV]
    disp([sprintf('Const-dN/N dE/E-cut (27): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					            % reduce number of simulation particles
  end
  if abs(cod) == 35				    		% narrow in on small z-width
    zw1 = beamline(j,2);			        % the full width in z to find peak current in [mm]
    zw2 = beamline(j,3);			        % the full width in z to find peak current in [mm]
	zbar = mean(z);							% z mean [mm]
    dz = z - zbar;							% subtract z mean [mm]
	i  = find(abs(dz)<zw1/2);				% find all particles in this width around the mean
    [Nz,Z] = hist(dz(i),Nbin);				% bin dz in this narrow window to fine peak in current
	[mx,ix] = max(Nz);						% find maximum current in this window
	i  = find(abs(dz-Z(ix))<zw2/2);			% find all particles in this width around the mean
    Ni = length(i);				    		% count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;						% rescale N-particles to reflect cuts
    d = d(i);					    		% reduce dE/E array inpose cuts
    z = z(i);					    		% reduce Z array inpose cuts
    E = E(i);					    		% reduce energy array inpose cuts
    Ebarcuts = mean(E);						% mean energy after cuts [GeV]
    disp([sprintf('z-window-cut (35): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;					    		% reduce number of simulation particles
  end
  if abs(cod) == 36                     	% USER'S Z-CUTS (36)
    zcuts = 1;					            % show Z-cuts on plots
    z1 = beamline(j,2);				        % minimum Z to allow through [m]
    z2 = beamline(j,3);				        % maximum Z "   "     "      [m]
    if z1 >= z2					            % bomb out if max<min (BT-file error)
      error(['Z-cuts (36) must have Z_min (col 2) < Z_max (col3) in ' fnfm])
    end
    i = find(z>z1 & z<z2);					% bomb out if cuts too tight
    if length(i) < 1
      error('Z-cuts (36) Zmin=%7.4f mm and Zmax=%7.4f mm threw out all particles',z1*1E3,z2*1E3)
    end
    Ni = length(i);							% count particles left after cuts
    Ne1 = Ne1*Ni/Nesim;						% rescale N-particles to reflect cuts
    d = d(i);								% reduce dE/E array inpose cuts
    z = z(i);								% reduce Z array inpose cuts
    E = E(i);								% reduce energy array inpose cuts
    Ebarcuts = mean(E);						% mean energy after cuts [GeV]
    disp([sprintf('Z-cut (36): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])  
    Nesim = Ni;								% reduce number of simulation particles
  end
  if abs(cod) == 37                         % USER'S constant-dN/N z cuts (37)
%    zcuts = 1;					            % show Z-cuts on plots
    dN_N = beamline(j,2);			        % fraction of max-z-amplitude particles to cut [ ]
    no_charge_loss = beamline(j,3);			% if==1, no real charge cut intended, just better binning
    [zsort,izsort] = sort(abs(z-mean(z)));	% sort the absolute value of zpos values (min to max)
    N1 = round(Nesim*dN_N);			        % throw out last N1 particles
    z(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of zpos
    d(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of dE/E
    E(izsort((Nesim-N1):Nesim)) = [];		% now throw them out of E
    Ni = length(z);				            % count particles left after cuts
    if no_charge_loss==0                    % if real charge cut intended...
      Ne1 = Ne1*Ni/Nesim;				        % ...rescale N-particles to reflect cuts
    end
    Ebarcuts = mean(E);				        % mean energy after cuts [GeV]
    disp([sprintf('Const-dN/N Z-cut (37): %6.3f',100*(1-Ni/Nesim)) '% of bunch'])
    Nesim = Ni;					            % reduce number of simulation particles
  end
  if abs(cod) == 44             % add a temporal modulation
    mod_amp = beamline(j,2);	% modulation rel. amplitude (typically 0.02 - 0.05)
    mod_lam = beamline(j,3);	% modulation wavelength [m]
    z = z + mod_amp*mod_lam/2/pi*cos(2*pi*z/mod_lam);
    sigz = std(z);				% re-calc bunch length for next possible pass through wake calculations
  end
  if abs(cod) == 45             % add an energy modulation
    mod_amp = beamline(j,2);	% energy modulation relative amplitude (e.g., 0.001 for 0.1%) [ ]
    mod_lam = beamline(j,3);	% energy modulation wavelength [m]
    E = E.*(1 + mod_amp*sin(2*pi*z/mod_lam));	% modulate energy
    Ebar = mean(E);				% mean particle energy [GeV]
    Ebarcuts = Ebar;			% mean energy after cuts - same as Ebar here [GeV]
    d     = (E - Ebar)/Ebar;	% relative energy deviation w.r.t. mean energy [ ]
  end
  if abs(cod) == 6 				% BUNCH COMPRESSION (R56/m, T566/m, E/GeV, U5666/m)
    ecuts = 0;					% no dE/E cuts shown on plots
    zcuts = 0;					% no Z cuts shown on plots
    iswake = 0;					% turn off induced voltage plot
    R56  = beamline(j,2);       % R56 value [m]
    T566 = beamline(j,3);       % T566 value [m] (=-3*R56/2 for non-quad system)
    U5666= beamline(j,5);       % U5666 value [m] (=2*R56 for non-quad system)
    E56  = beamline(j,4);       % Nominal energy of compressor [GeV]
    if E56 < 0.020				% need positive, reasonable nominal R56-energy [GeV]
      fprintf(['WARN: Compressor section (6) of R56=%7.4f m has nominal-energy too small (not ultra-relativistic) in ' fnfm '\n'],R56)
    end
    dd   = (E-E56)/E56;			% relative energy error w.r.t. nominal compressor energy
    z    = R56*dd + T566*dd.^2 + ...
           U5666*dd.^3 + z;		% compress or anti-compress bunch [m]
    sigz = std(z);				% re-calc bunch length for next possible pass through wake calculations
  end
  if abs(cod) == 7 				% BUNCH COMPRESSION CHICANE (R56/m, E/GeV [T566=-1.5*R56, U5666=2*R56])
    ecuts    = 0;					% no dE/E cuts shown on plots
    zcuts    = 0;					% no Z cuts shown on plots
    iswake   = 0;					% turn off induced voltage plot
    R56      = beamline(j,2);       % R56 value [m]
    dR56_R56 = beamline(j,4);       % relative R56 jitter 
    eps      = dR56_R56;            
    if R56>0
      error('R56 for chicane is always <0...  quitting')
    end
    if dR56_R56>0
      disp('Switched-on jitter on R56 in chicane')
    end
    T566 = -1.5*R56;            % T566 value [m]
    U5666=  2.0*R56;            % U5666 value [m]
    E56  = beamline(j,3);       % Nominal energy of compressor [GeV]
    if E56 < 0.020				% need positive, reasonable nominal R56-energy [GeV]
      error(['Chicane section (7) of R56=%7.4f m has nominal-energy too small (not ultra-relativistic) in ' fnfm],R56)
    end
    dd   = (E-E56)/E56;			% relative energy error w.r.t. nominal compressor energy
    z    = R56*(1-eps)*dd + T566*(1-eps)*dd.^2 + ...    % modified by P. Craievich 02/05/06 (relative R56 jitter)
           U5666*(1-eps)*dd.^3 + z ...
           + eps*R56/2;		    % compress or anti-compress bunch [m] 
    sigz = std(z);				% re-calc bunch length for next possible pass through wake calculations
  end
  if abs(cod) == 8 				% octupole ?
	disp('WARN: Octupole element is not reliable yet...')
    ecuts = 0;					% no dE/E cuts shown on plots
    zcuts = 0;					% no Z cuts shown on plots
    iswake = 0;					% turn off induced voltage plot
    K3   = beamline(j,2);       % octupole MAD k-value [m^-4]
    Enom = beamline(j,3);       % Nominal energy in octupole [GeV]
    Leff = beamline(j,4);       % effective magnetic length of octupole [m]
    eta  = beamline(j,5);       % dispersion in octupole [m]
    U5666 = K3*Leff*eta^4/6;
	dd   = (E-Enom)/Enom;		% relative energy error w.r.t. nominal compressor energy
    z    = U5666*dd.^3 + z;		% distort bunch as per octupole, assuming in chicane center [m]
  end
  if abs(cod) == 22			    % INCOHERENT ENERGY SPREAD ADDITION
    ecuts  = 0;				    % no dE/E cuts shown on plots
    zcuts  = 0;				    % no Z cuts shown on plots
    iswake = 0;				    % turn off induced voltage plot
    id_rms = beamline(j,2);     % rms incoherent relative energy spread to be added in quadrature [ ]
	d = d + id_rms*randn(length(d),1);	% incread dE/E by the incoherent addition [ ]
    E = Ebar*(1 + d);			% load energy array [GeV]
  end
  if abs(cod) == 98         % second order energy distribution
      r16 = beamline(j,2);
      t166= beamline(j,3);
      beta= beamline(j,4);
      emit= beamline(j,5);
      gam = Ebar*1000/0.510998928;
      sig = 1000*sqrt(beta*emit/gam);
      XX = r16*d + t166*d.^2 + sig*randn(length(d),1);
  end
  
  %ii = find(z);
  %z_bar = 1E3*mean(z(ii));	% mean z-pos AFTER CUTS [mm]
  
  if cod < 0 || cod == 99    % plot after each negative code point in beamline
                                                % end nargout < 1 (i.e., end plots section)
      jc = jc + 1;                                  % count output locations (where cod < 0 | cod == 99)
      dE_Ej(1:length(d),jc) = d(:);
      zposj(1:length(z),jc) = z(:);
      Ebarj(jc) = Ebar;
      %if nargout > 3				                % if FWHM output parameters wanted...
        zmm     = z*1E3;				            % convert to [mm]
        dpct    = d*100;				            % convert to [%]
		ii = find(zmm);
        [Nz1,Z] = hist(zmm(ii),Nbin);			    % bin the Z-distribution
        HIST_Z(1:Nbin,jc) = Nz1;                           % I'll take that too
        AXIS_Z(1:Nbin,jc) = Z;                             % With the axis
        ZFWmmj(jc) = FWHM(Z,Nz1,0.5);		        % calc. Z-FWHM [mm]
        dZ = mean(diff(Z));    	           	        % Z bin size [mm]
        I  = Nz1*(Ne1/1E10/Nesim)*elec*cspeed/dZ;	% convert N to peak current [kA]
        if gzfit
          [If,q,dq] = gauss_fit(Z,I,1E-3*ones(size(I)),0);
          sigzGj(jc) = q(4);				        % gaussian fit sigma_Z [mm]
          I_pkfj(jc) = max(If);
        end
 		ii = find(dpct);
        [Nd,D]  = hist(dpct(ii),Nbin);			    % bin the dE/E-distribution
        HIST_D(1:Nbin,jc) = Nd;                            % I'll take that too
        AXIS_D(1:Nbin,jc) = D;                             % With the axis
        dFWpctj(jc) = FWHM(D,Nd,0.5);		        % calc. dE/E-FWHM [%]
        %z_barj(jc) = z_bar;
        Ebarcutsj(jc) = Ebarcuts;
        I_pk1 = max(I);
        I_pkj(jc) = I_pk1;
        if gdfit
          [yf,q,dq] = gauss_fit(D,Nd,1E-3*ones(size(Nd)),0);
          sigEGj(jc) = q(4);				        % gaussian fit sigma_dE/E0 [%]
        end  
        %fcutj(jc)    = 1 - Ne1/Ne0;					% fraction of particles cut [ ]
        fcutj(jc) = Nesim;
      %end
  end                                               % end cod < 0 | cod == 99 stuff
  
  if abs(cod) == 99
    break
  end
end                                                 % end loop over all beamline section

if nargout == 1
    LT_OUTPUT.Z.DIST = zposj;
    LT_OUTPUT.Z.HIST = HIST_Z;
    LT_OUTPUT.Z.AXIS = AXIS_Z;
    %LT_OUTPUT.Z.AVG  = z_barj;
    LT_OUTPUT.Z.FWHM = ZFWmmj;
    if gzfit
        LT_OUTPUT.Z.SIG  = sigzGj;
        LT_OUTPUT.I.SIG  = I_pkfj;
    end

    LT_OUTPUT.E.DIST = dE_Ej;
    LT_OUTPUT.E.HIST = HIST_D;
    LT_OUTPUT.E.AXIS = AXIS_D;
    LT_OUTPUT.E.AVG  = Ebarj;
    LT_OUTPUT.E.FWHM = dFWpctj;
    if gdfit
        LT_OUTPUT.E.SIG  = sigEGj;
    end

    LT_OUTPUT.I.PART = fcutj;
    LT_OUTPUT.I.PEAK = I_pkj;
    
    if exist('XX','var')
        LT_OUTPUT.X.DIST = XX;
        [Nx,X]  = hist(XX,Nbin);
        LT_OUTPUT.X.HIST = Nx;
        LT_OUTPUT.X.AXIS = X;
    end
end

end_time = get_time;
disp(' ')
disp(['LiTrack ended:' end_time])
