  global PARAM; % Initial beam and machine parameters  
  
  % FACET energy set points. We are pretty sure about these, I think. . .
  E0        = PARAM.ENRG.E0;    % GeV ... initial energy
  E1        = PARAM.ENRG.E1;    % GeV ... energy at LBCC
  E2        = PARAM.ENRG.E2;    % GeV ... energy at FACET
  
  % NRTL compressor klystron and R56 #s
  NRTL_ampl = PARAM.NRTL.AMPL;  % AMPL DR13 11 VDES
  NRTL_phas = PARAM.NRTL.PHAS;  % on the zero-crossing
  NRTL_leff = PARAM.NRTL.LEFF;  % cavity length
  NRTL_R56  = PARAM.NRTL.R56;   % This is design val
  NRTL_T566 = PARAM.NRTL.T566;  % Design val?
  NRTL_ELO  = PARAM.NRTL.ELO;   % NRTL low energy cut
  NRTL_EHI  = PARAM.NRTL.EHI;   % NRTL high energy cut
  
  % Phase and length of 02-10
  LONE_leff = PARAM.LONE.LEFF;  % Length of LI02-LI10 (m)
  LONE_phas = PARAM.LONE.PHAS;  % Chirp phase
  LONE_gain = PARAM.LONE.GAIN;  % energy gain in LI02-LI10 (GeV)
  LONE_ampl = PARAM.LONE.FBAM;  % feedback amplitude (GV)
  
  % S10 chcn #s
  LI10_R56  = PARAM.LI10.R56;   % Measured val?
  %LI10_T566 = PARAM.LI10.T566;  % Measured val?
  LI10_ISR  = PARAM.LI10.ISR;   %
  LI10_ELO  = PARAM.LI10.ELO;   % S20 low energy cut
  LI10_EHI  = PARAM.LI10.EHI;   % S20 high energy cut
  
  % Energy gain and length of 02-10
  LTWO_leff = PARAM.LTWO.LEFF;  % Length of LI02-LI10 (m)
  LTWO_phas = PARAM.LTWO.PHAS;  % Chirp phase
  LTWO_gain = PARAM.LTWO.GAIN;
  LTWO_ampl = PARAM.LTWO.FBAM;  % feedback amplitude (GV)
  
  % S20 chcn R56 #s
  LI20_NLO  = PARAM.LI20.NLO;   % low energy notch edge
  LI20_NHI  = PARAM.LI20.NHI;   % high energy notch edge
  LI20_R56  = PARAM.LI20.R56;   % Measured val?
  LI20_T566 = PARAM.LI20.T566;  % Measured val?
  LI20_ISR  = PARAM.LI20.ISR;   %
  LI20_ELO  = PARAM.LI20.ELO;   % S20 low energy cut
  LI20_EHI  = PARAM.LI20.EHI;   % S20 high energy cut
  LI20_R16  = PARAM.LI20.R16;   % Dispersion
  LI20_T166 = PARAM.LI20.T166;  % Chromatic dispersion
  LI20_beta = PARAM.LI20.BETA;  % beta function
  LI20_emit = PARAM.LI20.EMIT;  % emittance
  
  % Initial beam parameters
  inp       = 'G';		        % gaussian Z and dE/E (see sigz0 =..., sigd0 =...)
  
  % 6mm bunches coming out of the ring teensy energy spread.
  sigz0     = PARAM.INIT.SIGZ0;	% rms bunch length used when inp=G or U above [m]
  sigd0     = PARAM.INIT.SIGD0;	% rms relative energy spread used when inp=G or U above [ ]
  z0_bar    = PARAM.INIT.Z0BAR; % axial offset of bunch [m] (used also with file input - mean of file removed first)
  d0_bar    = PARAM.INIT.D0BAR; % relative energy offset of bunch [ ]  (used also with file input - mean of file removed first)
  
  % 200K sim particles = 100K electrons per sim particle
  Nesim     = PARAM.INIT.NESIM;	% number of particles to generate for simulation when inp=G or U (reasonable: ~1000 to ~100000)
  Ne        = PARAM.INIT.NPART; % number of particles initially in bunch
  
  % The Holtzapple skew. Someday they'll name a skew after me. . .
  asym      = PARAM.INIT.ASYM;	% for inp='M' or 'G': sets rise/fall time width (-1<asym<1)
  
  % Our beam has no tail? That's a tall tale! Jesus I hope no one reads this. . .
  tail      = PARAM.INIT.TAIL;  % for inp='M' or 'G': sets rise/fall time width (0<=tail<1)
  cut       = PARAM.INIT.CUT;   % for inp='G': sets rise/fall time width (0.5<=cut<inf)
   
  % Other stuff
  splots    = PARAM.SIMU.PLOT;	% if =1, use small plots and show no wakes (for publish size plots)
  plot_frac = PARAM.SIMU.FRAC;  % fraction of particles to plot in the delta-z scatter-plots (0 < plot_frac <= 1)
  Nbin      = PARAM.SIMU.BIN;   % number of bins for z-coordinate (and dE/E for plots)
  gzfit     = PARAM.SIMU.ZFIT;  % if ==1: fit Z-distribution to gaussian (defaults to no-fit if 'gzfit' not provided)
  gdfit     = PARAM.SIMU.DFIT;  % if ==1: fit dE/E-distribution to gaussian (defaults to no-fit if 'gdfit' not provided)
  contf     = PARAM.SIMU.CONT;  % if ==1: get color contour image of z-d space (defaults to scatter plot if not provided)
  
  % S-band wavelength
  lambdaS   = 2.99792458e8/2856e6;

% The follwing array of file names, "wake_fn", is the point-charge wakefield filename(s) to be used.  The pointer
% to the used filename appears in the 5th column (wake ON/OFF) of the 'beamline' array below.  A "zero" (i.e. 0)
% means no wake used, and a value of j (e.g. 1,2,...) directs the calculation to use the jth point-charge wakefield
% file (i.e. the jth row of "wake_fn") for that part of the beamline.

wake_fn = ['slac.dat         '
           'slacx.dat        '
           'SlacL.dat        '
           'Slac_cu_rw.dat   '
           'SS_12700um_rw.dat'
           'Al_12700um_rw.dat'];		% name of point-charge wakefield file(s) ["slac.dat" is default if filename not given]

comment = 'FACET in Li20';	% text comment which appears at bottom of plots

% CODES:       |
%	           |		1				2				3				4				5				6
%==============|==============================================================================================
% compressor   |	code= 6           R56/m          T566/m          E_nom/GeV       U5666/m            -
% chicane      |	code= 7           R56/m         E_nom/GeV           -               -               -
% acceleration |	code=10       tot-energy/GeV    phase/deg        lambda/m   wake(ON=1,2**/OFF=0)  Length/m
% acceleration |	code=11     dEacc(phi=0)/GeV    phase/deg        lambda/m   wake(ON=1,2**/OFF=0)  Length/m
% E-spread add |	code=22          rms_dE/E           -               -               -               -
% E-window cut |	code=25         dE/E_window         -               -               -               -
% E-cut limits |	code=26          dE/E_min        dE/E_max           -               -               -
% con-N E-cut  |	code=27            dN/N             -               -               -               -
% Z-cut limits |	code=36            z_min           z_max            -               -               -
% con-N z-cut  |	code=37            dN/N             -               -               -               -
% modulation   |	code=44           mod-amp        lambda/m           -               -               -
% chromatic map|    code=98            R16/mm         T166/mm         beta/m       emittance/m          -
% STOP	       |	code=99             -               -               -               -               -
%=============================================================================================================
% CODE<0 makes a plot here, CODE>0 gives no plot here.



beamline = [
       11		0              0             lambdaS    0           0             % S-band
       11		NRTL_ampl      NRTL_phas     lambdaS    1           NRTL_leff     % Compressor cavity AMPL DR13 13 VDES
       26	    NRTL_ELO       NRTL_EHI      0          0           0             % Approximate energy acceptance of NRTL
       -6		NRTL_R56       NRTL_T566     E0         0           0             % Design NRTL ~0.603, BDES to KMOD for E-164 gives 0.588
       11		LONE_gain      LONE_phas     lambdaS    1           LONE_leff     % 2-6, nominal 9GeV, no feedback
       13       E1             LONE_ampl     -90        90          lambdaS       % Energy feedback to set 9GeV in chicane
       7	    LI10_R56       E1            0          0           0             % 2nd half of the chicane. Design was -0.0745, as built -0.076
       22       LI10_ISR       0             0          0           0             % Approximate SR growth in E-spread from chicane
       26       LI10_ELO       LI10_EHI      0          0           0             % Momentum Slits in FACET
       37		0.01           1             0          0           0             % Clip any rediculously long tails
      -11       LTWO_gain      LTWO_phas     lambdaS    1           LTWO_leff     % Boost to 23 GeV. 868m w/LCLS-II mods from P. Emma email 4-FEB-2011
       13       E2             LTWO_ampl     -90        90          lambdaS       % Energy feedback to set 20.35 GeV in S20 chicane
       6		LI20_R56       LI20_T566     E2         0           0             % FACET 'dogleg' like chicane
       22       LI20_ISR       0             0          0           0             % Approximate SR growth in E-spread from dogleg
       37       0.01           1             0          0           0             % Clip any rediculously long tails
       26       LI20_ELO       LI20_EHI      0          0           0             % Momentum Slits in FACET
       98       LI20_R16       LI20_T166     LI20_beta  LI20_emit   0             % Create chromatically dispersed distribution
       99	    0              0             0          0           0             % End
       ];
   


% Sign conventions used:
% =====================
%
% phase = 0 is beam on accelerating peak of RF (crest)
% phase < 0 is beam ahead of crest (i.e. bunch-head sees lower RF voltage than tail)
% The bunch head is at smaller values of z than the tail (i.e. head toward z<0)
% With these conventions, the R56 of a chicane is < 0 (and R56>0 for a FODO-arc) - note T566>0 for both
% 
% * = Note, the Energy-window cut (25-code) floats the center of the given FW window in order center it
%     on the most dense region (i.e. maximum number of particles).
% **  1:=1st wake file (e.g. S-band) is used, 2:=2nd wake file (e.g. X-band)