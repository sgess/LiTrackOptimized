% Simulation parameters
PARAM.SIMU.PLOT  = 0;       % Small plots?
PARAM.SIMU.FRAC  = 0.05;    % Fraction of sim particles to plot
PARAM.SIMU.BIN   = 256;     % Number of historgram bins
PARAM.SIMU.ZFIT  = 0;       % Longitudinal gauss fit
PARAM.SIMU.DFIT  = 0;       % Energy gauss fit
PARAM.SIMU.CONT  = 0;       % Contour plot

% Machine energy constraints
PARAM.ENRG.E0    = 1.19;    % Energy from ring (GeV)
PARAM.ENRG.E1    = 9.0;     % Energy at S10 (GeV)
PARAM.ENRG.E2    = 20.35;   % Energy at S20 (GeV)

% Beam initial conditions
PARAM.INIT.SIGZ0 = 7.00E-3; % RMS bunch length (m)
PARAM.INIT.SIGD0 = 8.00E-4; % RMS energy spread
PARAM.INIT.Z0BAR = 0;       % Z offset
PARAM.INIT.D0BAR = 0;       % Energy offset
PARAM.INIT.NESIM = 2E5;     % Number of simulated macro particles
PARAM.INIT.NPART = 2.00E10; % Number of electrons per bunch
PARAM.INIT.ASYM  = -0.150;  % The Holtzapple skew
PARAM.INIT.TAIL  = 0;       % Not sure what this is
PARAM.INIT.CUT   = 6;       % Not sure what this is

% NRTL bunch compressor
PARAM.NRTL.AMPL  = 0.0400;  % RTL compressor ampl (MV)
PARAM.NRTL.PHAS  = 90.00;   % RTL compressor phase (deg)
PARAM.NRTL.LEFF  = 2.1694;  % RTL cavity length (m)
PARAM.NRTL.R56   = 0.6026;  % RTL chicane R56 (m)
PARAM.NRTL.T566  = 1.053;   % RTL chicane T566 (m)
PARAM.NRTL.ELO   = -0.0300; % RTL lower momentum cut (GeV)
PARAM.NRTL.EHI   = 0.0300;  % RTL upper momentum cut (GeV)

% LI02-LI10 acceleration
PARAM.LONE.LEFF  = 809.5;   % Length of LI02-LI10 (m)
PARAM.LONE.CHRP  = 3.0536;  % chirp in 2-10 (GeV)
PARAM.LONE.RAMP  = 0.00;    % Ramped phase
PARAM.LONE.DECK  = -20.30;  % Staggered chirp
PARAM.LONE.PHAS  = PARAM.LONE.RAMP+PARAM.LONE.DECK;  % Total phase
PARAM.LONE.FBAM  = 0.235;   % feedback amplitude at S10 (GV)
PARAM.LONE.GAIN  = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

% LI10 bunch compressor
PARAM.LI10.R56   = -0.075786;% Sector 10 chicane R56 (m)
PARAM.LI10.ISR   = 5.9E-5;   % ISR energy spread from bends
PARAM.LI10.ELO   = -0.040;   % low energy cut
PARAM.LI10.EHI   = 0.040;    % high energy cut

% LI11-LI19 acceleration
PARAM.LTWO.LEFF  = 848;     % Length of LI02-LI10 (m)
PARAM.LTWO.CHRP  = 0;       % chirp in 2-10 (GeV)
PARAM.LTWO.PHAS  = PARAM.LONE.RAMP;    % 11-20 phase
PARAM.LTWO.FBAM  = 1.88;    % feedback amplitude at S20 (GV)
PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain

% LI20 bunch compressor
PARAM.LI20.NLO   = 0;       % notch low
PARAM.LI20.NHI   = 0;       % notch hi
PARAM.LI20.R56   = 0.0050;  % Sector 20 chicane R56 (m)
PARAM.LI20.T566  = 0.100;   % Sector 20 chicane T566 (m) % = 100 mm for R56 = 5mm from YS
PARAM.LI20.ISR   = 0.8E-5;  % ISR energy spread from bends
PARAM.LI20.ELO   = -0.035;  % RTL lower momentum cut (GeV)
PARAM.LI20.EHI   = 0.035;   % RTL upper momentum cut (GeV)
PARAM.LI20.R16   = 85;      % Dispersion at YAG
PARAM.LI20.T166  = 0.00;    % Second order dispersion at YAG
PARAM.LI20.BETA  = 3.0;     % Beta function at YAG
PARAM.LI20.EMIT  = 100e-6;  % Emittance in S20