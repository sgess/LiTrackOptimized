PARAM.MACH.RAMP  = 0;       % phase ramp

PARAM.ENRG.E0    = 1.19;    % Energy from ring (GeV)
PARAM.ENRG.E1    = 9.0;     % Energy at S10 (GeV)
PARAM.ENRG.E2    = 20.35;   % Energy at S20 (GeV)
%PARAM.ENRG.E2    = 23.00;   % Energy at S20 (GeV)

PARAM.INIT.SIGZ0 = 5.6E-3;  % RMS bunch length (m)
PARAM.INIT.SIGD0 = 7.39E-4; % RMS energy spread
PARAM.INIT.NESIM = 2E5;     % Number of simulated macro particles
PARAM.INIT.ASYM  = -0.245;  % The Holtzapple skew
PARAM.INIT.TAIL  = 0;       % Not sure what this is
PARAM.INIT.CUT   = 6;       % Not sure what this is

PARAM.NRTL.AMPL  = 0.0408;  % RTL compressor amplitude (GV)
PARAM.NRTL.PHAS  = 90;      % RTL compressor phase (deg)
PARAM.NRTL.LEFF  = 2.13;    % RTL cavity length (m)
PARAM.NRTL.R56   = 0.603;   % RTL chicane R56 (m)
PARAM.NRTL.T566  = 1.0535;  % RTL chicane T566 (m)
PARAM.NRTL.ELO   = -0.025;  % RTL lower momentum cut (GeV)
PARAM.NRTL.EHI   = 0.025;   % RTL upper momentum cut (GeV)

PARAM.LONE.LEFF  = 809.5;   % Length of LI02-LI10 (m)
%PARAM.LONE.PHAS  = -21.2;   % 2-10 phase for 'uniform' lattice
PARAM.LONE.PHAS  = -11.5275;% 2-10 phase for 'decker' lattice
%PARAM.LONE.GAIN  = 8.44;       % egain in 2-10, automatically set if 0 (GeV)
PARAM.LONE.GAIN  = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS);
PARAM.LONE.FBAM  = 0.235;   % feedback amplitude at S10 (GV)

PARAM.LI10.R56   = -0.076;  % Sector 10 chicane R56 (m)
%PARAM.LI10.T566  = 0.10;    % Sector 10 chicane T566 (m)
PARAM.LI10.ISR   = 5.9E-5;  % ISR energy spread from bends

%PARAM.LTWO.LEFF  = 868;     % Length of LI02-LI10 (m)
PARAM.LTWO.LEFF  = 848;     % Length of LI02-LI10 (m)
PARAM.LTWO.PHAS  = 0;       % 11-20 phase
PARAM.LTWO.FBAM  = 1.88;    % feedback amplitude at S20 (GV)

PARAM.LI20.R56   = 0.0030;   % Sector 20 chicane R56 (m)
PARAM.LI20.T566  = 0.10;    % Sector 20 chicane T566 (m)
PARAM.LI20.ISR   = 0.8E-5;  % ISR energy spread from bends
PARAM.LI20.ELO   = -0.03;   % RTL lower momentum cut (GeV)
PARAM.LI20.EHI   = 0.03;    % RTL upper momentum cut (GeV)

