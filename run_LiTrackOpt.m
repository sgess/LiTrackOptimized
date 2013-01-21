clear all;

% Create Parameter struct
global PARAM;

% Set Paramete guess values
sim_params;

%Run LiTrack
OUT = LiTrackOpt('FACETpar');

%Plot Output
plot(OUT.X.AXIS,OUT.X.HIST);