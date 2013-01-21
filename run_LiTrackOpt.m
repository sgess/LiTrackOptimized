clear all;

% Load sample spectra
load('data_samples.mat');
data_spectrum = SPECTRA(:,1)/sum(SPECTRA(:,1));

% Create Parameter struct
global PARAM;

% Set Paramete guess values
sim_params;

% Run LiTrack
OUT = LiTrackOpt('FACETpar');

% Interpolate simulated spectrum 
sim_spectrum = interpSim(OUT,spectrum_axis);

residual = sum((sim_spectrum - data_spectrum).^2);

% Plot Output
%plot(spectrum_axis,sim_spectrum,'color','b');
%hold on;
%plot(spectrum_axis,data_spectrum,'color','g');
plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b');
legend('DATA','SIMULATION','fontsize',14);
xlabel('X (mm)','fontsize',14);
text(-3.5,5e-3,['Residual = ' num2str(residual,'%0.2e')],'fontsize',14);

