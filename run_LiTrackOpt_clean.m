clear all;

fontsize = 14;

show = 1;

% Load spectrum data
load('data_samples_1108.mat');
data_spectrum = spec/sum(spec);
spectrum_axis = spectrum_axis/1000;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

% Initialize parameters
init_clean_params;


% Initialize ES
[w, dt]   = init_ES(nPar);      % Get frequencies and time step
ESsteps   = 3000;               % Set number of sim steps
alpha     = 2000;               % ES parameter
gain      = 2;                  % ES parameter
cost      = zeros(1,ESsteps);   % Cost
Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost
residual  = zeros(1,ESsteps);   % Chi2 difference between spectra

if show; figure(1); end;
tic;
for j=1:ESsteps;
    
    disp(j);
    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    
    % Calculate particle fraction
    Part_frac(j) = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;
    
    % Interpolate simulated spectrum
    sim_spectrum = interpSimSpec(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);
    
    % Calculate residual
    residual(j) = sum(data_spectrum.*(sim_spectrum - data_spectrum).^2);
    
    % Set Cost as the value of the residual + particle fraction
    cost(j) = 20 + log(residual(j)) + 0.001*Part_frac(j);
    
    % Set parameters
    set_clean_params;
    
    if show
        
        figure(1);
        subplot(2,2,1);
        plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b','linewidth',2);
        axis([-4 4 0 6e-3]);
        xlabel('X (mm)','fontsize',12);
        title('Bunch Spectra','fontsize',10);
        legend('DATA','SIM');
        
        subplot(2,2,2);
        plot(1:ESsteps,residual,'color','r','linewidth',2);
        axis([0 ESsteps 0 7e-6]);
        xlabel('Step','fontsize',12);
        title('Residual','fontsize',10);
        
        subplot(2,2,3);
        plot(1:ESsteps,1-Part_frac,'color','g','linewidth',2);
        axis([0 ESsteps 0.75 1.1]);
        xlabel('Step','fontsize',12);
        title('Particle Fraction','fontsize',10);
        
        subplot(2,2,4);
        plot(1:ESsteps,cost,'color','b','linewidth',2);
        axis([0 ESsteps 0 9]);
        xlabel('Step','fontsize',12);
        title('Cost','fontsize',10);
        
    end
    
end
toc;


% Plot Output
if show
    figure(2)
    plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b');
    legend('DATA','ES-SIMULATION');
    xlabel('X (mm)','fontsize',14);
    text(-3.5,5e-3,['Residual = ' num2str(residual(j-1),'%0.2e')],'fontsize',14);
end