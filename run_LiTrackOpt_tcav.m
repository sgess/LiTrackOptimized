clear all;

fontsize = 14;

show = 1;

% Load spectrum data
load('data_samples_1499.mat');
data_spectrum = SPECTRA(:,1)/sum(SPECTRA(:,1));
spectrum_axis = spectrum_axis/1000;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

% Set number of sim steps
ESsteps   = 30000;

% Initialize parameters
par_lims_tcav;
param_tcav;
params  = zeros(nPar,ESsteps);
pscaled = zeros(nPar,ESsteps); 
pCurrent = zeros(nPar,1);
pCurrent(1)  = PARAM.INIT.SIGZ0;    % Bunch Length
pCurrent(2)  = PARAM.INIT.SIGD0;    % Initial Energy Spread
pCurrent(3)  = PARAM.INIT.NPART;    % Number of Particles
pCurrent(4)  = PARAM.INIT.ASYM;     % Initial Gaussian Asymmetry
pCurrent(5)  = PARAM.NRTL.AMPL;     % Amplitude of RF Compressor
pCurrent(6)  = PARAM.NRTL.PHAS;     % RF Compressor Phase
pCurrent(7)  = PARAM.NRTL.R56;      % RTL compression
pCurrent(8)  = PARAM.NRTL.T566;     % RTL second order compression
pCurrent(9)  = decker;              % 2-10 Phase
pCurrent(10) = ramp;                % Ramp Phase
pCurrent(11) = PARAM.LI20.ELO;      % Low Energy Cutoff
pCurrent(12) = PARAM.LI20.EHI;      % High Energy Cutoff
pCurrent(13) = PARAM.LI20.BETA;     % Beta Function
pCurrent(14) = PARAM.LI20.R16;      % Dispersion
pCurrent(15) = PARAM.LI20.T166;     % 2nd Order Dispersion
pCurrent(16) = delta;               % Energy offset

% Record initial values
pInit = pCurrent;


% Initialize ES
[w, dt]   = init_ES(nPar);      % Get frequencies and time step
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
    %residual(j) = sum((sim_spectrum - data_spectrum).^2);

    % Set Cost as the value of the residual + particle fraction
    %cost(j) = 20 + log(residual(j)) + 0.001*Part_frac(j);
    cost(j) = 20 + log(residual(j));
    
    % ES Calc
    pLast = 2*(pCurrent-Cent)./Diff;
    pNext = pLast'+dt*cos(w*j*dt+gain*cost(j)).*(alpha*w).^0.5;
    pNext(pNext < -1) = -1;
    pNext(pNext >  1) =  1;
    pCurrent = Diff.*pNext'/2+Cent;
    
    
    % Update Params
    PARAM.INIT.SIGZ0 = pCurrent(1);   % Bunch Length
    PARAM.INIT.SIGD0 = pCurrent(2);   % Initial Energy Spread
    PARAM.INIT.NPART = pCurrent(3);   % Number of Particles
    PARAM.INIT.ASYM  = pCurrent(4);   % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL  = pCurrent(5);   % Amplitude of RF Compressor
    PARAM.NRTL.PHAS  = pCurrent(6);   % RF Compressor Phase
    PARAM.NRTL.R56   = pCurrent(7);   % RTL Compression
    PARAM.NRTL.T566  = pCurrent(8);   % RTL Second order compression
    decker           = pCurrent(9);   % 2-10 Phase
    ramp             = pCurrent(10);  % Ramp Phase
    PARAM.LI20.ELO   = pCurrent(11);  % Low Energy Cutoff
    PARAM.LI20.EHI   = pCurrent(12);  % High Energy Cutoff
    PARAM.LI20.BETA  = pCurrent(13);  % Beta Function
    PARAM.LI20.R16   = pCurrent(14);  % Dispersion
    PARAM.LI20.T166  = pCurrent(15);  % 2nd Order Dispersion
    delta            = pCurrent(16);  % Energy offset
    
    % Set dependent params
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    
    % Record evolving params
    params(:,j)  = pCurrent;
    pscaled(:,j) = pNext;
    
    
    if show
        
        figure(1);
        subplot(2,2,1);
        plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b','linewidth',2);
        axis([-5 5 0 3.5e-3]);
        xlabel('X (mm)','fontsize',12);
        title('Bunch Spectra','fontsize',10);
        legend('DATA','SIM');
        
        subplot(2,2,2);
        plot(1:ESsteps,residual,'color','r','linewidth',2);
        %axis([0 ESsteps 0 2e-4]);
        axis([0 ESsteps 0 4e-6]);
        xlabel('Step','fontsize',12);
        title('Residual','fontsize',10);
        
        subplot(2,2,3);
        plot(1:ESsteps,1-Part_frac,'color','g','linewidth',2);
        axis([0 ESsteps 0.75 1.1]);
        xlabel('Step','fontsize',12);
        title('Particle Fraction','fontsize',10);
        
        subplot(2,2,4);
        plot(1:ESsteps,cost,'color','b','linewidth',2);
        axis([0 ESsteps 0 12]);
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
    text(-3.5,5e-3,['Residual = ' num2str(residual(j),'%0.2e')],'fontsize',14);
    
    stupidest_plotter;
end