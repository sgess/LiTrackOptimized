clear all;

fontsize = 14;

% Load sample spectra
load('data_samples.mat');
data_spectrum = SPECTRA(:,44)/sum(SPECTRA(:,44));

% Create Parameter struct
global PARAM;

% Parameter limit file
par_limits;

% Set Parameter guess values
sim_params;

w0=1500;
w00=50000;

w=zeros(1,18);

pr = primes(w0);
lpr = length(pr);

% The w(i) values are w(i) = sqrt(pi), where p1,...,p18 are the 18 primes less than w0. 
% This is a easy quick way to get pretty good "independence" between the
% varius frequencies.
for j2=18:-1:1;
    w(j2)=w00*(pr(lpr-5*j2).^0.5);
end

% ES Time Step Size, choose dt small enough so that it takes 20 steps for
% the highest frequency cos(w(18)n dt) to complete one full oscillation
dt=(2*pi)/(20*max(w));


% Total Number of Extremum Seeking Steps
ESsteps = 5000;

% Scale for New Cost Part
eps = 1e-3;

% ES Time, a purely digital entity
EST = ESsteps*dt;

% alpha is, in a way, the size of the perturbation, maybe want different values
% for different parameters, depending how sensitive they are
alpha = 2000*ones(1,18);
%alpha(2)=100;

% gain is the gain of each parameter's ES loop, maybe want different values
% for different parameters, depending how sensitive they are
gain = 80000*ones(1,18);

% Save end-of-simulation parameters to use as initial conditions for the
% next simulation
oldparams=zeros(1,18);

% Vector of 18 parameters that we will optimize

params=zeros(18,ESsteps);
pscaled=zeros(18,ESsteps);
cost=zeros(1,ESsteps);
avecost=zeros(1,ESsteps);

% Average of each of the 18 parameters and the cost
ave=zeros(18,ESsteps);
avec=zeros(1,ESsteps);

% Average derivatives of each of the 18 parameters and the cost
aved=zeros(18,ESsteps);
avecd=zeros(1,ESsteps);

    params(1,1) = PARAM.INIT.SIGZ0;     % Bunch Length
    params(2,1) = PARAM.INIT.SIGD0;     % Initial Energy Spread
    params(3,1) = PARAM.INIT.NPART;     % Number of Particles
    params(4,1) = PARAM.INIT.ASYM;      % Initial Gaussian Asymmetry
    params(5,1) = PARAM.NRTL.AMPL;      % Amplitude of RF Compressor
    params(6,1) = PARAM.NRTL.PHAS;      % RF Compressor Phase
    params(7,1) = PARAM.NRTL.ELO;       % Low Energy Cutoff
    params(8,1) = PARAM.NRTL.EHI;       % High Energy Cutoff
    params(9,1) = decker;               % 2-10 Phase
    params(10,1) = ramp;                % Ramp Phase
    params(11,1) = PARAM.LI10.ELO;      % Low Energy Cutoff
    params(12,1) = PARAM.LI10.EHI;      % High Energy Cutoff
    params(13,1) = PARAM.LI20.ELO;      % Low Energy Cutoff
    params(14,1) = PARAM.LI20.EHI;      % High Energy Cutoff
    params(15,1) = PARAM.LI20.BETA;     % Beta Function
    params(16,1) = PARAM.LI20.R16;      % Dispersion
    params(17,1) = PARAM.LI20.T166;     % 2nd Order Dispersion
    params(18,1) = delta;               % Energy offset
    
    % Initialize the simulation using the parameters it ended with last
    % time as the new initial parameters.
    
    % params(:,1)=oldparams(1,:);

tic



for j=1:ESsteps-1;
    j
    
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    Part_frac(j) = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;
    
    % Interpolate simulated spectrum
    sim_spectrum = interpSim(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);
    
    % Calculate residual
    residual = sum((sim_spectrum - data_spectrum).^2);
    
    % Set Cost as the value of the residual
    cost(j) = residual + eps*Part_frac(j);
    
    
    
    
    
    pscaled(:,j)=(params(:,j)-Cent)./Diff;
    pscaled(:,1)=0;
    
    
    % Smooth out parameters
    %if j>40;
    %    for jv = 1:18;
    %        ave(jv,j) = mean(pscaled(jv,j-40:j));
    %    end
    %    
    %    avec(j) = mean(cost(j-40:j));
    %end
    
    % Estimate slopes
    %if j>80;
    %    for jv = 1:18;
    %        aved(jv,j) = (mean(ave(jv,j-40:j))-mean(ave(jv,j-80:j-40)))/(80*dt);
    %    end
    %    
    %    avecd(:,j) = (mean(avec(j-40:j))-mean(avec(j-80:j-40)))/(80*dt); 
    %end
    
    
    
    for k = 1:18;
        pscaled(k,j+1)=pscaled(k,j)+dt*cos(w(k)*j*dt+gain(k)*cost(j))*(alpha(k)*w(k))^0.5;
        if pscaled(k,j+1) < -1;
            pscaled(k,j+1) = -1;
        else if pscaled(k,j+1) > 1;
                pscaled(k,j+1) = 1;
            end
        end
    end
    
    params(:,j+1)=Diff.*pscaled(:,j+1)+Cent;
    
    
    
    PARAM.INIT.SIGZ0 = params(1,j+1);           % Bunch Length
    PARAM.INIT.SIGD0 = params(2,j+1);           % Initial Energy Spread
    PARAM.INIT.NPART = params(3,j+1);           % Number of Particles
    PARAM.INIT.ASYM = params(4,j+1);            % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL = params(5,j+1);            % Amplitude of RF Compressor
    PARAM.NRTL.PHAS = params(6,j+1);            % RF Compressor Phase
    PARAM.NRTL.ELO = params(7,j+1);             % Low Energy Cutoff
    PARAM.NRTL.EHI = params(8,j+1);             % High Energy Cutoff
    decker = params(9,j+1);                     % 2-10 Phase
    ramp = params(10,j+1);                      % Ramp Phase
    PARAM.LI10.ELO = params(11,j+1);            % Low Energy Cutoff
    PARAM.LI10.EHI = params(12,j+1);            % High Energy Cutoff
    PARAM.LI20.ELO = params(13,j+1);            % Low Energy Cutoff
    PARAM.LI20.EHI = params(14,j+1);            % High Energy Cutoff
    PARAM.LI20.BETA = params(15,j+1);           % Beta Function
    PARAM.LI20.R16 = params(16,j+1);            % Dispersion
    PARAM.LI20.T166 = params(17,j+1);           % 2nd Order Dispersion
    delta = params(18,j+1);                     % Energy offset
       
end
toc

% Save parameters to use as new starting point
oldparams(1,:)=params(:,j-1);


% Plot Output
figure(1)
plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b');
legend('DATA','ES-SIMULATION');
xlabel('X (mm)','fontsize',14);
text(-3.5,5e-3,['Residual = ' num2str(residual,'%0.2e')],'fontsize',14);