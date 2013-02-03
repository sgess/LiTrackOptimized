clear all;

fontsize = 14;

% Load sample profile
x = randn(1,100000);
y = randn(1,50000);
d = [25*x-100, 25*y+250];
[n,ax]=hist(d,128);
ax_min = ax(1);
dax = (ax(2) - ax(1))/1000;
pMin = -0.5;
pMax = 1.0;

samp_ax = (ax-ax_min)'/1000;
prof = n'/sum(n);

% load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

% Parameter limit file
%par_limits;
%new_par_lims;
%newer_par_lims;
notch_par_lims;

% Set Parameter guess values
%sim_params;
notch_params;

n_par = 18;

w0=1000;
w00=5000;

w=zeros(1,n_par);

pr = primes(w0);
lpr = length(pr);

% The w(i) values are w(i) = sqrt(pi), where p1,...,p17 are the 17 primes less than w0. 
% This is a easy quick way to get pretty good "independence" between the
% varius frequencies.
for j2=n_par:-1:1;
    w(j2)=w00*(pr(lpr-5*j2).^0.5);
end

% ES Time Step Size, choose dt small enough so that it takes 20 steps for
% the highest frequency cos(w(17)n dt) to complete one full oscillation
dt=(2*pi)/(8*w(n_par));


% Total Number of Extremum Seeking Steps
ESsteps = 3000;

% ES Time, a purely digital entity
EST = ESsteps*dt;

% alpha is, in a way, the size of the perturbation, maybe want different values
% for different parameters, depending how sensitive they are
alpha = 2000*ones(1,n_par);

% gain is the gain of each parameter's ES loop, maybe want different values
% for different parameters, depending how sensitive they are
gain = 2*ones(1,n_par);


% Vector of 17 parameters that we will optimize

params=zeros(n_par,ESsteps);
pscaled=zeros(n_par,ESsteps);
cost=zeros(1,ESsteps);
Part_frac=zeros(1,ESsteps);
residual=zeros(1,ESsteps);

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
    params(15,1) = PARAM.LI20.R56;
    params(16,1) = PARAM.LI20.NLO;
    params(17,1) = PARAM.LI20.NHI;
    params(18,1) = dZ;
tic

figure(1);
%subplot(2,2,1);

for j=1:ESsteps-1;
    j
    
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    PARAM.LI20.T566  = p(1)*PARAM.LI20.R56^2 + p(2)*PARAM.LI20.R56 + p(3);
    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    Part_frac(j) = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;
    
    
    % Interpolate simulated spectrum
    sim_profile = interpSimProf(OUT,PARAM.SIMU.BIN,2,dax,pMin,pMax,dZ);
    
    % Calculate residual
    %residual(j) = sum((sim_spectrum - data_spectrum).^2);
    residual(j) = sum(prof.*(sim_profile - prof).^2);
    
    % Set Cost as the value of the residual
    %cost(j) = residual;
    %cost(j) = 14 + log(residual(j)) + 0.001*Part_frac(j);
    %cost(j) = 20 + log(residual(j)) + 0.001*Part_frac(j);
    cost(j) = 20 + log(residual(j));
    
    pscaled(:,j)=2*(params(:,j)-Cent)./Diff;
    
    for k = 1:n_par;
        pscaled(k,j+1)=pscaled(k,j)+dt*cos(w(k)*j*dt+gain(k)*cost(j))*(alpha(k)*w(k))^0.5;
        if pscaled(k,j+1) < -1;
            pscaled(k,j+1) = -1;
        elseif pscaled(k,j+1) > 1;
            pscaled(k,j+1) = 1;
        end
        if pscaled(16,j+1) > pscaled(17,j+1)
            pscaled(16,j+1) = pscaled(17,j+1) - 0.0000001;
        end
    end
    
    params(:,j+1)=Diff.*pscaled(:,j+1)/2+Cent;
    
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
    PARAM.LI20.R56 = params(15,j+1);
    PARAM.LI20.NLO = params(16,j+1);
    PARAM.LI20.NHI = params(17,j+1);
    dZ = params(18,j+1);
    
    figure(1);
    subplot(2,2,1);
    plot(samp_ax,prof,'g',samp_ax,sim_profile,'b','linewidth',2);
    %axis([-4 4 0 4e-3]);
    xlabel('X (mm)','fontsize',12);
    title('Bunch Spectra','fontsize',10);
    legend('DATA','SIM');
    
    subplot(2,2,2);
    plot(1:ESsteps,residual,'color','r','linewidth',2);
    %axis([0 ESsteps 0 7e-6]);
    xlabel('Step','fontsize',12);
    title('Residual','fontsize',10);
    
    subplot(2,2,3);
    plot(1:ESsteps,1-Part_frac,'color','g','linewidth',2);
    %axis([0 ESsteps 0.75 1.1]);
    xlabel('Step','fontsize',12);
    title('Particle Fraction','fontsize',10);
    
    subplot(2,2,4);
    plot(1:ESsteps,cost,'color','b','linewidth',2);
    %axis([0 ESsteps 0 9]);
    xlabel('Step','fontsize',12);
    title('Cost','fontsize',10);
    
    %saveas(gca,['/Users/sgess/Desktop/plots/MOVIES/ES/short/k1_' num2str(j,'%03d') '.png']);
%     figure(1);
%     plot(1:ESsteps,cost,'b',1:ESsteps,residual*1e6,'r',1:ESsteps,(1-Part_frac),'g');
%     axis([0 ESsteps 0 9]);
    
    
    
    
end
toc


% Plot Output
figure(1)
plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b');
legend('DATA','ES-SIMULATION');
xlabel('X (mm)','fontsize',14);
text(-3.5,5e-3,['Residual = ' num2str(residual(j-1),'%0.2e')],'fontsize',14);
