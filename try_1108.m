%clear all;

%load('retry_1108.mat');
%load('retry_1103.mat');
spec_axis = DATA.AXIS.xx/1000;
spec_thing = DATA.YAG.spectrum(:,59);

%addpath(genpath('LiTrack'));
fontsize = 14;


show = 1;
nOut = 2;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

% Set number of sim steps
ESsteps   = 10000;

% Initialize parameters
par_lims_retry1108;
param_tcav;
params  = zeros(nPar,ESsteps);
pscaled = zeros(nPar,ESsteps); 

use_new = 0;

if use_new
    
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
    pCurrent(10) = l_two;              % 11-20 Phase
    pCurrent(11) = ramp;                % Ramp Phase
    pCurrent(12) = PARAM.LI20.BETA;     % Beta Function
    pCurrent(13) = PARAM.LI20.R16;      % Dispersion
    pCurrent(14) = PARAM.LI20.T166;     % 2nd Order Dispersion
    
    PARAM.LONE.PHAS = decker+ramp;  % Total PhasepCurrent(14)
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    
    PARAM.LTWO.PHAS = l_two+ramp;  % Total Phase
    PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain

else
    
    PARAM.INIT.SIGZ0 = pCurrent(1);   % Bunch Length
    PARAM.INIT.SIGD0 = pCurrent(2);   % Initial Energy Spread
    PARAM.INIT.NPART = pCurrent(3);   % Number of Particles
    PARAM.INIT.ASYM  = pCurrent(4);   % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL  = pCurrent(5);   % Amplitude of RF Compressor
    PARAM.NRTL.PHAS  = pCurrent(6);   % RF Compressor Phase
    PARAM.NRTL.R56   = pCurrent(7);   % RTL Compression
    PARAM.NRTL.T566  = pCurrent(8);   % RTL Second order compression
    decker           = pCurrent(9);   % 2-10 Phase
    l_two            = pCurrent(10);  % 11-20 Phase
    ramp             = pCurrent(11);  % Ramp Phase    
    PARAM.LI20.BETA  = pCurrent(12);  % Beta Function
    PARAM.LI20.R16   = pCurrent(13);  % Dispersion
    PARAM.LI20.T166  = pCurrent(14);  % 2nd Order Dispersion

    % Set dependent params
    PARAM.LONE.PHAS = decker+ramp;  % Total PhasepCurrent(14)
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    
    PARAM.LTWO.PHAS = l_two+ramp;  % Total Phase
    PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain
    
end

% Record initial values
pInit = pCurrent;
    
% Initialize ES
[w, dt]   = init_ES(nPar);      % ES frequencies and time step
alpha     = 500;               % ES parameter
gain      = 5*1600e-11;        % ES parameter
cost      = zeros(1,ESsteps);   % ES cost
Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost
residual  = zeros(1,ESsteps);   % Chi2 difference between spectra


if show; figure(1); end;

j = 0;
k = 0;


% Calculate axes
xx = spec_axis';
Lineout = spec_thing;

Line_minBG = Lineout-Lineout(1);
line_x  = xx;
x_avg = mean(line_x);
[MaxLine,max_ind] = max(Line_minBG);
SumLine = sum(Line_minBG);
center = sum(line_x.*Line_minBG)/sum(Line_minBG);




while j <= ESsteps
    
    display(j);
       
    del = 100*center/pCurrent(13);

    % Calculate delta axis
    dd = 100*(line_x-center)/pCurrent(13);
    
    
    %%%%%%%%%%%%%%%%%%%%
    % Simulation Stuff %
    %%%%%%%%%%%%%%%%%%%%
    
    % Run LiTrack
    try
        OUT = LiTrackOpt('FACETpar');
    catch err
        display('LiTracked failed');
        pCurrent = params(:,j-1);
    end
    
    % Calculate particle fraction
    Part_frac(j+1) = 1 - OUT.I.PART(nOut)/PARAM.INIT.NESIM;
    
    % Interpolate simulated spectrum
    SimDisp = interpSimX(OUT,line_x,PARAM.SIMU.BIN,center-x_avg);
    SumX = sum(SimDisp);
    normX = SumLine/SumX;
    ProfXLi = normX*SimDisp;
    [MaxSim,sim_ind] = max(ProfXLi);
    
    % Get bunch profile
    dZ = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
    zzLi = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
    ProfZLi = [0; OUT.Z.HIST(:,nOut); 0];
    
    % Calculate residual
    %residual(j) = sum(Line_minBG.*(ProfXLi' - Line_minBG).^2);
    residual(j+1) = sum(Line_minBG.*(ProfXLi - Line_minBG).^2);

    % Set Cost as the value of the residual + particle fraction
    cost(j+1) = residual(j+1);
    
    % ES Calc
    pLast = 2*(pCurrent-Cent)./Diff;
    pNext = pLast'+dt*cos(w*(j+1)*dt+gain*cost(j+1)).*(alpha*w).^0.5;
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
    l_two            = pCurrent(10);  % 11-20 Phase
    ramp             = pCurrent(11);  % Ramp Phase    
    PARAM.LI20.BETA  = pCurrent(12);  % Beta Function
    PARAM.LI20.R16   = pCurrent(13);  % Dispersion
    PARAM.LI20.T166  = pCurrent(14);  % 2nd Order Dispersion

    % Set dependent params
    PARAM.LONE.PHAS = decker+ramp;  % Total PhasepCurrent(14)
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    
    PARAM.LTWO.PHAS = l_two+ramp;  % Total Phase
    PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain
    
    % Record evolving params
    params(:,j+1)  = pCurrent;
    pscaled(:,j+1) = pNext;
    

    
    if show
        
        if j==0
        

        
            figure(1);
            h1 = subplot(2,2,1);
            plot(line_x,Line_minBG,'g',line_x,ProfXLi,'b','linewidth',2);
            spectra = get(gca,'Children');
            ax =  get(gca);
            xlabel('X (mm)','fontsize',12);
            title('Bunch Spectra','fontsize',10);
            legend('SIM','DATA');

            subplot(2,2,2);
            plot(1:ESsteps,residual,'color','r','linewidth',2);
            res_plot = get(gca,'Children');
            xlabel('Step','fontsize',12);
            title('Residual','fontsize',10);
            
            subplot(2,2,3);
            plot(1:ESsteps,1-Part_frac,'color','g','linewidth',2);
            part_plot = get(gca,'Children');
            axis([0 ESsteps 0.75 1.1]);
            xlabel('Step','fontsize',12);
            title('Particle Fraction','fontsize',10);

            subplot(2,2,4);
            plot(zzLi,ProfZLi,'color','r','linewidth',2);
            prof_plot = get(gca,'Children');
            xlabel('Z (mm)');
            title('Longitudinal Beam Profile');
        
        else

            
            set(spectra(1),'XData',line_x,'YData',Line_minBG);
            set(spectra(2),'XData',line_x,'YData',ProfXLi);
            %ax_max =  max(MaxLine,MaxSim)+10;
            %new_ax = [LineLim(3) LineLim(4) 0 ax_max];
            %axis(h1,new_ax);
            
            set(res_plot,'YData',residual);
            set(part_plot,'YData',1-Part_frac);
            set(prof_plot,'XData',zzLi,'YData',ProfZLi);
            pause(0.0001);
            
        end
            
            
    end
    
    j = j + 1;
    
    
end
