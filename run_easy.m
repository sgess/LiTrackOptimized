clear all;
%load('concat_1111.mat');
%load('concat_1103.mat');
%load('../DATA/nas/nas-li20-pm01/E200/2013/20130428/E200_10794/slim.mat');
load('/Users/sgess/Desktop/data/2013/slims/slim_10794.mat');
%load('/Users/sgess/Desktop/data/2013/slims/slim_10915.mat');
%spec_axis = cat_dat.yag_ax;
%spec_thing = cat_dat.YAG_SPEC(:,1);
spec_axis = data.YAG.axis;
%spec_thing = data.YAG.spectra(:,1);
spec_thing = data.YAG.spectra(:,1);
%spec_thing = mean(cat_dat.YAG_SPEC,2);
show = 1;
nOut = 3;

%Load wakefield data
global A;
A = load('slac.dat');

%Create Parameter struct
global PARAM;
param_04_16_13;
PARAM.LI20.R16 = 90;
PARAM.LI20.BETA = 3;

pars_init = [0.0068;        0.00084;         2.05e10;       -0.147];
sens_init = [0.10;          0.05;            0.08;           0.05];
name_init = {'INIT SIGZ0'; 'INIT SIGD0'; 'INIT NPART'; 'INIT ASYM'};

pars_nrtl = [0.0395;        90.3;         0.6026;       1.053];
sens_nrtl = [0.03;          0.005;           0.03;       0.300];
name_nrtl = {'NRTL AMPL'; 'NRTL PHAS'; 'NRTL R56'; 'NRTL T566'};

% pars_init = [0.0066;        2.1e10];
% sens_init = [0.2;               0.1];
% name_init = {'INIT SIGZ0'; 'INIT NPART'};
% 
% pars_nrtl = [0.0403;        90.30];
% sens_nrtl = [0.06;           0.01];
% name_nrtl = {'NRTL AMPL'; 'NRTL PHAS'};

pars_lone = [-20.3;              0];
sens_lone = [0.03;              -5];
name_lone = {'LONE DECK'; 'LONE RAMP'};

%pars_ltwo = [0.0];
%sens_ltwo = [-5];
%name_ltwo = {'LTWO PHAS'};
% 
% pars_li20 = [4;             84;             -100;        0.030;   -0.030];
% sens_li20 = [0.5;           0.01;            1;         0.5;     0.5];
% name_li20 = {'LI20 BETA'; 'LI20 R16'; 'LI20 T166'; 'LI20 EHI'; 'LI20 ELO'};
% 
pars_li20 = [90; 0;               ];
sens_li20 = [0.05; -500;              ];
name_li20 = {'LI20 R16'; 'LI20 T166';};

%pars_li20 = [0.005;           0.1];
%sens_li20 = [.1;              -0.2];
%name_li20 = {'LI20 R56';  'LI20 T566'}; 
 
pars = [pars_init; pars_nrtl; pars_lone; pars_li20];
sens = [sens_init; sens_nrtl; sens_lone; sens_li20];
name = [name_init; name_nrtl; name_lone; name_li20];

% pars = [pars_init; pars_nrtl; pars_lone; pars_ltwo];
% sens = [sens_init; sens_nrtl; sens_lone; sens_ltwo];
% name = [name_init; name_nrtl; name_lone; name_ltwo];

% pars = [pars_init; pars_nrtl; pars_lone; pars_ltwo; pars_li20];
% sens = [sens_init; sens_nrtl; sens_lone; sens_ltwo; sens_li20];
% name = [name_init; name_nrtl; name_lone; name_ltwo; name_li20];

%pars = pars_li20;
%sens = sens_li20;
%name = name_li20;

% Set number of sim steps
ESsteps   = 10000;

nPar = length(pars);
[Cent, Diff, lo_lims, hi_lims] = SetParLims(pars,sens);
SetPars(pars, name, nPar);

pInit    = pars;
pCurrent = pars;

params  = zeros(nPar,ESsteps);
pscaled = zeros(nPar,ESsteps);

% Initialize ES
[w, dt]   = init_ES(nPar);      % ES frequencies and time step
alpha     = 1000;                % ES parameter
gain      = 4.5e-1;               % ES parameter
cost      = zeros(1,ESsteps);   % ES cost
Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost
I_peak    = zeros(1,ESsteps);
residual  = zeros(1,ESsteps);   % Chi2 difference between spectra

if show; figure(1); end;

j = 0;

% Calculate axes
xx = spec_axis;
Lineout = spec_thing;

Line_minBG = Lineout-Lineout(1);
line_x  = xx;
x_avg = mean(line_x);
[MaxLine,max_ind] = max(Line_minBG);
SumLine = sum(Line_minBG);
center = sum(line_x.*Line_minBG)/sum(Line_minBG);
Line_minBG = Line_minBG/MaxLine;



while j <= ESsteps
    
    display(j);
    
    %%%%%%%%%%%%%%%%%%%%
    % Simulation Stuff %
    %%%%%%%%%%%%%%%%%%%%
    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    
    % Calculate particle fraction
    Part_frac(j+1) = 1 - OUT.I.PART(nOut)/PARAM.INIT.NESIM;
    I_peak(j+1)    = OUT.I.PEAK(nOut);
    
    % Interpolate simulated spectrum
    SimDisp = interpSimX(OUT,line_x,PARAM.SIMU.BIN,center-x_avg);
    SumX    = sum(SimDisp);
    normX   = SumLine/SumX;
    ProfXLi = normX*SimDisp;
    [MaxSim,sim_ind] = max(ProfXLi);
    ProfXLi = ProfXLi/MaxSim;
    
    % Get bunch profile
    dZ      = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
    zzLi    = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
    ProfZLi = [0; OUT.Z.HIST(:,nOut); 0];
    
    % Calculate residual
    %residual(j+1) = sum(Line_minBG.*(ProfXLi - Line_minBG).^2);
    residual(j+1) = sum((ProfXLi - Line_minBG).^2);

    % Set Cost as the value of the residual + particle fraction
    %cost(j+1) = residual(j+1)+1000000/OUT.I.PEAK(3);
    %cost(j+1) = 1e4*(23 - OUT.I.PEAK(3));
    cost(j+1) = residual(j+1);
    
    % ES Calc
    pLast = 2*(pCurrent-Cent)./Diff;
    pNext = pLast'+dt*cos(w*(j+1)*dt+gain*cost(j+1)).*(alpha*w).^0.5;
    pNext(pNext < -1) = -1;
    pNext(pNext >  1) =  1;
    pCurrent = Diff.*pNext'/2+Cent;
    
    % Update Params
    SetPars(pCurrent, name, nPar);
    
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
            plot(1:ESsteps,I_peak,'color','g','linewidth',2);
            part_plot = get(gca,'Children');
            xlabel('Step','fontsize',12);
            title('Peak Current','fontsize',10);

            subplot(2,2,4);
            plot(zzLi,ProfZLi,'color','r','linewidth',2);
            prof_plot = get(gca,'Children');
            xlabel('Z (mm)');
            title('Longitudinal Beam Profile');
        
        else
            
            set(spectra(1),'XData',line_x,'YData',Line_minBG);
            set(spectra(2),'XData',line_x,'YData',ProfXLi);
            set(res_plot,'YData',residual);
            set(part_plot,'YData',I_peak);
            set(prof_plot,'XData',zzLi,'YData',ProfZLi);
            pause(0.0001);
            
        end
            
            
    end
    
    j = j + 1;
       
end