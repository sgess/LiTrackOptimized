clear all;

show = 0;
nOut = 3;

global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

%param_tcav;
param_04_16_13;
PARAM.LI20.R16   = 85;
PARAM.LI20.BETA  = 4.0;
PARAM.LI20.T166  = 0;

pars_init = [0.0068;        0.0008;         2.1e10;       -0.15];
sens_init = [0.2;           0.2;            0.1;           0.3];
name_init = {'INIT SIGZ0'; 'INIT SIGD0'; 'INIT NPART'; 'INIT ASYM'};

pars_nrtl = [0.0400;        90.30;         0.602;       1.3];
sens_nrtl = [0.06;           0.01;            0.02;       0.1];
name_nrtl = {'NRTL AMPL'; 'NRTL PHAS'; 'NRTL R56'; 'NRTL T566'};

pars_lone = [-21.6];
sens_lone = [0.3];
name_lone = {'LONE PHAS'};

pars_ltwo = [-5];
sens_ltwo = [1];
name_ltwo = {'LTWO PHAS'};

pars = [pars_init; pars_nrtl; pars_lone; pars_ltwo];
sens = [sens_init; sens_nrtl; sens_lone; sens_ltwo];
name = [name_init; name_nrtl; name_lone; name_ltwo];

% Set number of sim steps
ESsteps   = 500;

nPar = length(pars);
[Cent, Diff, lo_lims, hi_lims] = SetParLims(pars,sens);
SetPars(pars, name, nPar);

pInit    = pars;
pCurrent = pars;

params  = zeros(nPar,ESsteps);
pscaled = zeros(nPar,ESsteps);

% Initialize ES
[w, dt]   = init_ES(nPar);      % ES frequencies and time step
alpha     = 500;                % ES parameter
gain      = 4.5e-5;               % ES parameter
cost      = zeros(1,ESsteps);   % ES cost
Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost


if show; figure(1); end;

% load data, get axis
load('/Users/sgess/Desktop/data/2013/slims/slim_10794.mat');
line_x = data.YAG.axis;
x_avg = mean(line_x);
specs = data.YAG.spectra(:,data.YAG.good_shot);
n_shots = sum(data.YAG.good_shot);

for i = 1:n_shots
    
    % Calculate axes
    Lineout = specs(:,i);
    Line_minBG = Lineout-Lineout(1);
    [MaxLine,max_ind] = max(Line_minBG);
    SumLine = sum(Line_minBG);
    center = sum(line_x.*Line_minBG)/sum(Line_minBG);
    I_peak    = zeros(1,ESsteps);
    residual  = zeros(1,ESsteps);   % Chi2 difference between spectra
    
    pInit = pCurrent;
    
    params  = zeros(nPar,ESsteps);
    pscaled = zeros(nPar,ESsteps);
    
    j = 0;
    
    while j <= ESsteps
        
        display(['Shot ' num2str(i) 'ESstep ' num2str(j)]);
        
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
        
        % Get bunch profile
        dZ      = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
        zzLi    = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
        ProfZLi = [0; OUT.Z.HIST(:,nOut); 0];
        
        % Calculate residual
        %residual(j+1) = sum(Line_minBG.*(ProfXLi - Line_minBG).^2);
        residual(j+1) = sum((ProfXLi - Line_minBG).^2);
        
        if residual(j+1) < 2e4
            shot_par(i) = PARAM;
            ress(i) = residual(j+1);
            ipk(i) = I_peak(j+1);
            break;
        end
        
        % Set Cost as the value of the residual + particle fraction
        %cost(j+1) = residual(j+1)+1000000/OUT.I.PEAK(3);
        %cost(j+1) = 1e4*(20 - OUT.I.PEAK(3));
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
    [r,b] = min(residual);       
    if b == 1
        SetPars(pInit,name, nPar);
        ipk(i) = I_peak(1);
    else
        SetPars(params(:,b-1),name, nPar);
        ipk(i) = I_peak(b-1);
    end
    shot_par(i) = PARAM;
    ress(i) = r;
end