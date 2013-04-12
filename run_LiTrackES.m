clear all;

addpath(genpath('LiTrack'));
fontsize = 14;

% Lock SIOC PVs
ismine = round(1000*rand());
lcaput('SIOC:SYS1:ML00:AO815',ismine);

% Trig PV
trigPV='TRIG:LI02:813';

% get machine data
state = facet_LongMachine();
MACH = MACH_ANA(state);

useCurrent = 1;
show = 1;
n_avg = 1;
n_samp = 10;
nOut = 3;

% Load wakefield data
global A;
A = load('slac.dat');

% Create Parameter struct
global PARAM;

% Set number of sim steps
ESsteps   = 10000;

% Initialize parameters
par_lims_03_30_13;
param_03_30_13;
params  = zeros(nPar,ESsteps);
pscaled = zeros(nPar,ESsteps); 
pCurrent = zeros(nPar,1);
pCurrent(1)  = 6.0e-3;    % Bunch Length
pCurrent(2)  = PARAM.INIT.SIGD0;    % Initial Energy Spread
pCurrent(3)  = 9.4e9;    % Number of Particles
pCurrent(4)  = PARAM.INIT.ASYM;     % Initial Gaussian Asymmetry
pCurrent(5)  = 0.038;     % Amplitude of RF Compressor
pCurrent(6)  = PARAM.NRTL.PHAS;     % RF Compressor Phase
pCurrent(7)  = PARAM.NRTL.R56;      % RTL compression
pCurrent(8)  = PARAM.NRTL.T566;     % RTL second order compression
pCurrent(9)  = -17.3;              % 2-10 Phase
pCurrent(10)  = 0;              % 11-20 Phase
pCurrent(11) = 0;                % Ramp Phase
pCurrent(12) = PARAM.LI20.BETA;     % Beta Function
pCurrent(13) = 115;      % Dispersion
pCurrent(14) = 100;     % 2nd Order Dispersion

if useCurrent
    pCurrent = lcaGetSmart(strcat('SIOC:SYS1:ML00:AO8',{'01','02','03','04','05','06','07','08','09','10','11','12','13','14'}));
end

% Record initial values
pInit = pCurrent;

% Record other interesting shit
interesting_shit = zeros(14,ESsteps);
profiles = zeros(PARAM.SIMU.BIN+2,round(ESsteps/n_samp));
profaxes = zeros(PARAM.SIMU.BIN+2,round(ESsteps/n_samp));
yagspecs = zeros(1392,round(ESsteps/n_samp));
Li_specs = zeros(1392,round(ESsteps/n_samp));
specaxes = zeros(1392,round(ESsteps/n_samp));
index = zeros(round(ESsteps/n_samp));

    
% Initialize ES
[w, dt]   = init_ES(nPar);      % ES frequencies and time step
alpha     = 1000;               % ES parameter
gain      = 10*1600e-10;        % ES parameter
cost      = zeros(1,ESsteps);   % ES cost
Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost
residual  = zeros(1,ESsteps);   % Chi2 difference between spectra

% Period of oscillation
%period = round(20*max(w)/min(w));
period = 56;

if show; figure(1); figure(2); end;

j = 0;
k = 0;

while j <= ESsteps
    
    display(j);
    
    %%%%%%%%%%%%%%%
    % Image Stuff %
    %%%%%%%%%%%%%%%

    % Test YAG to make sure it is not being used in a DAQ
    daq_stat = lcaGetSmart('YAGS:LI20:2432:ENABLE_DAQ');

    % Get YAG resolution
    res = lcaGetSmart('YAGS:LI20:2432:RESOLUTION');

    % Get YAG dispersion
    disp = lcaGetSmart('SIOC:SYS1:ML00:AO755');

    % Get YAG lineout limits
    LineLim = lcaGetSmart({'SIOC:SYS1:ML00:AO751' 'SIOC:SYS1:ML00:AO752' 'SIOC:SYS1:ML00:AO753' 'SIOC:SYS1:ML00:AO754'});
    
    % Get beam rate
    rate = lcaGetSmart('EVNT:SYS1:1:BEAMRATE');
    
    % Get averaging and sampling
    n_avg = lcaGetSmart('SIOC:SYS1:ML00:AO824');
    n_samp = lcaGetSmart('SIOC:SYS1:ML00:AO825');
    
    % Get 2-9 status
    trigVal=control_klysStatGet(trigPV);
    if j == 0; rate = 1; trigVal = 0; end;
    
    % Get YAG image
    if mod(j,n_samp) == 0
        if strcmp(daq_stat,'Disabled') && rate ~= 0 && trigVal == 0
            for a=1:n_avg
                im_dat = profmon_grab('YAGS:LI20:2432');
                if a == 1;
                    img = im_dat.img;
                else
                    img = im_dat.img + img;
                end
                pause(1/rate);
            end
            img = img/n_avg;
        end
    end
       
    % Calculate axes
    xx = im_dat.res(1)/1000*((1:im_dat.roiXN)-im_dat.centerX+im_dat.roiX);
    yy = -im_dat.res(2)/1000*((1:im_dat.roiYN)-im_dat.centerY+im_dat.roiY);

    % Match lineout limits to axes
    if LineLim(1) < yy(end); LineLim(1) = yy(end)+0.1; end;
    if LineLim(2) > yy(1); LineLim(2) = yy(1)-0.1; end;
    if LineLim(3) < xx(1); LineLim(3) = xx(1)+0.1; end;
    if LineLim(4) > xx(end); LineLim(4) = xx(end)-0.1; end;

    % Map lineout to pixels
    PixLim(1) = round(-LineLim(1)*1000/im_dat.res(2))+im_dat.centerY-im_dat.roiY;
    PixLim(2) = round(-LineLim(2)*1000/im_dat.res(2))+im_dat.centerY-im_dat.roiY;
    PixLim(3) = round(LineLim(3)*1000/im_dat.res(1))+im_dat.centerX-im_dat.roiX;
    PixLim(4) = round(LineLim(4)*1000/im_dat.res(1))+im_dat.centerX-im_dat.roiX;

    % Make x axis have even number of pixels
    if mod(PixLim(4) - PixLim(3),2) == 0
        PixLim(4) = PixLim(4)+1;
    end

    % Calculate lineout
    Lineout = mean(img(PixLim(2):PixLim(1),PixLim(3):PixLim(4)),1);
    Line_minBG = Lineout-Lineout(1);
    line_x  = xx(PixLim(3):PixLim(4));
    x_avg = mean(line_x);
    [MaxLine,max_ind] = max(Line_minBG);
    SumLine = sum(Line_minBG);
    center = sum(line_x.*Line_minBG)/sum(Line_minBG);
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
    if mod(j,n_samp) == 0
        profiles(:,j/n_samp + 1) = ProfZLi;
        profaxes(:,j/n_samp + 1) = zzLi;
    end
    
    % Calculate residual
    %residual(j) = sum(Line_minBG.*(ProfXLi' - Line_minBG).^2);
    residual(j+1) = sum((ProfXLi' - Line_minBG).^2);

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
    
    % Update interesting shit
    interesting_shit(1,j+1) = datenum(clock);                   % Step time
    interesting_shit(2,j+1) = FWHM(dd,Line_minBG);              % Measured FWHM
    interesting_shit(3,j+1) = sqrt(sum((dd.^2)...               % Measured RMS
        .*Line_minBG)/sum(Line_minBG));
    interesting_shit(4,j+1) = del;                              % Measured Energy offset
    interesting_shit(5,j+1) = OUT.E.FWHM(nOut);                 % Calc'd FWHM
    interesting_shit(6,j+1) = sqrt(sum((OUT.E.AXIS(:,nOut)...   % Calc'd RMS
        .^2).*OUT.E.HIST(:,nOut))/sum(OUT.E.HIST(:,nOut)));     
    interesting_shit(7,j+1) = OUT.Z.FWHM(nOut)/2.35;            % Calc's sigma_z
    
    machine_info = lcaGetSmart({'OTRS:LI20:3075:BLEN';'DR13:TORO:40:DATA';'DR13:AMPL:11:VACT';'SIOC:SYS1:ML00:AO854';'SIOC:SYS1:ML00:CALC029';'BLEN:LI20:3014:BRAW';});
    interesting_shit(8,j+1)  = machine_info(1);                 % TCAV BLEN
    interesting_shit(9,j+1)  = machine_info(2);                 % NRTL TORO
    interesting_shit(10,j+1) = machine_info(3);                 % NRTL AMPL
    interesting_shit(11,j+1) = machine_info(4);                 % NDR BLM
    interesting_shit(12,j+1) = machine_info(5);                 % NRTL PHAS
    interesting_shit(13,j+1) = machine_info(6);                 % S20 PYRO
    interesting_shit(14,j+1) = im_dat.pulseId;                  % Pulse ID
    
    % Record 1 in 10 profiles and spectra
    if mod(j-9,n_samp) == 0
        profiles(:,(j-9)/n_samp + 1) = ProfZLi;
        profaxes(:,(j-9)/n_samp + 1) = zzLi;
        ll = length(Line_minBG);
        yagspecs(1:ll,(j-9)/n_samp + 1) = Line_minBG;
        Li_specs(1:ll,(j-9)/n_samp + 1) = ProfXLi;
        specaxes(1:ll,(j-9)/n_samp + 1) = line_x;
        index((j-9)/n_samp + 1) = j;
    end
    
    % Average output over ES period
    if j <= period
        par_out = params(:,j+1);
        par_out(3) = par_out(3)/1e10;
        int_out = interesting_shit(2:7,j+1);
        int_out(6) = 1000*int_out(6);
    else
        par_out = mean(params(:,(j-period):(j+1)),2);
        par_out(3) = par_out(3)/1e10;
        int_out = mean(interesting_shit(2:7,(j-period):(j+1)),2);
        int_out(6) = 1000*int_out(6);
    end
        
    
    mine = lcaGet('SIOC:SYS1:ML00:AO815');
    if mine == ismine
        lcaPutSmart(strcat('SIOC:SYS1:ML00:AO8',{'01','02','03','04','05','06','07','08','09','10','11','12','13','14'})',par_out);
        lcaPutSmart(strcat('SIOC:SYS1:ML00:AO8',{'16','17','18','19','20','21','22','23'})',...
            [int_out;residual(j+1);1-Part_frac(j+1)]);
    end
    
    if show
        
        if j==0
        
            figure(2);
            image(xx,yy,im_dat.img,'CDataMapping','scaled');
            ol_yag = get(gca,'Children');
            axis xy;
            axis image;
            caxis([0 256]);
            rectangle('Position',[LineLim(3),LineLim(1),LineLim(4)-LineLim(3),LineLim(2)-LineLim(1)],...
                'edgecolor','r','linewidth',1,'linestyle','--');
            xlabel('x (mm)');
            ylabel('y (mm)');
            y_ti = title(['YAGS:LI20:2432 ' datestr(clock)]);
        
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
            set(ol_yag,'CData',im_dat.img);
            %title(y_ti, 'title', ['YAGS:LI20:2432 ' datestr(clock)]);
            
            set(spectra(1),'XData',line_x,'YData',Line_minBG);
            set(spectra(2),'XData',line_x,'YData',ProfXLi);
            ax_max =  max(MaxLine,MaxSim)+10;
            new_ax = [LineLim(3) LineLim(4) 0 ax_max];
            axis(h1,new_ax);
            
            set(res_plot,'YData',residual);
            set(part_plot,'YData',1-Part_frac);
            set(prof_plot,'XData',zzLi,'YData',ProfZLi);
            
            
        end
            
            
    end
    
    j = j + 1;
    if j == ESsteps
        display('Saving. . .');
        dataDate = clock;
        saveFull = datestr(dataDate,'yyyy_mm_dd_HHMMSS');
        dataRoot=fullfile(getenv('MATLABDATAFILES'),'data');
        dataYear=datestr(dataDate,'yyyy');
        dataMon=datestr(dataDate,'yyyy-mm');
        dataDay=datestr(dataDate,'yyyy-mm-dd');
        pathName=fullfile(dataRoot,dataYear,dataMon,dataDay);
        if ~exist(pathName,'dir'), try mkdir(pathName);catch end, end
        fileName = ['ESID_' num2str(ismine,'%02d') '_It_' num2str(k,'%01d') '_' saveFull '.mat'];
        
        clear('OUT');
        save(fullfile(pathName,fileName));
        
        params  = zeros(nPar,ESsteps);
        pscaled = zeros(nPar,ESsteps);
        
        YAG_FWHM = zeros(1,ESsteps);
        YAG_RMS  = zeros(1,ESsteps);
        YAG_DELT = zeros(1,ESsteps);
        LiT_FWHM = zeros(1,ESsteps);
        LiT_RMS  = zeros(1,ESsteps);
        LiT_GAUS = zeros(1,ESsteps);
        
        sim_time = zeros(1,ESsteps);

        cost      = zeros(1,ESsteps);   % Cost
        Part_frac = zeros(1,ESsteps);   % Fraction of Particles lost
        residual  = zeros(1,ESsteps);   % Chi2 difference between spectra
        
        k = k + 1;
        j = 0;
        
        
    end
    
end
