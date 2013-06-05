function SetPars(pars, names, nPar)

global PARAM;

REGIONS = {'SIMU';
           'ENRG';
           'INIT';
           'NRTL';
           'LONE';
           'LI10';
           'LTWO';
           'LI20'};
       
PARAMETERS = {'PLOT';
              'FRAC';
              'BIN';
              'ZFIT';
              'DFIT';
              'CONT';
              'E0';
              'E1';
              'E2';
              'SIGZ0';
              'SIGD0';
              'Z0BAR';
              'D0BAR';
              'NESIM';
              'NPART';
              'ASYM';
              'TAIL';
              'CUT';
              'AMPL';
              'PHAS';
              'RAMP';
              'DECK';
              'LEFF';
              'R56';
              'T566';
              'ELO';
              'EHI';
              'FBAM';
              'GAIN';
              'ISR';
              'NLO';
              'NHI';
              'R16';
              'T166';
              'BETA';
              'EMIT'};

remain = names;
for j = 1:2
    
    [FIELD(:,j),remain] = strtok(remain);
    
end
    
for i = 1:nPar
    reg = strcmp(FIELD(i,1),REGIONS);
    reg_exist = sum(reg);
    if ~reg_exist
        error(['Invalid region "' FIELD(i,1) '" for parameter ' num2str(i) ' with name "' names(i) '"']);
    end
    
    par = strcmp(FIELD(i,2),PARAMETERS);
    par_exist = sum(par);
    if ~par_exist
        error(['Invalid region "' FIELD(i,2) '" for parameter ' num2str(i) ' with name "' names(i) '"']);
    else
        PARAM.(char(FIELD(i,1))).(char(FIELD(i,2))) = pars(i); 
    end
end

PARAM.LONE.PHAS = PARAM.LONE.RAMP+PARAM.LONE.DECK;
PARAM.LTWO.PHAS = PARAM.LONE.RAMP;
PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
PARAM.LTWO.GAIN = (PARAM.ENRG.E2 - PARAM.ENRG.E1)/cosd(PARAM.LTWO.PHAS); % Energy gain