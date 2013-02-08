pscaled(:,j)=2*(params(:,j)-Cent)./Diff;
    

pscaled(:,j+1)=pscaled(:,j)+dt*cos(w*j*dt+gain*cost(j))*(alpha*w)^0.5;
        if pscaled(k,j+1) < -1;
            pscaled(k,j+1) = -1;
        else if pscaled(k,j+1) > 1;
                pscaled(k,j+1) = 1;
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
    PARAM.LI20.BETA = params(15,j+1);           % Beta Function
    PARAM.LI20.R16 = params(16,j+1);            % Dispersion
    PARAM.LI20.T166 = params(17,j+1);           % 2nd Order Dispersion
    delta = params(18,j+1);                     % Energy offset
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain