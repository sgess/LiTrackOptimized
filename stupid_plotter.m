% Plot all parameters and cost

figure(2)
for jp1=1:2;
    for jp2=1:9;
        subplot(9,2,jp2+(jp1-1)*9);
            plot(pscaled(jp2+(jp1-1)*9,1:j-1))
            %hold on
            %plot(ave(jp2+(jp1-1)*9,1:j-1),'r')
    end
end

% Plot slope of all parameters and cost

%figure(3)
%for jp1=1:2;
%    for jp2=1:9;
%        subplot(9,2,jp2+(jp1-1)*9);
%            plot(aved(jp2+(jp1-1)*9,120:j-1))
%    end
%end

% Plot all scaled parameters together

figure(4)
plot(cost(jbest,1:j-1))
%hold on
%plot(avec(1:j-1),'r')

% Plot Cost

figure(5)
plot(pscaled(:,1:j-1)')


% Find the best setting of the parameters
[y,mj] = min(cost(1,2:j-1));

aveparams=zeros(1,18);

aveparams(1,:)=params(:,mj);


    PARAM.INIT.SIGZ0 = aveparams(1);           % Bunch Length
    PARAM.INIT.SIGD0 = aveparams(2);           % Initial Energy Spread
    PARAM.INIT.NPART = aveparams(3);           % Number of Particles
    PARAM.INIT.ASYM = aveparams(4);            % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL = aveparams(5);            % Amplitude of RF Compressor
    PARAM.NRTL.PHAS = aveparams(6);            % RF Compressor Phase
    PARAM.NRTL.ELO = aveparams(7);             % Low Energy Cutoff
    PARAM.NRTL.EHI = aveparams(8);             % High Energy Cutoff
    decker = aveparams(9);                     % 2-10 Phase
    ramp = aveparams(10);                      % Ramp Phase
    PARAM.LI10.ELO = aveparams(11);            % Low Energy Cutoff
    PARAM.LI10.EHI = aveparams(12);            % High Energy Cutoff
    PARAM.LI20.ELO = aveparams(13);            % Low Energy Cutoff
    PARAM.LI20.EHI = aveparams(14);            % High Energy Cutoff
    PARAM.LI20.BETA = aveparams(15);           % Beta Function
    PARAM.LI20.R16 = aveparams(16);            % Dispersion
    PARAM.LI20.T166 = aveparams(17);           % 2nd Order Dispersion
    delta = aveparams(18);                     % Energy offset

    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    
    % Interpolate simulated spectrum
    sim_spectrum = interpSim(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);
    
    % Calculate residual
    residual = sum((sim_spectrum - data_spectrum).^2);
    
    % Set Cost as the value of the residual
    cost(jbest,j+1) = residual;
    

% Compare to the initial conditions

    PARAM.INIT.SIGZ0 = paramsorig(1);           % Bunch Length
    PARAM.INIT.SIGD0 = paramsorig(2);           % Initial Energy Spread
    PARAM.INIT.NPART = paramsorig(3);           % Number of Particles
    PARAM.INIT.ASYM = paramsorig(4);            % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL = paramsorig(5);            % Amplitude of RF Compressor
    PARAM.NRTL.PHAS = paramsorig(6);            % RF Compressor Phase
    PARAM.NRTL.ELO = paramsorig(7);             % Low Energy Cutoff
    PARAM.NRTL.EHI = paramsorig(8);             % High Energy Cutoff
    decker = paramsorig(9);                     % 2-10 Phase
    ramp = paramsorig(10);                      % Ramp Phase
    PARAM.LI10.ELO = paramsorig(11);            % Low Energy Cutoff
    PARAM.LI10.EHI = paramsorig(12);            % High Energy Cutoff
    PARAM.LI20.ELO = paramsorig(13);            % Low Energy Cutoff
    PARAM.LI20.EHI = paramsorig(14);            % High Energy Cutoff
    PARAM.LI20.BETA = paramsorig(15);           % Beta Function
    PARAM.LI20.R16 = paramsorig(16);            % Dispersion
    PARAM.LI20.T166 = paramsorig(17);           % 2nd Order Dispersion
    delta = paramsorig(18);                     % Energy offset

    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain

    
    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    
    % Interpolate simulated spectrum
    sim_spectrum_original = interpSim(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);
    
    % Calculate residual
    residual_original = sum((sim_spectrum_original - data_spectrum).^2);
    
   
    
% Plot Output
figure(6)
plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum_original,'r',spectrum_axis,sim_spectrum,'b');
legend('DATA','INITIAL CONDITIONS','ES-OPTIMIZED SIMULATION');
xlabel('X (mm)','fontsize',14);
text(-3.5,5.5e-3,['Original Residual = ' num2str(residual_original,'%0.2e')],'fontsize',14);
text(-3.5,5e-3,['Residual after ES = ' num2str(residual,'%0.2e')],'fontsize',14);



figure(9)
for jpc=1:bsteps;
    subplot(bsteps,1,jpc)
    plot(cost(jpc,1:ESsteps-1))
end


