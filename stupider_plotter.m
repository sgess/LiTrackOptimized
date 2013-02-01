res = residual(1:j-1);
cst = cost(1:j-1);

figure(2)
for jp1=1:2;
    for jp2=1:9;
        subplot(9,2,jp2+(jp1-1)*9);
            plot(params(jp2+(jp1-1)*9,1:j-1));
            title(par_name(jp2+(jp1-1)*9));
    end
end



[a,b]=min(res);
[c,d]=min(cst);

k = b;


PARAM.INIT.SIGZ0 = params(1,k);           % Bunch Length
PARAM.INIT.SIGD0 = params(2,k);           % Initial Energy Spread
PARAM.INIT.NPART = params(3,k);           % Number of Particles
PARAM.INIT.ASYM = params(4,k);            % Initial Gaussian Asymmetry
PARAM.NRTL.AMPL = params(5,k);            % Amplitude of RF Compressor
PARAM.NRTL.PHAS = params(6,k);            % RF Compressor Phase
PARAM.NRTL.ELO = params(7,k);             % Low Energy Cutoff
PARAM.NRTL.EHI = params(8,k);             % High Energy Cutoff
decker = params(9,k);                     % 2-10 Phase
ramp = params(10,k);                      % Ramp Phase
PARAM.LI10.ELO = params(11,k);            % Low Energy Cutoff
PARAM.LI10.EHI = params(12,k);            % High Energy Cutoff
PARAM.LI20.ELO = params(13,k);            % Low Energy Cutoff
PARAM.LI20.EHI = params(14,k);            % High Energy Cutoff
PARAM.LI20.BETA = params(15,k);           % Beta Function
PARAM.LI20.R16 = params(16,k);            % Dispersion
PARAM.LI20.T166 = params(17,k);           % 2nd Order Dispersion
delta = params(18,k);                     % Energy offset

PARAM.LONE.PHAS = decker+ramp;  % Total Phase
PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain


% Run LiTrack
OUT = LiTrackOpt('FACETpar');
pft = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;


% Interpolate simulated spectrum
sim_spectrum = interpSim(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);

% Calculate residual
rft = sum(data_spectrum.*(sim_spectrum - data_spectrum).^2);

% Set Cost as the value of the residual
%cost(j) = residual;
cft = 14 + log(rft) + 0.0001*pft;

figure(1)
plot(spectrum_axis,data_spectrum,'g',spectrum_axis,sim_spectrum,'b');
legend('DATA','ES-SIMULATION');
xlabel('X (mm)','fontsize',14);
text(-3.5,5e-3,['Residual = ' num2str(rft,'%0.2e')],'fontsize',14);

figure(3)
plot(OUT.Z.AXIS(:,2),OUT.Z.HIST(:,2));

figure(4)
plot(cost(1:j));

figure(5)
plot(residual(1:j));

figure(6)
plot(Part_frac(1:j));