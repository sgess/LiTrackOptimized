% Number of spectrums to plot
nspec=500;

% The spectrums
simspectrum_holder=zeros(560,nspec);
cost_holder=zeros(1,nspec);
params_holder=zeros(18,nspec);

% Chunks to break up total steps into, to plot total of nspec
chn=floor(ESsteps/nspec);

% ES Step Number
esn=zeros(1,nspec);

% Fast stuff?
global A;
A = load('slac.dat');


% Create a spectrum/cost/parameters at each step
for jsm=1:nspec;
    
    jsm
    
    % ES step number
    esn(jsm)=1+(jsm-1)*chn;
    
    
    % Record every nspecth parameter
    params_holder(:,jsm)=params(:,2+(jsm-1)*chn);
    
    
    % Set the paraeters to those values
    PARAM.INIT.SIGZ0 = params_holder(1,jsm);                    % Bunch Length
    PARAM.INIT.SIGD0 = params_holder(2,jsm);                    % Initial Energy Spread
    PARAM.INIT.NPART = params_holder(3,jsm);                    % Number of Particles
    PARAM.INIT.ASYM = params_holder(4,jsm);                     % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL = params_holder(5,jsm);                     % Amplitude of RF Compressor
    PARAM.NRTL.PHAS = params_holder(6,jsm);                     % RF Compressor Phase
    PARAM.NRTL.ELO = params_holder(7,jsm);                      % Low Energy Cutoff
    PARAM.NRTL.EHI = params_holder(8,jsm);                      % High Energy Cutoff
    decker = params_holder(9,jsm);                              % 2-10 Phase
    ramp = params_holder(10,jsm);                               % Ramp Phase
    PARAM.LI10.ELO = params_holder(11,jsm);                     % Low Energy Cutoff
    PARAM.LI10.EHI = params_holder(12,jsm);                     % High Energy Cutoff
    PARAM.LI20.ELO = params_holder(13,jsm);                     % Low Energy Cutoff
    PARAM.LI20.EHI = params_holder(14,jsm);                     % High Energy Cutoff
    PARAM.LI20.BETA = params_holder(15,jsm);                    % Beta Function
    PARAM.LI20.R16 = params_holder(16,jsm);                     % Dispersion
    PARAM.LI20.T166 = params_holder(17,jsm);                    % 2nd Order Dispersion
    delta = params_holder(18,jsm);                              % Energy offset

    
    % Random crazy Spencer stuff
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain


    % Run LiTrack
    OUT = LiTrackOpt('FACETpar');
    
    % Interpolate simulated spectrum
    simspectrum_holder(:,jsm) = interpSim(OUT,spectrum_axis,PARAM.SIMU.BIN,delta,PARAM.LI20.R16);
    
    % Calculate residual
    residual = sum((simspectrum_holder(:,jsm) - data_spectrum).^2);
    
    % Set Cost as a function of the residual
    cost_holder(jsm) = residual + eps*Part_frac(1+(jsm-1)*chn);
    
    
    % Create the movie frames
    
    %fid = figure;
    %hold on
    
    % Spectrums
    %plot(spectrum_axis,data_spectrum,'g',spectrum_axis,simspectrum_holder(:,1),'r',spectrum_axis,simspectrum_holder(:,jsm),'b');
    %legend('DATA','INITIAL CONDITIONS','ES-SPECTRUM');
    %xlabel('X (mm)','fontsize',14);

    %hold off
    
    %fid2 = figure;
    %hold on
    % Parameters
    %subplot(2,1,1)
    %%title('Scaled Parameters','fontsize',14)
    %xlabel('ES Step Number','fontsize',14);
    
    % Cost
    %subplot(2,1,2)
    %plot(cost_holder(1:jsm));
    
    
    
    
    
end
    


%% Set up the movies.

% Spectrum Movie
writerObj = VideoWriter('spect_Feb_5.avi'); % Name it.
writerObj.FrameRate = 30; % How many frames per second.
open(writerObj); 
 
for i=1:nspec;      
    % We just use pause but pretend you have some really complicated thing here...
    i
    % Spectrums
    fid=figure;
    plot(spectrum_axis,data_spectrum,'g',spectrum_axis,simspectrum_holder(:,1),'r',spectrum_axis,simspectrum_holder(:,i),'b');
    legend('DATA','INITIAL CONDITIONS','ES-SPECTRUM');
    xlabel('X (mm)','fontsize',14);
    
    %pause(0.1);
    figure(fid); % Makes sure you use your desired frame.
    
    
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 
    close all
end
hold off
close(writerObj); % Saves the movie.
