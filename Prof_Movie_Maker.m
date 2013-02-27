% Number of spectrums to plot
nspec=100;

% Chunks to break up total steps into, to plot total of nspec
chn=floor(ESsteps/nspec);

% ES Step Number
esn=zeros(1,nspec);

% Fast stuff?
global A;
A = load('slac.dat');


% Spectrum Movie
writerObj = VideoWriter('Beam_MORPH2_Feb_9.avi'); % Name it.
writerObj.FrameRate = 30; % How many frames per second.
open(writerObj); 

% Create a spectrum/cost/parameters at each step
for jsm=1:nspec;
    
    jsm
    
    % ES step number
    esn(jsm)=1+(jsm-1)*chn;
    

    % Set the paraeters to those values
    PARAM.INIT.SIGZ0 = params(1,1+(jsm-1)*chn);           % Bunch Length
    PARAM.INIT.SIGD0 = params(2,1+(jsm-1)*chn);           % Initial Energy Spread
    PARAM.INIT.NPART = params(3,1+(jsm-1)*chn);           % Number of Particles
    PARAM.INIT.ASYM = params(4,1+(jsm-1)*chn);            % Initial Gaussian Asymmetry
    PARAM.NRTL.AMPL = params(5,1+(jsm-1)*chn);            % Amplitude of RF Compressor
    PARAM.NRTL.PHAS = params(6,1+(jsm-1)*chn);            % RF Compressor Phase
    PARAM.NRTL.ELO = params(7,1+(jsm-1)*chn);             % Low Energy Cutoff
    PARAM.NRTL.EHI = params(8,1+(jsm-1)*chn);             % High Energy Cutoff
    decker = params(9,1+(jsm-1)*chn);                     % 2-10 Phase
    ramp = params(10,1+(jsm-1)*chn);                      % Ramp Phase
    PARAM.LI10.ELO = params(11,1+(jsm-1)*chn);            % Low Energy Cutoff
    PARAM.LI10.EHI = params(12,1+(jsm-1)*chn);            % High Energy Cutoff
    PARAM.LI20.ELO = params(13,1+(jsm-1)*chn);            % Low Energy Cutoff
    PARAM.LI20.EHI = params(14,1+(jsm-1)*chn);            % High Energy Cutoff
    PARAM.LI20.R56 = params(15,1+(jsm-1)*chn);
    PARAM.LI20.NLO = params(16,1+(jsm-1)*chn);
    PARAM.LI20.NHI = params(17,1+(jsm-1)*chn);
    dZ = params(18,1+(jsm-1)*chn);
    
    
        
    PARAM.LONE.PHAS = decker+ramp;  % Total Phase
    PARAM.LONE.GAIN = (PARAM.ENRG.E1 - PARAM.ENRG.E0)/cosd(PARAM.LONE.PHAS); % Energy gain
    PARAM.LI20.T566  = p(1)*PARAM.LI20.R56^2 + p(2)*PARAM.LI20.R56 + p(3);
    
    % Run LiTrack
    %OUT = LiTrackOpt('FACETpar');
    OUT = LiTrackOpt('FACETconNOTCH');
    Part_frac(jsm) = 1 - OUT.I.PART(2)/PARAM.INIT.NESIM;
    
    
    % Interpolate simulated spectrum
    SIM = interpSimProf(OUT,PARAM.SIMU.BIN,2,AXIS,dZ);
    
    % Calculate residual
    %residual(j) = sum((sim_spectrum - data_spectrum).^2);
    %residual(j) = sum(PROF.*(SIM - PROF).^2);
    residual(jsm) = sum((SIM - PROF).^2);
    
    % Set Cost as the value of the residual
    cost(jsm) = residual(jsm);
    %cost(j) = 14 + log(residual(j)) + 0.001*Part_frac(j);
    %cost(j) = 20 + log(residual(j)) + 0.001*Part_frac(j);
    %cost(j) = 20 + log(residual(j));
    
    
    % Check the size of PROF, and cut off profile to have the same length
    dprof=zeros(1,length(PROF));
    dprof(1,:)=profile(1,1:length(PROF));
    
    
    % Morph into the new shape over cs time steps, then stop morphing
    if 1+(jsm-1)*chn>cs;
        if j<2*cs+1;
            j3=j3+1+(jsm-1)*chn;
        end
    end
    
    % Slowly transform back into the original shape, to try some
    % "trajectory tracking"
    if 1+(jsm-1)*chn>3*cs;
        if j3>1.5;
            j3=j3-0.5*(1+(jsm-1)*chn);
        end
    end
    
    
    PROF=((cs-j3)/(cs-1))*PROF+((j3-1)/cs)*dprof;
    
    
    
    
    fid=figure;
    subplot(2,2,1);
    plot(1000*AXIS,PROF,'g',1000*AXIS,SIM,'b','linewidth',2);
    axis([1000*pMin 1000*pMax 0 0.04]);
    xlabel('Z (\mum)','fontsize',12);
    title('Beam Profiles','fontsize',10);
    legend('MODEL','SIM');
    
    subplot(2,2,2);
    plot(1:1+(jsm-1)*chn,residual(1:1+(jsm-1)*chn),'color','r','linewidth',2);
    %axis([0 ESsteps 0 7e-6]);
    xlabel('Step','fontsize',12);
    title('Residual','fontsize',10);
    
    subplot(2,2,3);
    plot(1:1+(jsm-1)*chn,1-Part_frac(1:1+(jsm-1)*chn),'color','g','linewidth',2);
    %axis([0 ESsteps 0.75 1.1]);
    xlabel('Step','fontsize',12);
    title('Particle Fraction','fontsize',10);
    
    subplot(2,2,4);
    plot(1:1+(jsm-1)*chn,cost(1:1+(jsm-1)*chn),'color','b','linewidth',2);
    %axis([0 ESsteps 0 9]);
    xlabel('Step','fontsize',12);
    title('Cost','fontsize',10);
    
    figure(fid); % Makes sure you use your desired frame.
    
    
    %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end
 
    close all
    
end

hold off
close(writerObj); % Saves the movie.


