% Number of spectrums to plot
nspec=500;

% Chunks to break up total steps into, to plot total of nspec
chn=floor(ESsteps/nspec);

% ES Step Number
esn=zeros(1,nspec);

% Fast stuff?
global A;
A = load('slac.dat');


% Spectrum Movie
writerObj = VideoWriter('Beam_MORPH_Feb_10.avi'); % Name it.
writerObj.FrameRate = 30; % How many frames per second.
open(writerObj); 

% Create a spectrum/cost/parameters at each step
for jsm=1:nspec;
    
    jsm
    
    % ES step number
    esn(jsm)=1+(jsm-1)*chn;

    
    
    
    fid=figure;
    subplot(2,2,1);
    plot(1000*AXIS,PROF_track(1+(jsm-1)*chn,1:length(AXIS)),'g',1000*AXIS,sim_track(1+(jsm-1)*chn,1:length(AXIS)),'b','linewidth',2);
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


