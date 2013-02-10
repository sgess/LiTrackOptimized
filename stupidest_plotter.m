% Initial conditions
figure(3);

% Initial bunch length
p = 1;
scale = 1000;
subplot(2,2,1);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('\sigma_z (mm)','fontsize',14);
title(par_name(p),'fontsize',14);

% Initial Energy Spread
p = 2;
scale = 100;
subplot(2,2,2);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('\sigma_{\delta} (%)','fontsize',14);
title(par_name(p),'fontsize',14);

% Initial Number of particles
p = 3;
scale = 1e-10;
subplot(2,2,3);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Particles (10^{10})','fontsize',14);
title(par_name(p),'fontsize',14);

% Initial asymmetry
p = 4;
scale = 1;
subplot(2,2,4);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Asymmetry','fontsize',14);
title(par_name(p),'fontsize',14);








% NRTL conditions
figure(4);

% NRTL compressor amplitude
p = 5;
scale = 1000;
subplot(2,2,1);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Amplitude (MV)','fontsize',14);
title(par_name(p),'fontsize',14);

% NRTL compressor phase
p = 6;
scale = 1;
subplot(2,2,2);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Phase (Deg)','fontsize',14);
title(par_name(p),'fontsize',14);

% NRTL R56
p = 7;
scale = 1;
subplot(2,2,3);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('R_{56} (m)','fontsize',14);
title(par_name(p),'fontsize',14);

% NRTL T566
p = 8;
scale = 1;
subplot(2,2,4);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('T_{566} (m)','fontsize',14);
title(par_name(p),'fontsize',14);





% Linac phase conditions
figure(5);

% Decker phase
p = 9;
scale = 1;
subplot(2,2,1);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('2-10 Phase (Deg)','fontsize',14);
title(par_name(p),'fontsize',14);

% Ramp phase
p = 10;
scale = 1;
subplot(2,2,2);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Phase Ramp (Deg)','fontsize',14);
title(par_name(p),'fontsize',14);

% S20 Low
p = 11;
scale = 1;
subplot(2,2,3);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Energy Cut (%)','fontsize',14);
title(par_name(p),'fontsize',14);

% S20 High
p = 12;
scale = 1;
subplot(2,2,4);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Energy Cut (%)','fontsize',14);
title(par_name(p),'fontsize',14);






% YAG map
figure(6);

% Beta at YAG
p = 13;
scale = 1;
subplot(2,2,1);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('\beta (m)','fontsize',14);
title(par_name(p),'fontsize',14);

% Eta at YAG
p = 14;
scale = 1;
subplot(2,2,2);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('\eta (mm)','fontsize',14);
title(par_name(p),'fontsize',14);

% R166
p = 15;
scale = 1;
subplot(2,2,3);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('R166 (mm)','fontsize',14);
title(par_name(p),'fontsize',14);

% delta
p = 16;
scale = 100;
subplot(2,2,4);
plot(scale*params(p,:));
line([0 ESsteps],[scale*Low_limit(p) scale*Low_limit(p)],'color','r');
line([0 ESsteps],[scale*High_limit(p) scale*High_limit(p)],'color','r');
if Low_limit(p) < 0
    lo_ax = 1.1*Low_limit(p);
else
    lo_ax = 0.9*Low_limit(p);
end
if High_limit(p) < 0
    hi_ax = 0.9*High_limit(p);
else
    hi_ax = 1.1*High_limit(p);
end
axis([0 ESsteps scale*lo_ax scale*hi_ax]);
xlabel('Simulation Step','fontsize',14);
ylabel('Offset (%)','fontsize',14);
title(par_name(p),'fontsize',14);