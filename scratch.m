data.YAG.profile = zeros(130,237);
data.YAG.prof_ax = zeros(130,237);



data.YAG.match_spec = zeros(760,237);
j=1;
for i=1:250
    
    if data.YAG.good_shot(i)
        
        PARAM = data.YAG.param(j);
        OUT = LiTrackOpt('FACETpar');

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
        
        data.YAG.profile(:,i) = ProfZLi;
        data.YAG.prof_ax(:,i) = zzLi;
        data.YAG.match_spec(:,i) = ProfXLi;
        j=j+1;
    end
    
end
%%

data.YAG.prof2 = zeros(260,250);
l=min(min(data.YAG.prof_ax));
h=max(max(data.YAG.prof_ax));
PAX = linspace(l,h,260);
nOut = 3;
j=1;
for i=1:250
    
    if data.YAG.good_shot(i)
        
        PARAM = data.YAG.param(j);
        OUT = LiTrackOpt('FACETpar');
        
        % Get bunch profile
        dZ      = OUT.Z.AXIS(2,nOut) - OUT.Z.AXIS(1,nOut);
        zzLi    = [OUT.Z.AXIS(1,nOut)-dZ; OUT.Z.AXIS(:,nOut); OUT.Z.AXIS(end,nOut)+dZ];
        ProfZLi = [0; OUT.Z.HIST(:,nOut); 0];
        center = sum(zzLi.*ProfZLi)/sum(ProfZLi);
        bl = interp1(data.YAG.prof_ax(:,i)-center,data.YAG.profiles(:,i),PAX,'linear','extrap');
        
        data.YAG.prof2(:,i) = bl;
        j=j+1;
    end
    
end

%%


for i =1:237
sigz0(i) = data.YAG.param(i).INIT.SIGZ0;
sigd0(i) = data.YAG.param(i).INIT.SIGD0;
npart(i) = data.YAG.param(i).INIT.NPART;
asym(i) = data.YAG.param(i).INIT.ASYM;
ampl(i) = data.YAG.param(i).NRTL.AMPL;
phas(i) = data.YAG.param(i).NRTL.PHAS;
nr56(i) = data.YAG.param(i).NRTL.R56;
nt56(i) = data.YAG.param(i).NRTL.T566;
lone(i) = data.YAG.param(i).LONE.PHAS;
ltwo(i) = data.YAG.param(i).LTWO.PHAS;
end

%%

figure(2);
subplot(2,5,1);
hist(sigz0,20);
xlabel('\sigma_z (m)','fontsize',14);
subplot(2,5,2);
hist(sigd0,20);
xlabel('\delta_0','fontsize',14);
subplot(2,5,3);
hist(npart,20);
xlabel('e^-','fontsize',14);
subplot(2,5,4);
hist(asym,20);
xlabel('Asymmetry','fontsize',14);
subplot(2,5,5);
hist(ampl,20);
xlabel('Comp. Ampl. (GV)','fontsize',14);
subplot(2,5,6);
hist(phas,20);
xlabel('Comp. Phas. (deg)','fontsize',14);
subplot(2,5,7);
hist(nr56,20);
xlabel('NTRL R56 (m)','fontsize',14);
subplot(2,5,8);
hist(nt56,20);
xlabel('NTRL T556 (m)','fontsize',14);
subplot(2,5,9);
hist(lone,20);
xlabel('2-10 Phas. (deg)','fontsize',14);
subplot(2,5,10);
hist(ltwo,20);
xlabel('11-20 Phas. (deg)','fontsize',14);

%%
figure(1);
imagesc(pys,100*data.YAG.axis/85,data.YAG.spectra(:,py_ind));
ylabel('\delta (%)','fontsize',14);
xlabel('Pyro','fontsize',14);
figure(2);
imagesc(pys,data.YAG.PAX,data.YAG.prof2(:,py_ind));
ylabel('Z (mm)','fontsize',14);
xlabel('Pyro','fontsize',14);

figure(3);
[is,i_ind] = sort(data.YAG.I_peak);
imagesc(is,data.YAG.PAX,data.YAG.prof2(:,i_ind));
ylabel('Z (mm)','fontsize',14);
xlabel('I_{peak}','fontsize',14);
%%

% assuming x and y are correlated variables . . .
x = randn(10000,1);
y = 2*randn(10000,1);
xx = -5:0.2:5;
yy = -10:0.2:10;

dist = hist2(x,y,xx,yy);
imagesc(xx,yy,dist);
