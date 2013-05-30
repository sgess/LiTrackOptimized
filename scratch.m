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
                