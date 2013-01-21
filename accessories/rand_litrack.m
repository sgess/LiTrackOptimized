seed  = 1;
fn  = input('Litrack file name (e.g. lcls): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end
gum  = prompt('Use Gaussian, Uniform, or Maximum for peak current calculation: ','gum','u');

eval(fnf)
%beamline(:,1) = abs(beamline(:,1));     % turn off plots and sample markers
c  = 2.99792458E8;
elec = 1.60217733E-19;
n  = nprompt('How many randomized runs',10,1,10000);
clear blf

% SABER (Mar. 1, 2006):
blf( 1,:) = [ 2  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
blf( 2,:) = [ 2  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
blf( 3,:) = [ 5  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
blf( 4,:) = [ 5  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
blf( 5,:) = [ 9  2  0  0.0008];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
blf( 6,:) = [ 9  3  1  0.3000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
dz0rms  = 0.3E-12*c*1E3;		% initial 'timing' jitter (mm)
dQ_Qrms = 0.025;			    % initial relative charge jitter ( )

% LCLSB for 300-as pulses at 100 pC with spoiler (Aug. 13, 2004):
%blf( 1,:) = [ 2  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 5  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 5  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 6  2  0  0.0025];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 6  3  1  0.5000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 8  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 8  3  1  0.0700];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 0.5E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.02;				    % initial relative charge jitter ( )

% MIT Feb. 6, 2004 (1 nC, 25 ps, non-Parmela input - L1 and L2 phases < 0: mit_lit.m):
%blf( 1,:) = [ 2  2  0  0.00002];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.03000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 3  2  0  0.00002];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 3  3  1  0.00600];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 4  2  0  0.00006];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 4  3  1  0.03500];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 6  2  0  0.00005];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 6  3  1  0.02000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 9,:) = [ 8  2  0  0.00090];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf(10,:) = [ 8  3  1  0.50000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 0.3E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.01;				    % initial relative charge jitter ( )

% LCLS28OCT03:
%blf( 1,:) = [ 2  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 5  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 5  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 6  2  0  0.0025];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 6  3  1  0.5000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 8  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 8  3  1  0.0700];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 9,:) = [11  2  0  0.0008];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf(10,:) = [11  3  1  0.1500];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 0.5E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.02;				    % initial relative charge jitter ( )

% LCLS10JUN03:
%blf( 1,:) = [ 2  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 5  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 5  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 6  2  0  0.0025];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 6  3  1  0.3000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 8  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 8  3  1  0.0700];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 9,:) = [12  2  0  0.0008];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf(10,:) = [12  3  1  0.0500];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 0.5E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.02;				    % initial relative charge jitter ( )

% LCLS16JUL03:
%blf( 1,:) = [ 2  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 5  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 5  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 6  2  0  0.0025];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 6  3  1  0.5000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 8  2  0  0.0010];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 8  3  1  0.0700];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 9,:) = [12  2  0  0.0008];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf(10,:) = [12  3  1  0.1500];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 0.8E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.02;				    % initial relative charge jitter ( )

% TESLAXFEL (Jul-Aug. 2003):
%blf( 1,:) = [ 2  2  0  0.0008];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [ 2  3  1  0.0500];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [ 3  2  0  0.0030];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [ 3  3  1  0.0700];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [ 5  2  0  0.0020];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 6,:) = [ 5  3  1  0.0800];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 7,:) = [ 7  2  0  0.0020];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 8,:) = [ 7  3  1  0.1000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 9,:) = [ 9  2  0  0.0009];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf(10,:) = [ 9  3  1  1.0000];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 1.5E-12*c*1E3;		% initial 'timing' jitter (mm)
%dQ_Qrms = 0.10;				    % initial relative charge jitter ( )
%%dz0rms  = 0.5E-12*c*1E3;		% initial 'timing' jitter (mm)
%%dQ_Qrms = 0.02;				    % initial relative charge jitter ( )

% SPPS?:
%blf( 1,:) = [2  2  0  0.0010/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 1,:) = [2  3  1  0.1000/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 2,:) = [5  2  0  0.0010/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 3,:) = [5  3  1  0.1000/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 4,:) = [9  2  0  0.0010/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%blf( 5,:) = [9  3  1  0.1000/2];	% [1] row, [2] column, [3] rel(=0) or abs(=1), and [4] rms variation
%dz0rms  = 2.0E-12*c*1E3;			% initial 'timing' jitter (mm)
%dQ_Qrms = 0.02;			        	% initial relative charge jitter ( )

[nk,dum] = size(blf);

start_time = get_time;
disp(['START TIME ==>> ' start_time])

beamline0 = beamline;
[zpos0,dE_E0,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0,sigzG0] = ...
                     litrack(fn,seed,0,0,0,beamline0);
%[zpos0,dE_E0,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0] = ...
%                     LiTrack(fn,seed,0,0,0,beamline0);
nsamp = length(zpos0(1,:));

E_bar     = zeros(n,nsamp);
dFWpct    = zeros(n,nsamp);
ZFWmm     = zeros(n,nsamp);
z_bar     = zeros(n,nsamp);
E_barcuts = zeros(n,nsamp);
fcut      = zeros(n,nsamp);
sigzG     = zeros(n,nsamp);
sigE      = zeros(n,1);
sigz      = zeros(n,1);
dd        = zeros(n,1);
Ipk       = zeros(n,1);
psav      = zeros(n,nk);
delta     = randn(n,nk);
dz0in     = dz0rms*randn(n,1) + z0_bar*1E3;
dQ_Qin    = dQ_Qrms*randn(n,1);
Ne0       = Ne;
Nbin0     = Nbin;

i0 = find(zpos0(:,end));
Nesim0 = length(i0);
[N,Z] = hist(zpos0(i0,end),Nbin0);
%Nbin1 = round(Nbin0*0.10);
%Nbin2 = round(Nbin0*0.90);
%q0    = plot_parab(Z(Nbin1:Nbin2)*1E6,N(Nbin1:Nbin2),1,'z','N','mm',' ');
ii = find(dE_E0(:,end));
sigE0 = 100*std(dE_E0(ii,end));
ii = find(zpos0(:,end));
%sigz0 = 1E3*std(zpos0(ii,end));
sigz0 = sigzG0(end);				 % gaussian fit used for SPPS
if gum=='m'
  Ipk0 = 1E3*Ne0*elec*c*max(N)/sum(N)*(1-fcut0(end))/mean(diff(Z));
elseif gum=='u'
  Ipk0 = 1E3*Ne0*elec*c/sigz0*(1-fcut0(end))/sqrt(12);
else
  Ipk0 = 1E3*Ne0*elec*c/sigz0*(1-fcut0(end))/sqrt(2*pi);
end
for j = 1:n
  disp(' ')
  disp(sprintf('RUN # %4.0f ...',j))
  disp(' ')
  for k = 1:nk
    p0  = beamline0(blf(k,1),blf(k,2));
    p = p0 + blf(k,3)*blf(k,4)*delta(j,k) + ...
            (1-blf(k,3))*p0*blf(k,4)*delta(j,k);	% add absolute OR relative change
    psav(j,k) = p;
    beamline(blf(k,1),blf(k,2)) = p;
  end
  [zpos,dE_E,E_bar(j,:),dFWpct(j,:),ZFWmm(j,:),z_bar(j,:),E_barcuts(j,:),fcut(j,:),sigzG(j,:)] = ...
                       litrack(fn,seed,dz0in(j),dQ_Qin(j),0,beamline);
%  [zpos,dE_E,E_bar(j,:),dFWpct(j,:),ZFWmm(j,:),z_bar(j,:),E_barcuts(j,:),fcut(j,:)] = ...
%                       LiTrack(fn,seed,dz0in(j),dQ_Qin(j),0,beamline);
  ii = find(dE_E(:,end));
  sigE(j) = 100*std(dE_E(:,end));
  ii = find(zpos(:,end));
%  sigz(j) = 1E3*std(zpos(ii,end));
  sigz(j) = sigzG(j);               % gaussian fit used for SPPS

  dd(j) = 100*(E_barcuts(j,end)/E_barcuts0(end)-1);
  ii0 = find(zpos(:,end));
  [N,Z] = hist(zpos(ii0,end),Nbin0);
%  q = plot_parab(Z(Nbin1:Nbin2)*1E6,N(Nbin1:Nbin2),1,'z','N','mm',' ');
%  Ipk(j) = Ne0*elec*c*max(N)/sum(N)*(1-fcut(j,end))/mean(diff(Z));
  if gum=='m'
    Ipk(j) = 1E3*Ne0*(1+dQ_Qin(j))*elec*c*max(N)/sum(N)*(1-fcut(end))/mean(diff(Z));
  elseif gum=='u'
    Ipk(j) = 1E3*Ne0*(1+dQ_Qin(j))*elec*c/sigz(j)*(1-fcut(end))/sqrt(12);
  else
    Ipk(j) = 1E3*Ne0*(1+dQ_Qin(j))*elec*c/sigz(j)*(1-fcut(end))/sqrt(2*pi);
  end
  if j == 1
    sigz_min = sigz(j);
    sigz_max = sigz(j);
    E_min = E_barcuts(j,end);
    E_max = E_barcuts(j,end);
    beamline_minz = beamline;
    beamline_maxz = beamline;
    beamline_minE = beamline;
    beamline_maxz = beamline;
    dz0in_minz  = dz0in(j);
    dQ_Qin_minz = dQ_Qin(j);
    dz0in_maxz  = dz0in(j);
    dQ_Qin_maxz = dQ_Qin(j);
    dz0in_minE  = dz0in(j);
    dQ_Qin_minE = dQ_Qin(j);
    dz0in_maxE  = dz0in(j);
    dQ_Qin_maxE = dQ_Qin(j);
  end
  if sigz(j) < sigz_min
    beamline_minz = beamline;
    dz0in_minz  = dz0in(j);
    dQ_Qin_minz = dQ_Qin(j);
    sigz_min = sigz(j);
  end
  if sigz(j) > sigz_max
    beamline_maxz = beamline;
    dz0in_maxz  = dz0in(j);
    dQ_Qin_maxz = dQ_Qin(j);
    sigz_max = sigz(j);
  end
  if E_barcuts(j,end) < E_min
    beamline_minE = beamline;
    dz0in_minE  = dz0in(j);
    dQ_Qin_minE = dQ_Qin(j);
    E_min = E_barcuts(j,end);
  end
  if E_barcuts(j,end) > E_max
    beamline_maxE = beamline;
    dz0in_maxE  = dz0in(j);
    dQ_Qin_maxE = dQ_Qin(j);
    E_max = E_barcuts(j,end);
  end
  stop_time = get_time;
  save rand_litrack.mat
end

stop_time = get_time;
disp(['STOP TIME ==>> ' stop_time])
