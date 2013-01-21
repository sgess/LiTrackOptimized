npulses = nprompt('Number of pulses per scanned variable',3);            % ODD number of pulses to run

npulstart_time = get_time;
if rem(npulses,2)==0
  error('npulses must be an odd number, e.g. 5')
end

seed  = 1;              % does nothing
c     = 2.99792458E8;

fn    = input('Litrack file name (e.g. lcls): ','s');
fnf = [fn '_lit'];
if exist(fnf)~=2
  error([fnf ' LiTrack file does not exist'])
end

eval(fnf);
beamline0 = beamline;

[zpos0,dE_E0,E_bar0,dFWpct0,ZFWmm0,z_bar0,E_barcuts0,fcut0] = ...
                                  litrack(fn,seed,0,0,0,beamline0);
nlocs = length(E_bar0);             % number of locations along beamline where bunch length measurements are made
if nlocs ~= 4
  error('Need 4 minus-sign locations in LiTrack beamline for this scan')
end

sigzj      = zeros(npulses,nlocs);
sigdj      = zeros(npulses,nlocs);
E_barj     = zeros(npulses,nlocs);
dFWpctj    = zeros(npulses,nlocs);
ZFWmmj     = zeros(npulses,nlocs);
z_barj     = zeros(npulses,nlocs);
E_barcutsj = zeros(npulses,nlocs);
fcutj      = zeros(npulses,nlocs);

j = (npulses+1)/2;
for k = 1:nlocs
  iok = find(zpos0(:,k)~=0);
  sigzj(j,k) = std(zpos0(iok,k));
  sigdj(j,k) = std(dE_E0(iok,k));
end
E_barj(j,:)     = E_bar0;
dFWpctj(j,:)    = dFWpct0;
ZFWmmj(j,:)     = ZFWmm0;
z_barj(j,:)     = z_bar0;
E_barcutsj(j,:) = E_barcuts0;
fcutj(j,:)      = fcut0;

% [1] row, [2] column, [3] rel(=0) or abs(=1), [4] rms variation, [5] where in beamline minus-signs to observe E and sig_z, and
% [6] is scalar to get blfunt (units) right
blf( 1,:) = [1  2  0  0.0010  1  100];	% L0 relative voltage
blf( 2,:) = [1  3  1  0.1000  1    1];	% L0 phase
blf( 3,:) = [3  2  0  0.0010  2  100];	% L1 relative voltage
blf( 4,:) = [3  3  1  0.1000  2    1];	% L1 phase
blf( 5,:) = [4  2  0  0.0025  2  100];	% Lx relative voltage
blf( 6,:) = [4  3  1  0.3000  2    1];	% Lx phase
blf( 7,:) = [6  2  0  0.0007  3  100];	% L2 relative voltage
blf( 8,:) = [6  3  1  0.0700  3    1];	% L2 phase
blf( 9,:) = [9  2  0  0.0005  4  100];	% L3 relative voltage
blf(10,:) = [9  3  1  0.0700  4    1];	% L3 phase
blfstr = ['L0 Relative Voltage'
          'L0 Phase           '
          'L1 Relative Voltage'
          'L1 Phase           '
          'Lx Relative Voltage'
          'Lx Phase           '
          'L2 Relative Voltage'
          'L2 Phase           '
          'L3 Relative Voltage'
          'L3 Phase           '];
blfunt = ['%  '
          'deg'
          '%  '
          'deg'
          '%  '
          'deg'
          '%  '
          'deg'
          '%  '
          'deg'];
      
[nk,dum] = size(blf);

s = ((-(npulses-1)/2):(npulses-1)/2);
for kk = 1:nk;
  disp(' ')
  disp(['Scanning ' blfstr(kk,:) '...'])
  delta = zeros(nk,1);
  delta(kk) = -(npulses-1)/2;
  for j = 1:npulses
    disp(['Running pulse # ' int2str(j)])
    if j~=(npulses+1)/2
      beamline_new = litrack_beamline(blf,delta,beamline);
      [zpos,dE_E,E_bar,dFWpct,ZFWmm,z_bar,E_barcuts,fcut] = ...
                                        litrack(fn,seed,0,0,0,beamline_new);
      for k = 1:nlocs
        iok = find(zpos(:,k)~=0);
        sigzj(j,k) = std(zpos(iok,k));
        sigdj(j,k) = std(dE_E(iok,k));
      end
      E_barj(j,:)     = E_bar;
      dFWpctj(j,:)    = dFWpct;
      ZFWmmj(j,:)     = ZFWmm;
      z_barj(j,:)     = z_bar;
      E_barcutsj(j,:) = E_barcuts;
      fcutj(j,:)      = fcut;
    end
    delta(kk)       = delta(kk) + 1;
  end
  figure(kk)
  subplot(221)
  plot_polyfit(s*blf(kk,4)*blf(kk,6),E_barj(:,blf(kk,5)),1,2,blfstr(kk,:),'{\itE}',blfunt(kk,:),'GeV')
  enhance_plot('times',12)
  subplot(222)
  plot_polyfit(s*blf(kk,4)*blf(kk,6),sigzj(:,blf(kk,5))*1E3,1,2,blfstr(kk,:),'{\it\sigma_z}',blfunt(kk,:),'mm')
  enhance_plot('times',12)
  subplot(223)
  plot_polyfit(s*blf(kk,4)*blf(kk,6),E_barj(:,nlocs),1,2,blfstr(kk,:),'{\itE} in Undulator',blfunt(kk,:),'GeV')
  enhance_plot('times',12)
  subplot(224)
  plot_polyfit(s*blf(kk,4)*blf(kk,6),sigzj(:,nlocs)*1E3,1,2,blfstr(kk,:),'{\it\sigma_z} in Undulator',blfunt(kk,:),'mm')
  enhance_plot('times',12)
  stop_time = get_time;
  save feedback_litrack.mat
end
