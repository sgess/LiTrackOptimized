function PROF_SPEC = interpSimProf(OUT,nBin,iEL,PROF_AXIS,dZ)

%Extract axis info
Bins = length(PROF_AXIS);
dAX = PROF_AXIS(2) - PROF_AXIS(1);

% Identify Max and Min of Simulated energy distribution
z_max = OUT.Z.AXIS(nBin,iEL);
z_min = OUT.Z.AXIS(1,iEL);

% Find bin separation in number of sample bins
N = ceil((z_max - z_min)/dAX);

% Interpolate the simulated distribution onto the sample axis
zz = linspace(1,nBin,N);
SZ = interp1(OUT.Z.HIST(:,iEL),zz);
sim_sum = sum(SZ);

% Embed profile spectrum
[~,off_bin] = min(abs(PROF_AXIS - dZ));
max_bin = N + off_bin - 1;
min_bin = off_bin;
PROF_SPEC = zeros(1,Bins);
if max_bin > Bins
    warning('dZ = %0.4f is too high and moves profile out of range.',dZ);
    top_bin = Bins - off_bin + 1;
    min_bin = off_bin;
    PROF_SPEC(min_bin:Bins) = SZ(1:top_bin)/sim_sum;
elseif off_bin == 1
    warning('dZ = %0.4f is too low and moves profile out of range.',dZ);
    bot_bin = round((PROF_AXIS(1) - dZ)/dAX) + 1;
    max_bin = N - bot_bin + 1;
    PROF_SPEC(1:max_bin) = SZ(bot_bin:N)/sim_sum;
else
    PROF_SPEC(min_bin:max_bin) = SZ/sim_sum;
end
