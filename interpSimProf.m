function [PROF_SPEC, PROF_AXIS] = interpSimProf(OUT,nBin,iEL,dAX,pMin,pMax,dZ)

% Create Profile axis
N_lo = floor(pMin/dAX);
N_hi = ceil(pMax/dAX);
PROF_AXIS = dAX*(N_lo:N_hi);
Bins = length(PROF_AXIS);

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
offset = round(dZ/dAX);
max_bin = N + offset;
min_bin = offset;
PROF_SPEC = zeros(1:Bins);
if max_bin > N_hi
    warning('dZ = %0.4f is to high and moves profile out of range.',dZ);
    top_bin = max_bin - N_hi + 1;
    min_bin = max_bin - N;
    PROF_SPEC(min_bin:) = SZ(1:top_bin)/sim_sum;
elseif min_bin < N_lo
    warning('dZ = %0.4f is to low and moves profile out of range.',dZ);
    bot_bin = 1 - min_bin;
    max_bin = N - bot_bin + 1;
    PROF_SPEC(1:max_bin) = SZ(bot_bin:N)/sim_sum;
else
    PROF_SPEC(min_bin:max_bin) = SZ/sim_sum;
end
