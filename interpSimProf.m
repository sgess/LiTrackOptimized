function PROF_SPEC = interpSimProf(OUT,samp_axis,nBin,dZ)

% Identify Max and Min of Simulated energy distribution
z_max = OUT.Z.AXIS(nBin,2);
z_min = OUT.Z.AXIS(1,2);

% Find the Max and Min on the YAG energy axis
[~,iMax] = min(abs(z_max - samp_axis));
[~,iMin] = min(abs(z_min - samp_axis));
N = iMax - iMin + 1;

% Interpolate the simulated distribution onto the YAG axis
zz = linspace(1,nBin,N);
SZ = interp1(OUT.Z.HIST(:,2),zz);
sim_sum = sum(SZ);
sim_cen = round(sum((1:N).*SZ)/sim_sum);

% Embed onto spectrum axis
PIX = length(samp_axis);
PROF_SPEC = zeros(PIX,1);
bin = samp_axis(2) - samp_axis(1);
offset = round(dZ/bin);
max_bin = offset + PIX/2 - sim_cen + N - 1;
min_bin = offset + PIX/2 - sim_cen;
if max_bin > PIX
    warning('dZ = %0.4f is to high and moves profile out of range.',dZ);
    top_bin = N - (max_bin - PIX) + 1;
    min_bin = max_bin - N;
    PROF_SPEC(min_bin:PIX) = SZ(1:top_bin)/sim_sum;
elseif min_bin < 1
    warning('dZ = %0.4f is to low and moves profile out of range.',dZ);
    bot_bin = 1 - min_bin;
    max_bin = N - bot_bin + 1;
    PROF_SPEC(1:max_bin) = SZ(bot_bin:N)/sim_sum;
else
    PROF_SPEC(min_bin:max_bin) = SZ/sim_sum;
end
