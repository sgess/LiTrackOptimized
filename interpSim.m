function SIM_SPEC = interpSim(OUT,spectrum_axis,nBin)

% Identify Max and Min of Simulated energy distribution
x_max = OUT.X.AXIS(nBin);
x_min = OUT.X.AXIS(1);

% Find the Max and Min on the YAG energy axis
[~,iMax] = min(abs(x_max - spectrum_axis));
[~,iMin] = min(abs(x_min - spectrum_axis));
N = iMax - iMin + 1;

% Interpolate the simulated distribution onto the YAG axis
xx = linspace(1,nBin,N);
SX = interp1(OUT.X.HIST,xx);
sim_sum = sum(SX);
sim_cen = round(sum((1:N).*SX)/sim_sum);

% Embed onto spectrum axis
PIX = length(spectrum_axis);
SIM_SPEC = zeros(PIX,1);
SIM_SPEC(round(PIX/2-sim_cen):round(PIX/2-sim_cen+N-1)) = SX/sim_sum;