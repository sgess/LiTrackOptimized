function SIM_SPEC = interpSimSpec(OUT,spectrum_axis,nBin,delta,eta)

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
bin = spectrum_axis(2) - spectrum_axis(1);
offset = round(delta*eta/bin);
max_bin = offset + PIX/2 - sim_cen + N - 1;
min_bin = offset + PIX/2 - sim_cen;
if max_bin > PIX
    warning('Delta = %0.4f is to high and moves spectrum out of range.',delta);
    max_bin = PIX;
    min_bin = PIX - N + 1;
elseif min_bin < 1
    warning('Delta = %0.4f is to low and moves spectrum out of range.',delta);
    min_bin = 1;
    max_bin = N;
end
SIM_SPEC(min_bin:max_bin) = SX/sim_sum;