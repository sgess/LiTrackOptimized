%p_z = [	% 12% peak current jitter: for CDR (~DEC 01)
%	0.67 0.10	% L0 RF phase [deg]
%	0.17 0.10	% L1 RF phase [deg]
%	1.20 0.30	% Lx RF phase [deg]
%	0.23 0.07	% L2 RF phase [deg]
%	12.0 1.00	% L3 RF phase [deg]
%	0.29 0.10	% L0 RF voltage [%]
%	0.35 0.10	% L1 RF voltage [%]
%	1.70 0.25	% Lx RF voltage [%]
%	1.10 0.10	% L2 RF voltage [%]
%	12.0 1.00	% L3 RF voltage [%]
%	0.14 0.02	% BC1 chicanes [%]
%	0.78 0.02	% BC2 chicanes [%]
%	1.80 0.70	% gun timing [psec]
%	5.60 2.00	% bunch charge [%]
%									];

%p_E = [	% 0.10 % relative e- beam energy jitter: for CDR (~DEC 01)
%	3.70 0.10	% L0 RF phase [deg]
%	0.24 0.10	% L1 RF phase [deg]
%	10.0 0.80	% Lx RF phase [deg]
%	0.36 0.07	% L2 RF phase [deg]
%	0.48 0.07	% L3 RF phase [deg]
%	0.32 0.10	% L0 RF voltage [%]
%	0.34 0.10	% L1 RF voltage [%]
%	2.00 0.25	% Lx RF voltage [%]
%	0.64 0.07	% L2 RF voltage [%]
%	0.15 0.05	% L3 RF voltage [%]
%	0.16 0.02	% BC1 chicanes [%]
%	0.65 0.02	% BC2 chicanes [%]
%	1.30 0.83	% gun timing [psec]
%	46.0 5.00	% bunch charge [%]
%									];

%p_z = [	% 1-nC nominal case: 12% peak current jitter: new for 135 MeV inj., 26JUN03
%	0.71 0.10	% L0 RF phase [deg]
%	0.19 0.10	% L1 RF phase [deg]
%	1.55 0.50	% Lx RF phase [deg]
%	0.23 0.07	% L2 RF phase [deg]
%	8.89 0.15	% L3 RF phase [deg]
%	0.27 0.10	% L0 RF voltage [%]
%	0.27 0.10	% L1 RF voltage [%]
%	1.29 0.25	% Lx RF voltage [%]
%	1.12 0.10	% L2 RF voltage [%]
%	6.63 0.08	% L3 RF voltage [%]
%	0.14 0.01	% BC1 chicanes [%]
%	0.78 0.05	% BC2 chicanes [%]
%	3.99 0.80	% gun timing [psec]
%	5.64 2.00	% bunch charge [%]
%									];

%p_E = [	% 1-nC nominal case: 0.10 % relative e- beam energy jitter: new for 135 MeV inj., 26JUN03
%	4.05 0.10	% L0 RF phase [deg]
%	0.36 0.10	% L1 RF phase [deg]
%	5.34 0.50	% Lx RF phase [deg]
%	0.69 0.07	% L2 RF phase [deg]
%	0.35 0.15	% L3 RF phase [deg]
%	0.35 0.10	% L0 RF voltage [%]
%	0.32 0.10	% L1 RF voltage [%]
%	2.26 0.25	% Lx RF voltage [%]
%	1.29 0.10	% L2 RF voltage [%]
%	0.15 0.08	% L3 RF voltage [%]
%	0.16 0.01	% BC1 chicanes [%]
%	0.65 0.05	% BC2 chicanes [%]
%	1.69 0.80	% gun timing [psec]
%	53.0 2.00	% bunch charge [%]
%									];
									
p_z = [	% 200-pC case: 12% peak current jitter, L3-phase=-6.6 deg
	2.52 0.30	% L0 RF phase [deg]
	0.94 0.10	% L1 RF phase [deg]
	2.50 0.75	% Lx RF phase [deg]
	0.57 0.10	% L2 RF phase [deg]
	2.57 0.15	% L3 RF phase [deg]
	0.71 0.10	% L0 RF voltage [%]
	0.70 0.10	% L1 RF voltage [%]
	2.65 0.25	% Lx RF voltage [%]
	2.30 0.10	% L2 RF voltage [%]
	3.40 0.07	% L3 RF voltage [%]
	0.14 0.01	% BC1 chicanes [%]
	0.78 0.02	% BC2 chicanes [%]
	4.76 0.80	% gun timing [psec]
   14.81 5.00	% bunch charge [%]
									];
  
  p_E = [	% 200-pC case: 0.10 % relative e- beam energy jitter, L3-phase=-6.6 deg
	3.48 0.30	% L0 RF phase [deg]
	0.35 0.10	% L1 RF phase [deg]
	2.89 0.75	% Lx RF phase [deg]
	0.29 0.10	% L2 RF phase [deg]
	0.72 0.15	% L3 RF phase [deg]
	0.29 0.10	% L0 RF voltage [%]
	0.27 0.10	% L1 RF voltage [%]
	1.84 0.25	% Lx RF voltage [%]
	0.46 0.10	% L2 RF voltage [%]
	0.15 0.07	% L3 RF voltage [%]
	0.16 0.01	% BC1 chicanes [%]
	0.65 0.02	% BC2 chicanes [%]
	2.18 0.80	% gun timing [psec]
	92.1 5.00	% bunch charge [%]
  									];

	r_z = p_z(:,2)./p_z(:,1);
	r_E = p_E(:,2)./p_E(:,1);
	
	n_z = length(r_z);
	n_E = length(r_E);
	
	disp(' ')
	disp(sprintf('sigZ budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_z'*r_z),n_z))
	disp(sprintf('dE/E budget (should be < 1): %5.3f; (%3.0f parameters)',sqrt(r_E'*r_E),n_E))
	disp(' ')
