function [dE,zc,sigz,W] = quick_wake(z,L,Ne,Nbin,fn)
%===============================================================================
%        [dE,zc,sigz] = quick_wake(z[,L,Ne,Nbin,fn]);
%
%	Function to return the wakefield induced energy profile vs. z for
%	a set of given axial coordinates "z".
%
%  INPUTS:	z:		The internal axial coordinates, within the bunch, of
%					each electron with respect to any fixed point [m]
%			L:		(Optional, DEF=1 m) The length of the linac [m]
%			Ne:		(Optional, DEF=1  ) The number of electrons in the bunch
%			Nbin:   (Optional, DEF=100) The number of bins to use
%			fn:		(Optional, DEF='slac.dat') File name containing longitudinal
%					point wake (DEF='slac.dat')
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			sigz:	rms bunch length (standard deviation) [m]
%===============================================================================

% Number of simulated particles
nn = length(z);

% RMS longitudinal spread
sigz = std(z);

% Wakefile: A(:,1) = Z (m), A(:,2) = W(Z) (V/C/m)
A  = load(fn);
nA = length(A);

% Histogram particles into Nbins with zeros at ends for interpolation
[N,zc] = hist(z,Nbin-2);
% Bin spacing (m)
dzc = zc(2)-zc(1);
% Zero padded axis
zc = [zc(1)-dzc zc zc(Nbin-2)+dzc];
% Zero padded histogram
N = [0 N 0];

% max Z in wake file (last point usually projected estimate)
maxz_fn = A(nA-1);
if (zc(Nbin)-zc(1)) > maxz_fn
    disp(' ')
    disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and ' fn ' is inaccurate there'])
    disp(' ')
end

% delta vector
dE = zeros(Nbin,1);
% SI electric charge
e = 1.6022E-19;
% scale variable (mC/V ?)
scl = -L*(Ne/nn)*e*1E-6;

% wake vector
W = zeros(Nbin,Nbin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over beam profile bins %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Particles in bin j are acted on by wake
for j = 1:Nbin
    % Particles in bin k produce wake felt by particles in bin j
    for k = 1:Nbin
        
%         hold off;
%         pause;
%         semilogy(A(:,1),A(:,2));
%         axis([0.005 0.01 1e13 1e14]);
%         hold on;
        
        % For k > j, there is no wake because particles in bin k are behind bin j
        if k > j
            break
        end
        % Separation between bins (m)
        dz = dzc*(j - k);
        % Bin in wakefile that matches separation
        [ddz,ii] = min(abs(A(:,1)-dz));
        % Sign of bin for interpolation
        ddz = sign(A(ii,1)-dz)*ddz;
        if ddz > 0
            % Indices of interpolation
            i1 = ii - 1;
            i2 = ii;
            if i1 == 0
                error('Index into zeroth entry of wakefield array - try finer steps in wake file');
            end
            % Z points of interpolation
            dz1 = A(i1,1);
            dz2 = A(i2,1);
            % Wake points of interpolation
            W1  = A(i1,2);
            W2  = A(i2,2);
            
%             plot(dz1,W1,'r*');
%             plot(dz2,W2,'r*');
            
            % Wake value
            W(j,k)   = W2 - (W2-W1)*ddz/(dz2-dz1);
            
%             plot(dz2+ddz,W,'g*');
            
        elseif ddz < 0
            % Indices of interpolation
            i1 = ii;
            i2 = ii + 1;
            if i2 > length(A(:,1))
                error('WARN: Index to a point beyond wakefield array - try extending the wake file');
            end
            % Z points of interpolation
            dz1 = A(i1,1);
            dz2 = A(i2,1);
            % Wake points of interpolation
            W1  = A(i1,2);
            W2  = A(i2,2);
            
%             plot(dz1,W1,'m*');
%             plot(dz2,W2,'m*');
            
            % Wake value
            W(j,k)   = W1 - (W2-W1)*ddz/(dz2-dz1);
            
%             plot(dz1+ddz,W,'k*');

        else
            % use 1/2 wake for self loading of bin
            W(j,k)   = A(ii,2)/2;		
        end
        % The sum of the fields*number of particles 
        dE(j) = dE(j) + scl*N(k)*W(j,k);		
    end
end