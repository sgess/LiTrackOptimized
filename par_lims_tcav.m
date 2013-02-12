par_name = {'Initial bunch length';
            'Initial energy spread';
            'Initial number of particles';
            'Initial asymmetry';
            'NRTL Compressor Amplitude';
            'NRTL Compressor Phase';
            'NRTL R56';
            'NRTL T566';
            'Decker Phase';
            'Ramp Phase';
            'S20 Low-energy cutoff';
            'S20 High-energy cutoff';
            'S20 Beta';
            'YAG Dispersion';
            'YAG 2nd Order Dispersion';
            'Energy Offset';};

% Lower limit on params    
Low_limit = [4e-3;
             4e-4;
             1.7e10;
             -0.4;
             0.040;
             88;
             0.58;
             0.9;
             -20;
             -5;
             -0.04;
             0.02;
             4.5;
             80;
             -1000;
             -0.01;];
         
% Higher limit on params                 
High_limit= [10e-3;
             12e-4;
             1.9e10;
             0.05;
             0.048;
             91;
             0.63;
             1.1;
             -17;
             +5;
             -0.02;
             0.04;
             6.0;
             160;
             1000;
             0.01;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;

nPar = 16;

if length(par_name) ~= nPar
    error('Update number of parameters');
end