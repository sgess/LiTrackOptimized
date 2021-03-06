par_name = {'Initial bunch length';
            'Initial energy spread';
            'Initial number of particles';
            'Initial asymmetry';
            'NRTL Compressor Amplitude';
            'NRTL Compressor Phase';
            'NRTL Low-energy cutoff';
            'NRTL High-energy cutoff';
            'Decker Phase';
            'Ramp Phase';
            'S10 Low-energy cutoff';
            'S10 High-energy cutoff';
            'S20 Low-energy cutoff';
            'S20 High-energy cutoff';
            'S20 Beta';
            'YAG Dispersion';
            'YAG 2nd Order Dispersion';
            'Energy Offset';};

% Lower limit on params    
Low_limit = [5e-3;
             6e-4;
             1.7e10;
             -0.3;
             0.03;
             87;
             -0.04;
             0.02;
             -25;
             -5;
             -0.05;
             0.025;
             -0.04;
             0.02;
             3.5;
             100;
             -300;
             -0.01;];
         
% Higher limit on params                 
High_limit= [10e-3;
             12e-4;
             2.3e10;
             0.1;
             0.045;
             93;
             -0.02;
             0.04;
             -18;
             +5;
             -0.025;
             0.05;
             -0.02;
             0.04;
             6.5;
             135;
             2000;
             0.01;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;

nPar = 18;

if length(par_name) ~= nPar
    error('Update number of parameters');
end