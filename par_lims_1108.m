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
Low_limit = [5e-3;
             5e-4;
             1.7e10;
             -0.3;
             0.033;
             85;
             0.4821;
             0.753;
             -25;
             -5;
             -0.04;
             0.02;
             4.0;
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
             0.7231;
             1.3984;
             -18;
             +5;
             -0.02;
             0.04;
             6.0;
             130;
             2000;
             0.01;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;

nPar = 16;

if length(par_name) ~= nPar
    error('Update number of parameters');
end