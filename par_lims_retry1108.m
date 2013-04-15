par_name = {'Initial bunch length';
            'Initial energy spread';
            'Initial number of particles';
            'Initial asymmetry';
            'NRTL Compressor Amplitude';
            'NRTL Compressor Phase';
            'NRTL R56';
            'NRTL T566';
            '2-10 Phase';
            '11-20 Phase';
            'Ramp Phase';
            'S20 Beta';
            'YAG Dispersion';
            'YAG 2nd Order Dispersion'};

% Lower limit on params    
Low_limit = [7.5e-3;
             7e-4;
             1.7e10;
             -0.3;
             0.036;
             89;
             0.58;
             0.8;
             -27;
             -5;
             -5;
             4;
             100;
             -2000;];
         
% Higher limit on params                 
High_limit= [9e-3;
             9e-4;
             2.2e10;
             -0.05;
             0.040;
             91;
             0.63;
             1.2;
             -20;
             +5;
             +5;
             6.0;
             130;
             2000;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;

nPar = 14;

if length(par_name) ~= nPar
    error('Update number of parameters');
end