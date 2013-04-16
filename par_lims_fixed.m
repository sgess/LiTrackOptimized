par_name = {'Initial bunch length';
            'Initial energy spread';
            'Initial number of particles';
            'Initial asymmetry';
            'NRTL Compressor Amplitude';
            'NRTL Compressor Phase';
            '2-10 Phase';
            '11-20 Phase';
            'Ramp Phase';};

% Lower limit on params    
Low_limit = [7.0e-3;
             7e-4;
             1.6e10;
             -0.3;
             0.036;
             88;
             -27;
             -5;
             -5;];
         
% Higher limit on params                 
High_limit= [9e-3;
             9e-4;
             2.4e10;
             -0.05;
             0.042;
             93;
             -20;
             +5;
             +5;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;

nPar = 9;

if length(par_name) ~= nPar
    error('Update number of parameters');
end