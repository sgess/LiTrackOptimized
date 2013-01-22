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
            'YAG 2nd Order Dispersion';};

% Lower limit on params    
Low_limit = [5e-3;
             4e-3;
             1.8e10;
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
             2;
             90;
             0;];
         
% Higher limit on params                 
High_limit= [10e-3;
             11e-3;
             2.4e10;
             0.0;
             0.045;
             93;
             -0.02;
             0.04;
             -15;
             +5;
             -0.025;
             0.05;
             -0.02;
             0.04;
             6;
             130;
             2000;];
         
Cent=(High_limit+Low_limit)/2;
Diff=High_limit-Low_limit;