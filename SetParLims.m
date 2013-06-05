function [Cent, Diff, lo_lims, hi_lims] = SetParLims(pars, frac, lo_lims, hi_lims)

if nargin < 3
    
    if isequal(size(frac),[1 1]);
        
        lo_abs = (1-frac)*abs(pars);
        hi_abs = (1+frac)*abs(pars);
        
    elseif isequal(size(frac),size(pars))
        
        lo_abs = (1-frac).*abs(pars);
        hi_abs = (1+frac).*abs(pars);
    else
        error('Limit sensitive must be a scalar or vector of the same size as "pars"');
    end
    
    lo_lims = (pars > 0).*lo_abs - (pars < 0).*hi_abs;
    hi_lims = (pars > 0).*hi_abs - (pars < 0).*lo_abs;
    
    tf = frac < 0;
    lo_lims(tf) = frac(tf);
    hi_lims(tf) = -frac(tf);
    
else
    
    tf = sum(lo_lims > hi_lims);
    if ~tf; error('Some low limits greater high limits'); end;
    
end

Cent=(hi_lims+lo_lims)/2;
Diff=hi_lims-lo_lims;