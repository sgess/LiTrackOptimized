function f1 = par_plot(params,lo_lims,hi_lims,name,nPar,nmin,nmax,f1)

if nargin < 8
    f1 = figure;
else
    figure(f1);
end

if mod(nPar,2) == 0; k = 2; end
if mod(nPar,3) == 0; k = 3; end
if mod(nPar,4) == 0; k = 4; end
if mod(nPar,5) == 0; k = 5; end
if ~exist('k','var'); k = ceil(nPar/2); end
 

for i = 1:nPar 
    
    subplot(k,nPar/k,i);
    plot(params(i,nmin:nmax));
    line([nmin nmax],[lo_lims(i) lo_lims(i)],'color','r');
    line([nmin nmax],[hi_lims(i) hi_lims(i)],'color','r');
    title(name(i));
    
end