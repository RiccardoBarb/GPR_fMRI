

function start_parpool(Nworkers)
%Nworkers = 4;

p=(gcp('nocreate'));

if isempty(p)
    parpool(Nworkers);
    maxNumCompThreads(Nworkers);
else
    fprintf('Already connected to %d Workers \n close and restart them \n',p.NumWorkers);
    delete(p);                                                             %felix 210819 delete pools to try if they then release their momory
    parpool(Nworkers);                                                     %felix 210819 restart pools
    maxNumCompThreads(Nworkers);
end

%turn warnings off
pctRunOnAll warning off