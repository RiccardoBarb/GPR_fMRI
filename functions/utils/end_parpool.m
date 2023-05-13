function end_parpool
p=(gcp('nocreate'));
if ~isempty(p)
    delete(p);
end
end