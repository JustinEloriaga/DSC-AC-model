function [available,reason]=dsc_parallel_available()
%DSC_PARALLEL_AVAILABLE Read-only preflight, called only for parallel chains.
available=false;
reason='Parallel chains require an installed, licensed Parallel Computing Toolbox.';
try
    installed=~isempty(ver('parallel'));
    if installed&&license('test','Distrib_Computing_Toolbox')&& ...
            exist('parpool','file')~=0&&exist('parfeval','file')~=0
        available=true; reason='';
    end
catch problem
    reason=sprintf('%s %s',reason,problem.message);
end
end
