function r = concat_chains(r_chains)
% Concatenate MCMC draw structs returned by tvsvar_modified_msv2_gam2_gen*
% along the draw dimension. Each chain has identical field shapes.

if numel(r_chains) == 1
    r = r_chains{1};
    return
end

nc = numel(r_chains);
B_c     = cell(nc, 1);
V_c     = cell(nc, 1);
rt_c    = cell(nc, 1);
P_c     = cell(nc, 1);
sig2r_c = cell(nc, 1);
h_c     = cell(nc, 1);
sig2h_c = cell(nc, 1);
for c = 1:nc
    B_c{c}     = r_chains{c}.B;
    V_c{c}     = r_chains{c}.V;
    rt_c{c}    = r_chains{c}.r;
    P_c{c}     = r_chains{c}.P;
    sig2r_c{c} = r_chains{c}.sig2r;
    h_c{c}     = r_chains{c}.h;
    sig2h_c{c} = r_chains{c}.sig2h;
end

r.B     = cat(3, B_c{:});
r.V     = cat(3, V_c{:});
r.r     = cat(3, rt_c{:});
r.P     = cat(4, P_c{:});
r.sig2r = vertcat(sig2r_c{:});
r.h     = cat(3, h_c{:});
r.sig2h = vertcat(sig2h_c{:});
end
