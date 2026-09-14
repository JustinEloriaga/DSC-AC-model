function v = dsc_random_walk_draw(z,initial_scale,innovation_sd)
% Exact factor for Cov(v_t,v_u)=innovation_sd^2*(initial_scale+min(t,u)).
% Equivalent to dense Brownian Cholesky multiplication, in O(T) storage/work.
z=z(:);
z(1)=sqrt(initial_scale+1)*z(1);
v=innovation_sd*cumsum(z);
end
