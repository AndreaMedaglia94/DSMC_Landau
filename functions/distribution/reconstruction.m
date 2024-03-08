function f3D = reconstruction(vx, vy, vz, p_sim, p_phys)

% given vx, vy, vz, reconstruct the 3D distribution

counts = histcn([vx,vy,vz],p_sim.VEdges,p_sim.VEdges,p_sim.VEdges) ./ p_sim.N;

f3D    = counts ./ sum( counts .* p_sim.dV.^3, 'all' ) .* p_phys.rho ;


end