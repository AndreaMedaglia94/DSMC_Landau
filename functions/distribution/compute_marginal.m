function f1D = compute_marginal(vx, vy, vz, p_sim, p_phys)

% given the 3D distribution, compute the marginal in x

f3D = reconstruction(vx, vy, vz, p_sim, p_phys) ;

f1D = sum( sum( f3D , 3) , 2) .* p_sim.dV.^2 ;

end