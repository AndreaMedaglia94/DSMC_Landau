function [obs, distr] = set_observables(vx, vy, vz, p_sim, p_phys)

[obs{1,1}, obs{2,1}, obs{3,1}, obs{4,1}, obs{5,1}]= Observables(vx, vy, vz, p_sim);

distr{1,1}   = reconstruction(vx, vy, vz, p_sim, p_phys) ;
distr{2,1}   = compute_marginal(vx, vy, vz, p_sim, p_phys);

end