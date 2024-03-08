function [p_phys, p_sim, p_sch] = set_parameters

fprintf('Direct Simulation Monte Carlo for the space homogeneous Landau equation \n \n');

p_sch = set_scheme_parameters;

p_phys = set_physical_parameters;

p_sim = set_simulation_parameters(p_sch, p_phys);

end