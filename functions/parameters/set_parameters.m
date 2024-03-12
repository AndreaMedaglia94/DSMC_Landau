function [p_phys, p_sim, p_sch] = set_parameters

fprintf('----------------------------------------------------------------------- \n');
fprintf('Direct Simulation Monte Carlo for the space homogeneous Landau equation \n\n');
fprintf('Code written by Andrea Medaglia (last update 12/03/2024) \n');
fprintf('https://doi.org/10.1016/j.jcp.2024.112845 \n');
fprintf('----------------------------------------------------------------------- \n \n \n');
fprintf('----------------------------------------------------------------------- \n');

p_sch = set_scheme_parameters;

p_phys = set_physical_parameters;

p_sim = set_simulation_parameters(p_sch, p_phys);
fprintf('----------------------------------------------------------------------- \n \n \n');
fprintf('----------------------------------------------------------------------- \n');
end