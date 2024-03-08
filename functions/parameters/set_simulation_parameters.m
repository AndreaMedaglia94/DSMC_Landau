function p_sim = set_simulation_parameters(p_sch, p_phys)

% number of particles
p_sim.N  = 1e6; 


% % % time paramteres % % %

% final time of the simulation
p_sim.tf  = 5;

% time step
p_sim.dt  = 0.1; 

% parameter epsilon approximating the Boltzmann equation
p_sim.epsi = p_phys.rho * p_sim.dt; 

% total number of steps
p_sim.ntot   = ceil(p_sim.tf / p_sim.dt);

% if Bird -> set the small time update dt_c 
if strcmp(p_sch.coll, 'B')
    p_sim.dtc  = 2 * p_sim.epsi / ( p_phys.rho * p_sim.N ); 
end

% time step for the observables
p_sim.t_obs    = p_sim.dt; 
p_sim.time_obs = 0:p_sim.t_obs:p_sim.tf;

% time step for the plot
p_sim.t_plt  = 1; 


% % % parameters for initial data % % %

% temperature along the x-y-z axis
p_sim.Tx = 1 ;
p_sim.Ty = 1 ;
p_sim.Tz = 1 ; 

% total initial temperature
p_sim.Ttot = (p_sim.Tx + p_sim.Ty + p_sim.Tz) / 3 ;

% if BKW -> set the constant for the analytic solution
if strcmp(p_sch.test, 'BKW')
    p_sim.BKW_B = 1/8;
    p_sim.BKW_C = 2/5;
    p_sim.BKW_t = 0  ;
    if p_sim.Tx ~= p_sim.Ty || p_sim.Tx ~= p_sim.Tz || p_sim.Tz ~= p_sim.Ty
        fprintf('Error in the BKW test, temperature along axis are not equal \n');
        stop
    end
end


% % % parameters for the reconstruction % % %

% right boundary of the v-domain for the reconstruction
p_sim.L     = 5;

% number of bins per direction for the histogram reconstruction
p_sim.Nbins  = 50;
p_sim.Nexact = 100;

% edges and step for the histogram reconstruction 
p_sim.VEdges = linspace(-p_sim.L, p_sim.L, p_sim.Nbins);
p_sim.dV     = p_sim.VEdges(2) - p_sim.VEdges(1) ;
p_sim.VCells = p_sim.VEdges(1:end-1) + p_sim.dV/2 ;

% edges and step for the exact solution
p_sim.VEdges_exact = linspace(-p_sim.L, p_sim.L, p_sim.Nexact);
p_sim.dV_exact     = p_sim.VEdges_exact(2) - p_sim.VEdges_exact(1) ;
p_sim.VCells_exact = p_sim.VEdges_exact(1:end-1) + p_sim.dV_exact/2 ;
[p_sim.Vx_exact,p_sim.Vy_exact,p_sim.Vz_exact] = meshgrid(p_sim.VCells_exact,p_sim.VCells_exact,p_sim.VCells_exact);

end