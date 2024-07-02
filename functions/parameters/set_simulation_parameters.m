function p_sim = set_simulation_parameters(p_sch, p_phys)

% number of particles
p_sim.N  = 1e6; 


% % % time paramteres % % %

% final time of the simulation
if strcmp(p_sch.test, 'BKW')
    p_sim.tf  = 12;
elseif strcmp(p_sch.test, 'Trub')
    if strcmp(p_sch.pot, 'Maxwell')
        p_sim.tf  = 3;
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.tf  = 14;
    end
end

% time step
p_sim.dt  = 0.05; 

% parameter epsilon approximating the Boltzmann equation
if p_sch.timeorder == 2
    p_sim.epsi = p_phys.rho * p_sim.dt ./ (3-sqrt(5)) ; 
elseif p_sch.timeorder == 1
    p_sim.epsi = p_phys.rho * p_sim.dt  ; 
end

p_sim.dt_tilde = p_phys.rho .* p_sim.dt ./ p_sim.epsi ;

% total number of steps
p_sim.ntot   = ceil(p_sim.tf / p_sim.dt);

% time step for the observables
p_sim.t_obs    = p_sim.dt; 
p_sim.time_obs = 0:p_sim.t_obs:p_sim.tf;

% time step for the plot
p_sim.t_plt  = 1; 


% % % parameters for initial data % % %

% temperature along the x-y-z axis
if strcmp(p_sch.test, 'BKW')
    p_sim.Tx = 2 ;
    p_sim.Ty = 2 ;
    p_sim.Tz = 2 ; 

    p_sim.BKW_B = 1/8;
    p_sim.BKW_C = 2/5;
    p_sim.BKW_t = 0  ;
elseif strcmp(p_sch.test, 'Trub')
    if strcmp(p_sch.pot, 'Maxwell')
        p_sim.Tx = 0.085 ;
        p_sim.Ty = 0.085 ;
        p_sim.Tz = 0.04 ; 
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.Tx = 0.07 ;
        p_sim.Ty = 0.07 ;
        p_sim.Tz = 0.05 ; 
    end
end

% total initial temperature
p_sim.Ttot = (p_sim.Tx + p_sim.Ty + p_sim.Tz) / 3 ;

% % % parameters for the reconstruction % % %

% right boundary of the v-domain for the reconstruction
p_sim.L     = 5;

% number of bins per direction for the histogram reconstruction
p_sim.Nbins  = 100;
p_sim.Nexact = p_sim.Nbins;

% edges and step for the histogram reconstruction 
p_sim.VEdges = linspace(-p_sim.L, p_sim.L, p_sim.Nbins+1);
p_sim.dV     = p_sim.VEdges(2) - p_sim.VEdges(1) ;
p_sim.VCells = p_sim.VEdges(1:end-1) + p_sim.dV/2 ;

% edges and step for the exact solution
p_sim.VEdges_exact = linspace(-p_sim.L, p_sim.L, p_sim.Nexact+1);
p_sim.dV_exact     = p_sim.VEdges_exact(2) - p_sim.VEdges_exact(1) ;
p_sim.VCells_exact = p_sim.VEdges_exact(1:end-1) + p_sim.dV_exact/2 ;
[p_sim.Vx_exact,p_sim.Vy_exact,p_sim.Vz_exact] = meshgrid(p_sim.VCells_exact,p_sim.VCells_exact,p_sim.VCells_exact);

% % % other parameters % % %

% NonLinear equation for A
if strcmp(p_sch.kernel, 'D1')
    if strcmp(p_sch.pot, 'Maxwell')
        nonlineq   = @(x) (coth(x)-1./x-exp(-p_sim.epsi)) ;
        options    = optimset('TolX',1e-24);
        p_sim.A    = fzero(nonlineq, 1./p_sim.epsi,options);
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.nonlineq =  @(x,y) (coth(x)-1/x-exp(-y)) ;
    end
end

end