function [obs, distr] = Solve(vx, vy, vz, p_sch, p_sim, p_phys, obs, distr) 

fprintf('Solving the equation: \n');

% time counters
counter_obs = 1 ;
counter_plt = 1 ;

% time cycle
for n=1:p_sim.ntot
    
    fprintf('t=%f\n',n*p_sim.dt);

    if p_sch.timeorder == 2
        [vx, vy, vz] = LandauMidPoint(vx, vy, vz, p_sch, p_sim, p_phys);
    elseif p_sch.timeorder == 1
        [vx, vy, vz] = NanbuBabovski(vx, vy, vz, p_sch, p_sim, p_phys);
    end
    
    if mod(n*p_sim.dt,p_sim.t_obs)==0
        counter_obs = counter_obs + 1 ;
        [obs{1,counter_obs}, obs{2,counter_obs}, obs{3,counter_obs}, obs{4,counter_obs}, obs{5,counter_obs}]= Observables(vx, vy, vz, p_sim);
    end
    
    if mod(n*p_sim.dt,p_sim.t_plt)==0
        if strcmp(p_sch.test, 'BKW')
            counter_plt = counter_plt + 1 ;
            distr{1,counter_plt}   = reconstruction(vx, vy, vz, p_sim, p_phys) ;
            distr{2,counter_plt}   = compute_marginal(vx, vy, vz, p_sim, p_phys);
        end
    end

 end    
    
fprintf('Done \n');
fprintf('----------------------------------------------------------------------- \n\n\n');
end