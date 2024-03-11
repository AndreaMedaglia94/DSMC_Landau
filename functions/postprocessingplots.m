function postprocessingplots(obs,p_sim,p_phys,p_sch)

if strcmp(p_sch.test, 'BKW')
    M4_numerics = cell2array(obs,5) ;    
	[~, ~, M4] = BKW_analytic(p_sim, p_sim.time_obs) ;
    
    plot(p_sim.time_obs,M4,'k-')
    hold on
    plot(p_sim.time_obs,M4_numerics,'ro')
    
elseif strcmp(p_sch.test, 'Trub')

    DT = Trubnikov_analytic(p_sim, p_phys, p_sch, p_sim.time_obs);
    Txyz_numerics = cell2array(obs,3) ;

    plot(p_sim.time_obs, DT, 'k-')
    hold on
    plot(p_sim.time_obs, (Txyz_numerics(1,:)-Txyz_numerics(3,:))./(Txyz_numerics(1,1)-Txyz_numerics(3,1)), 'ro')
    hold off

end

end