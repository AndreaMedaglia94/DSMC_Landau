function [vx, vy, vz] = InitialData(p_sch, p_sim, p_phys)

fprintf('Initializing the particles: ');

if strcmp(p_sch.test, 'BKW')
    [vx, vy, vz] = InitialDataBKW(p_sim);
elseif strcmp(p_sch.test, 'Trub')
    [vx, vy, vz] = InitialDataTrubnikov(p_sim);
else
    fprintf('Error in the selection of the test to execute \n');
    stop
end

if strcmp(p_sch.init_conditions, 'YES')
    if strcmp(p_sch.test, 'BKW')
        f1D = compute_marginal(vx, vy, vz, p_sim, p_phys) ;
        [~, f_BKW_1D, ~] = BKW_analytic(p_sim, 0) ;
        
        plot(p_sim.VCells_exact, f_BKW_1D, 'k-')
        hold on
        plot(p_sim.VCells, f1D, 'ro')
        hold off
        legend('BKW', 'DSMC', 'interpreter', 'latex', 'Location','best')
        legend boxoff
        xlabel('$v_x$','interpreter', 'latex')
        title('Marginal $f(v,t)$, $t=0$','interpreter', 'latex')
        drawnow
        
    end
end

fprintf('Done \n');

end