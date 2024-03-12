function postprocessingplots(obs,distr,p_sim,p_phys,p_sch)

fprintf('Plotting the results: ');

if strcmp(p_sch.test, 'BKW')
    M1_numerics = cell2array(obs,1) ; 
    M2_numerics = cell2array(obs,2) ; 
    M4_numerics = cell2array(obs,5) ;    
	[~, ~, M4] = BKW_analytic(p_sim, p_sim.time_obs) ;
    
    % 1st & 2nd order moment
    figure(1)
    plot(p_sim.time_obs,sum(M1_numerics,1)./3,'ro','LineWidth',1.2,'MarkerSize',8)
    hold on
    plot(p_sim.time_obs, zeros(size(p_sim.time_obs)) ,'k-','LineWidth',1.5)
    plot(p_sim.time_obs,M2_numerics,'r*','LineWidth',1.2,'MarkerSize',8)
    plot(p_sim.time_obs, ones(size(p_sim.time_obs)) .* 3./2 ,'k--','LineWidth',1.5)
    hold off
    xlabel('$t$','Interpreter','latex','FontSize',15)
    title('First and Second order moment','Interpreter','latex','FontSize',15)
    legend('$\mathrm{M}1(t)$ DSMC', '$\mathrm{M}1(t)$ Exact', '$\mathrm{M}2(t)$ DSMC', '$\mathrm{M}2(t)$ Exact','Interpreter','latex','FontSize',15,'location','best')
    legend boxoff
    ylim([-0.5 2])

    % 4th order moment
    figure(2)
    plot(p_sim.time_obs,M4_numerics,'ro','LineWidth',1.2,'MarkerSize',8)
    hold on
    plot(p_sim.time_obs,M4,'k-','LineWidth',1.5)
    hold off
    xlabel('$t$','Interpreter','latex','FontSize',15)
    legend('$\mathrm{M}4(t)$ DSMC', '$\mathrm{M}4(t)$ Exact','Interpreter','latex','FontSize',15,'location','best')
    legend boxoff
    title('Fourth order moment','Interpreter','latex','FontSize',15)
    ylim([7.5 9.1])

    % Marginals at fixed times
    for i=0:p_sim.t_plt:p_sim.tf
        figure(i+3)

        f1D = distr{2,i+1};
        [~, f_BKW_1D, ~] = BKW_analytic(p_sim, i) ;
        
        plot(p_sim.VCells_exact, f_BKW_1D, 'k-','LineWidth',1.5)
        hold on
        plot(p_sim.VCells, f1D, 'ro','LineWidth',1.2,'MarkerSize',8)
        hold off
        legend('Exact BKW', 'DSMC', 'interpreter', 'latex', 'Location','northeast','FontSize',15)
        legend boxoff
        xlabel('$v_x$','interpreter', 'latex','FontSize',15)
        ylim([0 0.45])

        titl = sprintf('Marginal $f(v,t)$, $t=%d$', i);
        title(titl,'interpreter', 'latex','FontSize',15)

        drawnow
    end



elseif strcmp(p_sch.test, 'Trub')

    DT = Trubnikov_analytic(p_sim, p_phys, p_sch, p_sim.time_obs);
    Txyz_numerics = cell2array(obs,3) ;

    plot(p_sim.time_obs, DT, 'k-','LineWidth',1.5)
    hold on
    plot(p_sim.time_obs, (Txyz_numerics(1,:)-Txyz_numerics(3,:))./(Txyz_numerics(1,1)-Txyz_numerics(3,1)), 'ro','LineWidth',1.2,'MarkerSize',8)
    hold off
    xlabel('$t$','interpreter', 'latex','FontSize',15)
    ylabel('$\Delta T(t) / \Delta T(0)$','interpreter', 'latex','FontSize',15)
    legend('Trubnikov', 'DSMC', 'interpreter', 'latex', 'Location','northeast','FontSize',15)
    legend boxoff
 
    title('Trubnikov test','Interpreter','latex','FontSize',15)

end

fprintf('Done \n \n \n');

end