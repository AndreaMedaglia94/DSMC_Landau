function [mean, energy, temperature, temperaturetot, moment4] = Observables(vx, vy, vz, p_sim)

mean           = [sum(vx), sum(vy), sum(vz) ] ./ p_sim.N;

energy         = 0.5 .* sum( vx.^2 + vy.^2 + vz.^2) ./ p_sim.N;
 
temperature    = [sum(vx.^2), sum(vy.^2), sum(vz.^2) ] ./ p_sim.N;

temperaturetot = sum(temperature) ./3 ;

moment4        = sum( vx.^4 + vy.^4 + vz.^4) ./ p_sim.N;

end