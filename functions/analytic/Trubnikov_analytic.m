function DT = Trubnikov_analytic(p_sim, p_phys, p_sch, time)

if strcmp(p_sch.pot, 'Maxwell')
    tauT = 2/3 * p_phys.rho ;
elseif strcmp(p_sch.pot, 'Coulomb')
    tt   = pi * sqrt(2*p_phys.m) * 8 * p_phys.epsi0^2 * p_sim.Ttot^(3/2) / ( p_phys.Lambda * p_phys.rho * p_phys.e^4) ; 
    tauT = 5/8*sqrt(2*pi) * tt ;
end

DT = exp( - time ./ tauT ) ;

end

