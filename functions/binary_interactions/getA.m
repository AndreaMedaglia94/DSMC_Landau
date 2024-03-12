function A = getA(p_sim, one_over_tau, p_phys)

s   = p_sim.epsi ./ p_phys.rho .* one_over_tau ;
A   = zeros(p_sim.N/2, 1);

options = optimset('TolX',1e-24);

for i=1:p_sim.N/2
    if s(i) < 5e-3
        A(i) = 1/s(i) ;
    elseif s(i) > 5e3
        A(i)  = 1e-8 ;
    else
        fun   = @(x) p_sim.nonlineq(x,s(i));
        A(i)  = fzero(fun, 1/s(i), options ); 
    end
end

end