function [vx, vy, vz] = InitialDataTrubnikov(p_sim)

% sample a standard Gaussian in 3D
v_sample = randn(p_sim.N,3);

% apply a conservative method to have momentum u and energy E
u = 0 ; 
E = 1 ;

v_prime = sum(v_sample) / p_sim.N ;  
E_prime = sum(v_sample.^2) / p_sim.N;
tau = sqrt( (E_prime-v_prime.^2/2)/(E-u.^2/2) );
lambda = v_prime - tau * u ;
v_sample = (v_sample - lambda)./tau ;

% apply again to go to machine accuracy
v_prime = sum(v_sample) / p_sim.N ;  
E_prime = sum(v_sample.^2) / p_sim.N;
tau = sqrt( (E_prime-v_prime.^2/2)/(E-u.^2/2) );
lambda = v_prime - tau * u ;
v_sample = (v_sample - lambda)./tau ;

% rescale the sample with the right temperatures
vx   = sqrt(p_sim.Tx) .* v_sample(:,1);
vy   = sqrt(p_sim.Ty) .* v_sample(:,2);
vz   = sqrt(p_sim.Tz) .* v_sample(:,3);


end