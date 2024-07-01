function [vx, vy, vz] = InitialDataBKW(p_sim)

% sample initial data from the BKW solution by inverting the cumulative

% set some parameter to define the cumulative as if Ttot = 1 
K  = 1 * ( 1 - p_sim.BKW_C .* exp(-4 .* p_sim.BKW_B .* p_sim.BKW_t) );
norm = 1./(2.*pi.*K).^(3/2);
c2   = ( 1 - K ) ./ ( 2.*K.^2 ) ;
b  = 2 * K ;

% sample random numbers 
xi = rand(p_sim.N,1);
% and evaluate the inverse of the cumulative on this numbers
xx     = linspace(0,6,1000);
yy     = 4.*pi.*norm.*c2.*( 3./8.*sqrt(pi).*b.^(5./2).*erf(xx./sqrt(b)) - b./4.*exp(-xx.^2./b).*(3.*b.*xx + 2.*xx.^3) );
rr     = interp1(yy,xx,xi);

% sample random angles and rotate particles
theta = acos( 1 - 2 .* rand(p_sim.N,1) ) ;
phi   = rand(p_sim.N,1) .* 2.*pi;

vxr    = rr.*sin(theta).*cos(phi); 
vyr    = rr.*sin(theta).*sin(phi); 
vzr    = rr.*cos(theta); 

% rescale the particles with the total temperature
% it works only if Tx=Ty=Tz
vx = sqrt(p_sim.Ttot) .* vxr;
vy = sqrt(p_sim.Ttot) .* vyr;
vz = sqrt(p_sim.Ttot) .* vzr;

end