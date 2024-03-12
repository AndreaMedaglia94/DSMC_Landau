function [wwx,wwy,wwz]=NanbuBabovski(wx, wy, wz, p_sch, p_sim, p_phys)

% Random coupling of the particles
j  = randperm(p_sim.N); 
j1 = j(1:p_sim.N/2); 
j2 = j(p_sim.N/2+1:p_sim.N);

% Pre-interaction opinions on the projections
wxi = wx(j1,:); wxj=wx(j2,:);
wyi = wy(j1,:); wyj=wy(j2,:);
wzi = wz(j1,:); wzj=wz(j2,:);

qx  = wxi-wxj; 
qy  = wyi-wyj; 
qz  = wzi-wzj; 

% Compute cos theta
qnorm = sqrt(qx.^2 + qy.^2 + qz.^2) ;

if strcmp( p_sch.pot, 'Maxwell') 
    one_over_tau = ones(p_sim.N/2,1) ;
elseif strcmp( p_sch.pot, 'Coulomb') 
    one_over_tau = 4.*pi.*(p_phys.e.^2./(4.*pi.*p_phys.m./2.*p_phys.epsi0)).^2.*p_phys.rho.*p_phys.Lambda./(qnorm.^3) ;
end


if strcmp(p_sch.kernel, 'D1')
    if strcmp( p_sch.pot, 'Maxwell') 
        U   = rand(p_sim.N/2,1); 
        ct  = 1 ./ p_sim.A .* log( exp(-p_sim.A) + 2 .* U .* sinh(p_sim.A) );
    elseif strcmp( p_sch.pot, 'Coulomb') 
        A   = getA(p_sim, one_over_tau, p_phys);
        U   = rand(p_sim.N/2,1);  
        ct  = getCost(U,A);
    end
elseif strcmp(p_sch.kernel, 'D2')
    tau0    = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    nu_tau0 = (1 - 2 .* tau0) ;
    ct    = nu_tau0 .*( tau0<=1 ) + (-1) .* ( tau0>1 );
elseif strcmp(p_sch.kernel, 'D3')
    tau0  = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    ct    =  1 - 2.*tanh( tau0 ) ;
end

% Compute the collisions
st = sin(acos(ct));

qperp = sqrt(qy.^2 + qz.^2);
ep    = 2.*pi.*rand(p_sim.N/2,1);

hx = qperp .* cos(ep);
hy = -(qy.*qx.*cos(ep)+qnorm.*qz.*sin(ep))./qperp;
hz = -(qz.*qx.*cos(ep)-qnorm.*qy.*sin(ep))./qperp;

% compute the post-collisional velocities
wwxi=wxi-0.5.*(qx.*(1-ct)+hx.*st);
wwyi=wyi-0.5.*(qy.*(1-ct)+hy.*st);
wwzi=wzi-0.5.*(qz.*(1-ct)+hz.*st);

wwxj=wxj+0.5.*(qx.*(1-ct)+hx.*st);
wwyj=wyj+0.5.*(qy.*(1-ct)+hy.*st);
wwzj=wzj+0.5.*(qz.*(1-ct)+hz.*st);

wwx = [wwxi; wwxj];
wwy = [wwyi; wwyj];
wwz = [wwzi; wwzj];

end