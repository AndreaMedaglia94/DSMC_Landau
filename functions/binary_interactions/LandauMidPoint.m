function [wwx,wwy,wwz] = LandauMidPoint(wx, wy, wz, p_sch, p_sim, p_phys)

%%%%%%%%%%%%%%%%%%%%%%
%%% MIDPOINT STAGE %%%
%%%%%%%%%%%%%%%%%%%%%%

% Random coupling of the particles
j   = randperm(p_sim.N); 
wx  = wx(j);
wy  = wy(j);
wz  = wz(j);

% Pre-interaction opinions
wxi = wx(1:p_sim.N/2,:); wxj=wx(p_sim.N/2+1:p_sim.N,:);
wyi = wy(1:p_sim.N/2,:); wyj=wy(p_sim.N/2+1:p_sim.N,:);
wzi = wz(1:p_sim.N/2,:); wzj=wz(p_sim.N/2+1:p_sim.N,:);

% Binary interaction
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
%     A     = getA(par.tau,one_over_t,par.rho,N,M);
%     U     = rand(N/2,1); 
%     ct    = getCost(U,A,N,M);
elseif strcmp(p_sch.kernel, 'D2')
    tau0    = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    nu_tau0 = (1 - 2 .* tt) ;
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

wx_1 = [wwxi; wwxj ];
wy_1 = [wwyi; wwyj ];
wz_1 = [wwzi; wwzj ];

%%%%%%%%%%%%%%%%%%%%%%
%%%   n+1  STAGE   %%%
%%%%%%%%%%%%%%%%%%%%%%

% probabilities of interaction
p0 = 1 - p_sim.dt_tilde + p_sim.dt_tilde^2/2 ;                   
p1 = p_sim.dt_tilde - 3*p_sim.dt_tilde^2/2 + p_sim.dt_tilde^3/4 ;
p2 = p_sim.dt_tilde^2 - p_sim.dt_tilde^3/2 ;
p3 = p_sim.dt_tilde^3/4 ;

Ncol_1 = floor( p1 * p_sim.N / 2 ) *2 ;         % P(f_0,f_0)
Ncol_2 = floor( p2 * p_sim.N / 2 ) *2 ;         % P(f_0,f_1)
Ncol_3 = floor( p3 * p_sim.N / 2 ) *2 ;         % P(f_1,f_1)
Ncol_0 = p_sim.N - Ncol_1 - Ncol_2 - Ncol_3 ;   % no coll

% random coupling
% j   = randperm(p_sim.N); 
% wx  = wx(j); wx_1 = wx_1(j) ;
% wy  = wy(j); wy_1 = wy_1(j) ;
% wz  = wz(j); wz_1 = wz_1(j) ;

j0   = 1:Ncol_0 ;

j1   = Ncol_0+1:Ncol_1+Ncol_0 ;
j1_1 = j1(1:Ncol_1/2); j1_2 = j1(Ncol_1/2+1:Ncol_1);

j2   = Ncol_0+Ncol_1+1:Ncol_0+Ncol_1+Ncol_2 ;
j2_1 = j2(1:Ncol_2/2); j2_2 = j2(Ncol_2/2+1:Ncol_2);

j3   = Ncol_0+Ncol_1+Ncol_2+1:Ncol_0+Ncol_1+Ncol_2+Ncol_3 ;
j3_1 = j3(1:Ncol_3/2); j3_2 = j3(Ncol_3/2+1:Ncol_3);

% Pre-interaction opinions
% wx, wy, wz = f_0 
% wx_1, wy_1, wz_1 = f_1
wxi = [wx(j1_1,:); wx(j2_1,:);   wx_1(j3_1,:)] ;
wxj = [wx(j1_2,:); wx_1(j2_2,:); wx_1(j3_2,:)] ;

wyi = [wy(j1_1,:); wy(j2_1,:);   wy_1(j3_1,:)] ;
wyj = [wy(j1_2,:); wy_1(j2_2,:); wy_1(j3_2,:)] ;

wzi = [wz(j1_1,:); wz(j2_1,:);   wz_1(j3_1,:)] ;
wzj = [wz(j1_2,:); wz_1(j2_2,:); wz_1(j3_2,:)] ;

% Binary interaction
qx  = wxi-wxj; 
qy  = wyi-wyj; 
qz  = wzi-wzj; 

% Compute cos theta
qnorm = sqrt(qx.^2 + qy.^2 + qz.^2) ;

if strcmp( p_sch.pot, 'Maxwell') 
    one_over_tau = ones((Ncol_1+Ncol_2+Ncol_3)/2,1) ;
elseif strcmp( p_sch.pot, 'Coulomb') 
    one_over_tau = 4.*pi.*(p_phys.e.^2./(4.*pi.*p_phys.m./2.*p_phys.epsi0)).^2.*p_phys.rho.*p_phys.Lambda./(qnorm.^3) ;
end


if strcmp(p_sch.kernel, 'D1')
%     A     = getA(par.tau,one_over_t,par.rho,N,M);
%     U     = rand(N/2,1); 
%     ct    = getCost(U,A,N,M);
elseif strcmp(p_sch.kernel, 'D2')
    tau0    = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    nu_tau0 = (1 - 2 .* tt) ;
    ct    = nu_tau0 .*( tau0<=1 ) + (-1) .* ( tau0>1 );
elseif strcmp(p_sch.kernel, 'D3')
    tau0  = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    ct    =  1 - 2.*tanh( tau0 ) ;
end

% Compute the collisions
st = sin(acos(ct));

qperp = sqrt(qy.^2 + qz.^2);
ep    = 2.*pi.*rand((Ncol_1+Ncol_2+Ncol_3)/2,1);

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

wwx = [wx(j0,:); wwxi; wwxj];
wwy = [wy(j0,:); wwyi; wwyj];
wwz = [wz(j0,:); wwzi; wwzj];

end