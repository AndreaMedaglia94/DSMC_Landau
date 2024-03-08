function [wx_hat,wy_hat,wz_hat]=Bird(wx_hat,wy_hat,wz_hat,par,PHI,W,N,M,dt,dtc)

global opt

kk  = randperm(N); 
wx_hat = wx_hat(kk,:); wy_hat = wy_hat(kk,:); wz_hat = wz_hat(kk,:);

tc = 0 ;
while tc < dt
    j1 = randi(N); j2 = randi(N);
    while j1==j2
        j2=randi(N);
    end
    
    % Pre-interaction opinions on the projections
    wxi_hat = wx_hat(j1,:); wxj_hat=wx_hat(j2,:);
    wyi_hat = wy_hat(j1,:); wyj_hat=wy_hat(j2,:);
    wzi_hat = wz_hat(j1,:); wzj_hat=wz_hat(j2,:);

    qx_hat  = wxi_hat-wxj_hat; 
    qy_hat  = wyi_hat-wyj_hat; 
    qz_hat  = wzi_hat-wzj_hat; 

    % Reconstruction of the velocities 
    wxMi =   wxi_hat * PHI; wxMj =   wxj_hat * PHI; 
    wyMi =   wyi_hat * PHI; wyMj =   wyj_hat * PHI;
    wzMi =   wzi_hat * PHI; wzMj =   wzj_hat * PHI;

    % Compute cos theta
    qx = wxMi - wxMj; qy = wyMi - wyMj; qz = wzMi - wzMj;
    qtilde = sqrt(qx.^2 + qy.^2 + qz.^2) ;
    
    one_over_t = 4.*pi.*(par.e.^2./(4.*pi.*par.m./2.*par.epsi)).^2.*par.rho.*0.5./(qtilde.^3) ;
    
    if strcmp(opt, 'D1')
        A = getA(par.tau,one_over_t,par.rho,N,M);
        U = rand; 
        ct = getCost(U,A,N,M);
    elseif strcmp(opt, 'D2')
        tt    = par.tau./(2*par.rho) .* one_over_t;
        nu_tt = (1 - 2 .* tt) ;
        ct = nu_tt.*(tt<=1) + (-1).*(tt>1);
    elseif strcmp(opt, 'D3')
        tt    = par.tau./(2*par.rho) .* one_over_t;
        ct = -2.*tanh(tt)+1;
    end
    
    % Compute the collision
    st = sin(acos(ct));

    qperp = sqrt(qy.^2 + qz.^2);
    ep    = 2.*pi.*rand;

    hx = qperp .* cos(ep);
    hy = -(qy.*qx.*cos(ep)+qtilde.*qz.*sin(ep))./qperp;
    hz = -(qz.*qx.*cos(ep)-qtilde.*qy.*sin(ep))./qperp;

    Wx_hat = (hx.*st.*W')*PHI'./sum(W);
    Wy_hat = (hy.*st.*W')*PHI'./sum(W);
    Wz_hat = (hz.*st.*W')*PHI'./sum(W);

    Vx_hat = (qx.*ct.*W')*PHI'./sum(W);
    Vy_hat = (qy.*ct.*W')*PHI'./sum(W);
    Vz_hat = (qz.*ct.*W')*PHI'./sum(W);

    wx_hat(j1,:)=wxi_hat -0.5.*(qx_hat -Vx_hat +Wx_hat);
    wy_hat(j1,:)=wyi_hat -0.5.*(qy_hat -Vy_hat +Wy_hat);
    wz_hat(j1,:)=wzi_hat -0.5.*(qz_hat -Vz_hat +Wz_hat);

    wx_hat(j2,:)=wxj_hat +0.5.*(qx_hat -Vx_hat +Wx_hat);
    wy_hat(j2,:)=wyj_hat +0.5.*(qy_hat -Vy_hat +Wy_hat);
    wz_hat(j2,:)=wzj_hat +0.5.*(qz_hat -Vz_hat +Wz_hat);
    
    tc = tc + dtc ;
end

end