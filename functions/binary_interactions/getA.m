function [A] = getA(tau,one_over_t,rho,N,M)

global collision_type

if strcmp(collision_type, 'NB')
    s       = tau ./ rho .* one_over_t ;
    A       = zeros(N/2,M);
    
    options = optimset('TolX',1e-24);
    
    for i=1:N/2
        for k=1:M
            fun     = @(x) (coth(x)-1/x-exp(-s(i,k))) ;
            if s(i,k) < 5e-3
                A(i,k) = 1/s(i,k) ;
            elseif s(i,k) > 5e3
                A(i,k)  = 1e-8 ;
            else
                A(i,k)  = fzero(fun, 1/s(i,k), options ); 
            end
        end
    end

elseif strcmp(collision_type, 'Bird')
    s       = tau ./ rho .* one_over_t ;
    A       = zeros(1,M);
    
    options = optimset('TolX',1e-24);

    for k=1:M
        fun     = @(x) (coth(x)-1/x-exp(-s(k))) ;
        if s(k) < 1e-4
            A(k) = 1/s(k) ;
        elseif s(k) > 5e3
            A(k)  = 1e-11 ;
        else
            A(k)  = fzero(fun, 1/s(k), options ); 
        end
    end
   
end

end