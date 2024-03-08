function [ct] = getCost(U,A,N,M)

global collision_type

if strcmp(collision_type, 'NB')
    ct = zeros(N/2,M);
    for i=1:N/2
        for k=1:M
            if A(i,k) < 1e-8
                ct(i,k) = 2.*U(i) - 1 - (2.*U(i) - 1).^2./2.* A(i,k) ; 
            elseif A(i,k) > 400
                ct(i,k) = 1 ; 
            else
                ct(i,k) = log( U(i).*exp(A(i,k)) + (1-U(i)).*exp(-A(i,k)) ) ./ A(i,k)  ; 
            end
            
        end
    end
elseif strcmp(collision_type, 'Bird')
    ct = zeros(1,M);
    for k=1:M
        if A(k) < 1e-12
            ct(k) = 2.*U - 1 - (2.*U - 1).^2./2.* A(k) ; 
        elseif A(k) > 650
            ct(k) = 1 ; 
        else
            ct(k) = log( U.*exp(A(k)) + (1-U).*exp(-A(k)) ) ./ A(k)  ; 
        end
        
    end  
end

end