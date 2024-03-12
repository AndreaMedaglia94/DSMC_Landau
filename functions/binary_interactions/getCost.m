function [ct] = getCost(U,A)

ct = (A<1e-8) .* (2.*U-1-(2.*U-1).^2./2.*A) + (A>400).*1 + (A>1e-8).*(A<400).*log(U.*exp(A)+(1-U).*exp(-A))./A;

end