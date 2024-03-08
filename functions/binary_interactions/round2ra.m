function Y = round2ra(X)
% Round towards the nearest integer (random rounding / stochastic rounding)
%
% (c) 2013 Stephen Cobeldick
%
% Rounds the elements of X to the nearest integers. X may be an N-D matrix.
% Elements with a fraction of 0.5 round randomly to the nearest integer.
% For complex X, the imaginary and real parts are rounded independently.
%
% Syntax:
%  Y = round2ra(X)
%
% See also ROUND2EV ROUND2OD ROUND2DN ROUND2UP ROUND2ZE ROUND2SF ROUND2DP ROUND60063 DATEROUND ROUND

if isreal(X)
    Y = 0.5==mod(X,1);
    Z = double(Y);
    Z(Y) = 0.5-binornd(1,0.5,sum(Z(:)),1);
else
    Q = 0.5==mod(real(X),1);
    I = 0.5==mod(imag(X),1);
    R = double(Q);
    J = double(I);
    R(Q) = 0.5-binornd(1,0.5,sum(R(:)),1);
    J(I) = 0.5-binornd(1,0.5,sum(J(:)),1);
    Z = complex(R,J);
end
Y = round(X-Z);
%----------------------------------------------------------------------End!