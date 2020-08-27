function [G,GU] = fgrad(m,cost,u1,u2,la,rho,tau)
%FGRAD Wasserstein-BLASSO gradient for FFW (1D)
%   Detailed explanation goes here

debug = 0;

% scaling
f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);

% helpers
N = 1./Dnumel2(m);
I1 = zeros(2*m-1); I1(:,1) = 1;
I2 = zeros(2*m-1); I2(1,:) = 1;
I  = I1 + I2;
%
Y = u1.*I1 + u2.'.*I2;
%

% linear part
C    = tau * cost;
C(1) = C(1)+1;

% gradient part: terms depending on TVALS
GT = @(T,h) 1/2/la * Tprod2(m,N.*(T.*I-Y),h) - 1/rho*Tprod2(m,T,h);

% gradient
G  = @(U,h)   f0*( Tprod2(m,N.*C,h) + GT(Tproj2(m,U),h) + 1/rho*U*(U'*h) );
GU = @(T,U,h) f0*( Tprod2(m,N.*C,h) + GT(T,h) + 1/rho*U*(U'*h) );

if debug
    disp('debug mode!!');
end

end

