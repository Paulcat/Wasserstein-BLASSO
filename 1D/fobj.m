function F = fobj(m,cost,u1,u2,la,rho,tau)
%FOBJ Wasserstein-BLASSO objective for FFW (1D)
%   Detailed explanation goes here

debug = 0;

% scaling
f0 = 4*la / (norm(u1,'fro')^2 + norm(u2,'fro')^2);

% helper
d = length(m);
normT2 = @(T) sum( Dnumel2(m) .* abs(T).^2, 1:d);

% objective part: terms depending on TVALS
FT = @(T) real( tau * trace(cost'*T) + T(1) ) + ...
    1/4/la * ( norm(T(:,1)-u1,'fro')^2 + norm(T(1,:).'-u2,'fro')^2 ) + ...
    -1/2/rho * normT2(T);

% objective
F = @(U) f0 * ( FT(Tproj2(m,U)) + 1/2/rho * norm(U'*U,'fro')^2 );

if debug
    disp('debug mode!!');
    FT = tau*trace(cost'*T) + T(1);
    F = @(U) FT(Tproj2(m,U));
end

end

