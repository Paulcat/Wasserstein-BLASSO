function [mu,nu] = ffw_ls(co)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

debug = 0;

[cxx,cyy,cxy,cx,cy] = deal(co{:});

denom = 4*cxx*cyy - cxy^2;
num1  = cxy*cy - 2*cyy*cx;
num2  = cxy*cx - 2*cxx*cy;

if denom == 0
    if cxx == 0
        mu = 0;
        nu = max(0,-cy/(2*cyy)); % unbounded lineseach
    else
        mu = max(0,-cx/(2*cxx));
        nu = 0;
    end
else
    mu = num1/denom;
    nu = num2/denom;
    
    if ~(mu>=0) || ~(nu>=0)
        f1 = cyy*nu^2 + cy*nu; % mu = 0?
        f2 = cxx*mu^2 + cx*mu; % or nu = 0?
        
        I = f1 < f2; % if I then mu=0, else nu=0
        mu = (1-I)*max(0,mu);
        nu = I*max(0,nu);
    end
end

if debug
    A = [cxx, cxy/2; cxy/2, cyy];
    B = [cx, cy];
    warning off
    S = sqrtm(A);
    warning on
    %cvx_solver mosek
    cvx_precision high
    cvx_begin %quiet
    variable X(2)
    %
    X >= 0
    %sum(X) <= 1
    %
    minimize( norm(S*X,'fro')^2 + B*X )
    cvx_end
    
    err1 = abs(X(1)-mu)/abs(X(1));
    err2 = abs(X(2)-nu)/abs(X(2));
    fprintf('\t linesearch error 1: %d\n', err1);
    fprintf('\t linesearch error 2: %d\n', err2);
end

end

