function [eVecm,nit] = ffw_lmo(g,v0,options)
%FFW_LMO Linear minimization oracle over SDP cone
%   [EVEC,EVAL] = FFW_LMO(G,V0,options) produces the minimal eigenvalue
%   EVAL and a corresponding eigenvector EVEC of the gradient operator G,
%   using the Power Iterations algorithm.

% set options
maxit = getoptions(options,'maxiter',500);
tol   = getoptions(options,'tol',1e-9);

% initialize
v   = v0;
gv  = g(v);
ev  = v'*gv;
nit = 0;

while norm(gv-ev*v,'fro')/abs(ev) > tol && nit < maxit
    v  = gv / norm(gv,'fro');
    gv = g(v);
    ev = v'*gv;
    
    nit = nit+1;
end

if ev < 0 % then ev is the lowest eigenvalue
    eValm = ev;
    eVecm = v;
else
    e1 = ev;
    
    % re-initialize
    v  = v0;
    gv = g(v);
    ev = v'*gv;
    
    while norm(gv-ev*v,'fro')/abs(ev) > tol && nit < 2*maxit
        x  = gv - e1*v;
        v  = x/norm(x,'fro');
        gv = g(v);
        ev = v'*gv;
        
        nit = nit+1;
    end
    
    eValm = ev;
    eVecm = v;
end

%infos.niter = nit;


end

