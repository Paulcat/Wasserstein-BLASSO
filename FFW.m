function [U] = FFW(problem,options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



% load problem parameters
m  = problem.vardim;
M  = prod(m);
%
f   = problem.fobj;
g   = problem.grad;
gt  = problem.grad_pre;
lco = problem.ls;
%
options_lmo  = set_lmo_options (options);
options_bfgs = set_bfgs_options(options);

% set options
U0    = getoptions(options,'init',zeros(M,1));
Om    = getoptions(options,'Om',ones(M,1));
maxit = getoptions(options,'maxiter',20);
tol   = getoptions(options,'tol',1e-5);

% *** Initialization ***
niter = 0;
crit  = -Inf;
U     = U0;
v0    = ones(M,1)/sqrt(M);


% DISPLAY
fprintf('\n\n')
fprintf('--------------------- FFW Algorithm -------------------------\n')
fprintf('-------------------------------------------------------------\n');
fprintf('IT  OBJ \t GAP \t    PI   (TIME)  LS \t  BFGS (TIME)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('%-3i %-+4.4e\n',niter,f(U0))

while crit < 0 && niter < maxit
    
    % *** Linear Minimization Oracle ***
    %T = Tproj2(m,U);
    T = Tproj4(m,U);
    glmo = @(h) Om .* gt(T,U,Om.*h);
    %
    tic;
    [eVecm,nPI] = ffw_lmo(glmo,v0,options_lmo);
    time_lmo = toc;
    eVecm = Om .* eVecm;
    
    % stopping criterion
    ege  = eVecm'*g(U,eVecm);
    crit = ege;
    
    
    % *** Frank-Wolfe update (with linesearch) ***
    %t  = Tproj2(m,eVecm);
    t  = Tproj4(m,eVecm);
    co = lco(U,T,eVecm,t);
    [mu,nu] = ffw_ls(co);
    if mu==0
        U = sqrt(nu)*eVecm;
    elseif nu==0
        U = sqrt(mu)*U;
    else
        U = [sqrt(mu)*U, sqrt(nu)*eVecm];
    end
    
    
    % *** BFGS step ***
    if options_bfgs.on
        tic;
        [U,nBFGS] = ffw_bfgs(U,f,@(U)2*g(U,U),options_bfgs);
        time_bfgs = toc;
    end
    
    % update monitors
    niter = niter+1;
    
    % display
    fprintf('%-3i %-+4.4e  %-+4.2e  %-4i (%4.1f)  %-.1e  %-4i (%4.1f)\n', ...
        niter,f(U),crit,nPI,time_lmo,mu/nu,nBFGS,time_bfgs);
end


end

function opt_bfgs = set_bfgs_options(options)
opt_bfgs.on              = getoptions(options, 'bfgsOn', 1);
opt_bfgs.display         = 'off';                                    % off | final | iter | full | excessive
opt_bfgs.optTol          = 1e-16;                                    % first-order optimality (exitflag 1)
opt_bfgs.progTol         = getoptions(options, 'bfgsProgTol', 1e-8); % parameters change (exitflag 2)
opt_bfgs.MaxFunEvals     = 100000;
opt_bfgs.MaxIter         = getoptions(options, 'bfgsMaxIter', 500);  % (exitflag 0)
opt_bfgs.Method          = 'lbfgs';
opt_bfgs.DerivativeCheck = 'off';
opt_bfgs.Corr            = 15;
opt_bfgs.Damped          = 0;
opt_bfgs.numDiff         = 0; % use-provided gradient
end

function opt_lmo = set_lmo_options(options)
opt_lmo.tol     = getoptions(options, 'lmoTol', 1e-16);
opt_lmo.maxiter = getoptions(options, 'lmoMaxIter', 1000); 
end

