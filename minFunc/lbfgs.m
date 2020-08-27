function [x,f,exitflag,output] = lbfgs(funObj,x0,options,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    options = [];
end

% get parameters
[verbose,verboseI,debug,optTol,progTol,maxFunEvals,maxIter,Damped,...
    corrections,useMex,HvFunc,LS_init,LS_interp,LS_multi,c1,c2] = processInputOptions(options);

LBFGS = 5;

% initialize
p = length(x0);
d = zeros(p,1);
x = x0;
t = 1;

funEvalMultiplier = 1;

% evaluate initial point
[f,g] = funObj(x,varargin{:});
computeHessian = 0;
funEvals = 1;

% output log
if verboseI
    fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
end

% compute optimality of initial point
optCond = max(abs(g));

if nargout > 3
    % Initialize trace
    trace.fval = f;
    trace.funcCount = funEvals;
    trace.optCond = optCond;
end

% exit if initial point is optimal
if optCond <= optTol
    exitflag = 1;
    msg = 'Optimality Condition below optTol';
    if verbose
        fprintf('%s\n',msg)
    end
    if nargout > 3
        output = struct('iterations',0,'funcCount',1,...
            'firstorderopt',max(abs(g)),'message',msg,'trace',trace);
    end
    return;
end

for i=1:maxIter
    %LBFGS
    
    % ******** COMPUTE DESCENT DIRECTION *******************
    
    % update the direction and stepsizes
    if Damped
        if i==1
            d = -g; % initially use steepest descent direction
            old_dirs = zeros(length(g),0);
            old_stps = zeros(length(d),0);
            Hdiag = 1;
        else
            [old_dirs,old_stps,Hdiag] = dampedUpdate(g-g_old,t*d,corrections,...
                debug,old_dirs,old_stps,Hdiag);
            if useMex
                d = lbfgsC(-g,old_dirs,old_stps,Hdiag);
            else
                d = lbfgs(-g,old_dirs,old_stps,Hdiag);
            end
        end
    else
        if i==1
            d = -g; % initially use steepest descent direction
            S = zeros(p,corrections);
            Y = zeros(p,corrections);
            YS = zeros(corrections,1);
            lbfgs_start = 1;
            lbfgs_end = 0;
            Hdiag = 1;
        else
            [S,Y,YS,lbfgs_start,lbfgs_end,Hdiag,skipped] = lbfgsAdd(g-g_old,...
                t*d,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag,useMex);
            if debug && skipped
                fprintf('Skipped L-BFGS updated\n');
            end
            if useMex
                d = lbfgsProdC(g,S,Y,YS,int32(lbfgs_start),int32(lbfgs_end),Hdiag);
            else
                d = lbfgsProd(g,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag);
            end
        end
    end
    g_old = g;
    
    if ~isLegal(d)
        fprintf('Step direction is illegal!\n');
        pause;
        return
    end
    
    
    % *************** COMPUTE STEP LENGTH *******************
    
    % Directional derivative
    gtd = g'*d;
    
    % Check that progress can be made along direction
    if gtd > -progTol
        exitflag=2;
        msg = 'Directional Derivative below progTol';
        break;
    end
    
    % Select Initial Guess
    if i==1
       t = min(1,1/sum(abs(g)));
    else
        if LS_init == 0
            % Newton step
            t = 1;
        elseif LS_init == 1
            % Close to previous step length
            t = t*min(2,(gtd_old)/(gtd));
        elseif LS_init == 2
            % Quadratic Initialization based on {f,g} and previous f
            t = min(1,2*(f-f_old)/(gtd));
        elseif LS_init == 3
            % Double previous step length
            t = min(1,t*2);
        elseif LS_init == 4
            % Scaled step length if possible
            if isempty(HvFunc)
                % No user-supplied Hessian-vector function,
                % use automatic differentiation
                dHd = d'*autoHv(d,x,g,0,funObj,varargin{:});
            else
                % Use user-supplid Hessian-vector function
                dHd = d'*HvFunc(d,x,varargin{:});
            end

            funEvals = funEvals + 1;
            if dHd > 0
                t = -gtd/(dHd);
            else
                t = min(1,2*(f-f_old)/(gtd));
            end
        end
        
        if t <= 0
            t = 1;
        end
    end
    f_old = f;
    gtd_old = gtd;
    
    % Line Search
    [t,f,g,LSfunEvals] = WolfeLineSearch(x,t,d,f,g,gtd,c1,c2,LS_interp,LS_multi,25,progTol,0,0,1,funObj,varargin{:});
    funEvals = funEvals + LSfunEvals;
    %fprintf('step size: %d\n',t);
    x = x + t*d;
    
    % Compute optimality condition
    optCond = max(abs(g));
    
    % Output iteration information
    if verboseI
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,t,f,optCond);
    end
    
    if nargout > 3
        % Update Trace
        trace.fval(end+1,1) = f;
        trace.funcCount(end+1,1) = funEvals;
        trace.optCond(end+1,1) = optCond;
    end
    
    % Check Optimality condition
    if optCond <= optTol
        exitflag = 1;
        msg = 'Optimality Condition below optTol';
        break;
    end
    
    % *************** CHECK FOR LACK OF PROGRESS ******************
    if max(abs(t*d)) <= progTol
        %max(abs(t*d))
        exitflag = 2;
        msg = 'Step size below progTol';
        break
    end
    
    if abs(f-f_old) < progTol
        exitflag = 2;
        msg = 'Function value changing by less than progTol';
        break;
    end
    
    % *************** CHECK FOR GOING OVER ITERATION/EVALUATION LIMIT ****
    if funEvals*funEvalMultiplier >= maxFunEvals
        exitflag = 0;
        msg = 'Reached Maximum Number of Function Evaluation';
        break
    end
    
    if i == maxIter
        exitflag = 0;
        msg='Reached Maximum Number of Iterations';
        break;
    end
    
end

if verbose
    fprintf('%s\n',msg);
end
if nargout > 3
    output = struct('iterations',i,'funcCount',funEvals*funEvalMultiplier,...
        'firstorderopt',max(abs(g)),'message',msg,'trace',trace);
end

end

