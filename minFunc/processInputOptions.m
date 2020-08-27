function [verbose,verboseI,debug,optTol,progTol,maxFunEvals,maxIter,Damped,...
    corrections,useMex,HvFunc,LS_init,LS_interp,LS_multi,c1,c2] = processInputOptions(o)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

verbose  = 1;
verboseI = 1;
debug = 0;

o = toUpper(o);

LS_init = 0;
LS_interp = 2;
LS_multi = 0;
Damped = 0;
c2 = 0.9;

if isfield(o,'DISPLAY')
    switch(upper(o.DISPLAY))
        case 0
            verbose = 0;
            verboseI = 0;
        case 'FINAL'
            verboseI = 0;
        case 'OFF'
            verbose = 0;
            verboseI = 0;
    end
end

maxFunEvals = getOpt(o,'MAXFUNEVALS',1000);
maxIter = getOpt(o,'MAXITER',500);
optTol = getOpt(o,'OPTTOL',1e-5);
progTol = getOpt(o,'PROGTOL',1e-9);
corrections = getOpt(o,'CORRECTIONS',100);
corrections = getOpt(o,'CORR',corrections);
useMex = getOpt(o,'USEMEX',1);
HvFunc = getOpt(o,'HVFUNC',[]);
c1 = getOpt(o,'C1',1e-4);
c2 = getOpt(o,'C2',c2);

LS_init = getOpt(o,'LS_INIT',LS_init);
Damped = getOpt(o,'DAMPED',Damped);
LS_interp = getOpt(o,'LS_interp',LS_interp);
LS_multi = getOpt(o,'LS_multi',LS_multi);

end

function [v] = getOpt(options,opt,default)
if isfield(options,opt)
    if ~isempty(getfield(options,opt))
        v = getfield(options,opt);
    else
        v = default;
    end
else
    v = default;
end
end

function [o] = toUpper(o)
if ~isempty(o)
    fn = fieldnames(o);
    for i = 1:length(fn)
        o = setfield(o,upper(fn{i}),getfield(o,fn{i}));
    end
end
end
